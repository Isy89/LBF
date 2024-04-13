import datetime
import hashlib
import json
import logging
import pathlib
from itertools import combinations
from time import sleep

import click
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import scipy.stats
import seaborn
import statsmodels
import statsmodels.stats.multitest
import statsmodels.stats.weightstats
import stringdb
from PyComplexHeatmap import ClusterMapPlotter, HeatmapAnnotation, anno_simple, anno_scatterplot
from scikit_posthocs import posthoc_dunn
from scipy.stats import kruskal
from sklearn.preprocessing import MinMaxScaler, RobustScaler

import lbfextract.fextract
from lbfextract.transcription_factor_analysis.loaders import ResultsLoader
from lbfextract.transcription_factor_analysis.utils import remove_outliers, generate_time_stamp

logger = logging.getLogger(__name__)
matplotlib.use('Agg')


class DiffSignalAnalysis:
    def __init__(self,
                 df: pl.DataFrame | pd.DataFrame,
                 output_path: pathlib.Path,
                 signal_col_index: list[str],
                 outer_group_column: str,
                 inner_group_column: str,
                 value_column: str,
                 correction_method: str,
                 max_iter: int,
                 alpha: float,
                 overwrite: bool,
                 remove_outliers: bool,
                 save_indivitual_plots: bool, ):
        self.df = pl.from_pandas(df) if isinstance(df, pd.DataFrame) else df
        self.outer_group_column = outer_group_column
        self.inner_group_column = inner_group_column
        self.value_column = value_column
        self.correction_method = correction_method
        self.max_iter = max_iter
        self.alpha = alpha
        self.signal_col_index = signal_col_index
        self.results: pd.DataFrame | None = None
        self.output_path = output_path
        self.overwrite = overwrite
        self.time_stamp = generate_time_stamp()
        self.enrichment_results: dict[str, pd.DataFrame] | None = None
        self.remove_outliers = remove_outliers
        self.save_indivitual_plots = save_indivitual_plots
        self.data_for_plot = None

        if not self.output_path.exists():
            logger.info(
                f"output path: {self.output_path} does not exists. creating it with all missing parent folders.")
            self.output_path.mkdir(exist_ok=True, parents=True)

        self.run_id = self._get_run_id()

        if self.output_path.joinpath(self.get_output_folder()).exists():
            logger.info(f"run id: {self.run_id} already exists.")
            if self.overwrite:
                logger.info(f"overwriting run id: {self.run_id}")
            else:
                logger.info(f"not overwriting run id: {self.run_id}, to overwrite set overwrite=True")
                return
        else:
            self.output_path.joinpath(self.get_output_folder()).mkdir(exist_ok=True, parents=True)

    def _get_run_id(self):
        dict_for_hash = {
            "df": self.df,
            "outer_group_column": self.outer_group_column,
            "inner_group_column": self.inner_group_column,
            "value_column": self.value_column,
            "correction_method": self.correction_method,
            "max_iter": self.max_iter,
            "alpha": self.alpha,
            "signal_col_index": self.signal_col_index,
            "output_path": self.output_path,
        }
        h = hashlib.sha1()

        h.update(json.dumps(dict_for_hash, sort_keys=True, default=str).encode())
        return h.hexdigest()

    def plot_heatmap(self,
                     data,
                     clusters: dict,
                     split_key,
                     save: bool,
                     adjusted_pval: dict[pd.Series],
                     log2_fc: dict[pd.Series],
                     figsize=(10, 10),
                     filename: str = "heatmap_diff_active_signals.pdf",
                     output_path: pathlib.Path = None,
                     col_cluster=True,
                     row_cluster=True,
                     col_split_gap=0.8,
                     row_split_gap=0.8,
                     label='values',
                     row_dendrogram=True,
                     show_rownames=True,
                     show_colnames=True,
                     subplot_gap=5,
                     legend_vpad=15,
                     tree_kws=None,
                     verbose=0,
                     legend_gap=5,
                     legend_hpad=5,
                     cmap='RdYlBu_r',
                     xticklabels_kws={'labelrotation': -90, 'labelcolor': 'black'},
                     **kwargs
                     ):
        tree_kws = {'row_cmap': 'Set1'} if tree_kws is None else tree_kws
        self.data_for_plot = dict(
            data=data,
            clusters=clusters,
            split_key=split_key,
            save=save,
            adjusted_pval=adjusted_pval,
            log2_fc=log2_fc,
            figsize=figsize,
            filename=filename,
            output_path=output_path,
            col_cluster=col_cluster,
            row_cluster=row_cluster,
            col_split_gap=col_split_gap,
            row_split_gap=row_split_gap,
            label=label,
            row_dendrogram=row_dendrogram,
            show_rownames=show_rownames,
            show_colnames=show_colnames,
            subplot_gap=subplot_gap,
            legend_vpad=legend_vpad,
            tree_kws=tree_kws,
            verbose=verbose,
            legend_gap=legend_gap,
            legend_hpad=legend_hpad,
            cmap=cmap,
            xticklabels_kws=xticklabels_kws,
            **kwargs

        )

        output_path = output_path if output_path else pathlib.Path.cwd()

        fig = plt.figure(figsize=figsize)
        ha = HeatmapAnnotation(
            **{k: anno_simple(v, cmap='Set2') for k, v in clusters.items()},
            label_side="bottom",
            axis=0
        )

        col_annot = HeatmapAnnotation(
            **{f"adjp_{k}": anno_scatterplot(adjusted_pval[k], height=20) for k in adjusted_pval},
            **{f"log2_fc_{k}": anno_scatterplot(log2_fc[k], height=20) for k in log2_fc},
            hgap=3,
            legend=False
        )
        row_cluster = row_cluster if data.shape[1] > clusters[split_key].unique().shape[0] else False

        cm1 = ClusterMapPlotter(data,
                                left_annotation=ha,
                                top_annotation=col_annot,
                                col_cluster=col_cluster,
                                row_cluster=row_cluster,
                                row_split=clusters[split_key],
                                col_split_gap=col_split_gap,
                                row_split_gap=row_split_gap,
                                label=label,
                                row_dendrogram=row_dendrogram,
                                show_rownames=show_rownames,
                                show_colnames=show_colnames,
                                subplot_gap=subplot_gap,
                                tree_kws=tree_kws,
                                verbose=verbose,
                                legend_gap=legend_gap,
                                legend_hpad=legend_hpad,
                                legend_vpad=legend_vpad,
                                cmap=cmap,
                                xticklabels_kws=xticklabels_kws,
                                row_cluster_method="average",
                                **kwargs)
        if save:
            plt.savefig(output_path / filename, bbox_inches='tight', dpi=300)
        else:
            fig.show()

        return fig

    def preprocess_df(self):
        df = self.df.to_pandas()
        df = df[~df[self.outer_group_column].isin(["NA", "NaN", "None"])]
        df = df[~df[self.outer_group_column].isna()]
        self.df = df

    def extract_column_for_statistical_test(self, df):
        return df[[self.outer_group_column, self.inner_group_column, self.value_column]]

    def extract_column_for_plot(self, df):
        return df[[self.outer_group_column, self.inner_group_column, self.value_column, "sample_name"]]

    def normalize_between_samples(self, df, normalization_type=RobustScaler, normalization_type_kwgs=None):
        df = df.copy()
        normalization_type_kwgs = {} if normalization_type_kwgs is None else normalization_type_kwgs
        df_group_col = df[[self.outer_group_column, "sample_name"]].copy()
        df = df.pivot_table(index="sample_name",
                            columns=self.inner_group_column,
                            values=self.value_column)
        columns = df.columns.copy()
        index = df.index.copy()
        df = normalization_type(**normalization_type_kwgs).fit_transform(df)
        df = pd.DataFrame(df, columns=columns, index=index)
        df = df.reset_index().melt(id_vars=["sample_name"], var_name=self.inner_group_column,
                                   value_name=self.value_column)
        return pd.merge(df, df_group_col, on='sample_name', how='inner')

    @staticmethod
    def attach_group_column(x, group_name):
        return pd.DataFrame(
            {
                "group": x.apply(lambda x: group_name),
                "value": x
            }
        )

    def get_diff_active_gi_two_groups(self, genomic_intervals_groups, groups_names_combinations):
        group_1_name, group_2_name = groups_names_combinations[0]
        matrix_res = []
        for gi in genomic_intervals_groups.groups:
            gi_groups_gruped_by_group = genomic_intervals_groups.get_group(gi).groupby(self.outer_group_column)

            long_df = pd.concat([
                self.attach_group_column(remove_outliers(group[self.value_column]), name)
                if self.remove_outliers
                else self.attach_group_column(group[self.value_column], name)
                for name, group
                in gi_groups_gruped_by_group])

            groups_values = {name: group.value.values for name, group in long_df.groupby("group")}
            _, p_value = scipy.stats.mannwhitneyu(groups_values[group_1_name], groups_values[group_2_name])

            group_df_means = long_df.groupby("group").mean()
            group_df_stds = long_df.groupby("group").std()

            for comb in groups_names_combinations:
                matrix_res.append(
                    [
                        "-".join(comb),
                        group_df_means.loc[comb[0], "value"],
                        group_df_means.loc[comb[1], "value"],
                        group_df_stds.loc[comb[0], "value"],
                        group_df_stds.loc[comb[1], "value"],
                        gi,
                        p_value
                    ]
                )

        columns_names = [
            "group_pairs",
            "mean_group_1",
            "mean_group_2",
            "std_group_1",
            "std_group_2",
            self.inner_group_column,
            "p_value"
        ]

        df = pd.DataFrame(matrix_res, columns=columns_names)
        df[self.inner_group_column] = df[self.inner_group_column].str.replace("/", "_")
        return df

    def get_diff_active_gi_multiple_groups(self, genomic_intervals_groups, groups_names_combinations):
        matrix_res = []
        for gi in genomic_intervals_groups.groups:
            gi_groups_gruped_by_group = genomic_intervals_groups.get_group(gi).groupby(self.outer_group_column)

            long_df = pd.concat([
                self.attach_group_column(remove_outliers(group[self.value_column]), name)
                if self.remove_outliers
                else self.attach_group_column(group[self.value_column], name)
                for name, group
                in gi_groups_gruped_by_group])

            kruskal_p_value = kruskal(*[group.value for name, group in long_df.groupby("group")]).pvalue
            perform_dunn_test = True if kruskal_p_value <= self.alpha else False
            if perform_dunn_test:
                dunn_test_res = posthoc_dunn(long_df, val_col="value", group_col="group")

            group_df_means = long_df.groupby("group").mean()
            group_df_stds = long_df.groupby("group").std()

            for comb in groups_names_combinations:

                if perform_dunn_test:
                    p_value = dunn_test_res.loc[comb[0], comb[1]]
                else:
                    p_value = np.nan

                matrix_res.append(
                    [
                        "-".join(comb),
                        group_df_means.loc[comb[0], "value"],
                        group_df_means.loc[comb[1], "value"],
                        group_df_stds.loc[comb[0], "value"],
                        group_df_stds.loc[comb[1], "value"],
                        gi,
                        kruskal_p_value,
                        p_value
                    ]
                )

        columns_names = [
            "group_pairs",
            "mean_group_1",
            "mean_group_2",
            "std_group_1",
            "std_group_2",
            self.inner_group_column,
            "kruskal_p_value",
            "p_value"
        ]

        df = pd.DataFrame(matrix_res, columns=columns_names)
        df[self.inner_group_column] = df[self.inner_group_column].str.replace("/", "_")
        return df

    def run(self):
        self.preprocess_df()
        df_for_statistical_tests = self.extract_column_for_plot(self.df)

        # df_for_statistical_tests = self.normalize_between_samples(df_for_statistical_tests)
        genomic_intervals_groups = df_for_statistical_tests.groupby(self.inner_group_column)
        group_names = self.df[self.outer_group_column].unique()
        groups_names_combinations = sorted(combinations(group_names, 2))
        n_groups = len(group_names)

        if n_groups == 2:
            df = self.get_diff_active_gi_two_groups(genomic_intervals_groups, groups_names_combinations)
            self.results = self.get_global_correction_two_groups(df)
        elif n_groups > 2:
            df = self.get_diff_active_gi_multiple_groups(genomic_intervals_groups, groups_names_combinations)
            self.results = self.get_global_correction_multiple_groups(df)
        else:
            raise ValueError(f"Number of groups should be 2 or more found {len(group_names)}")

        diff_active_region_per_group_pair = self.results.query("Rejected == True").groupby("group_pairs")

        self.enrichment_results = {
            f"{group_pair}": self._get_enrichment(
                diff_active_region_per_group_pair.get_group(group_pair)[self.inner_group_column]
                    .to_list()
            )
            for group_pair in diff_active_region_per_group_pair.groups
        }
        return self

    def get_global_correction_multiple_groups(self, df_):
        df_ = df_.copy()

        unique_kruskal = df_.loc[:, [self.inner_group_column, "kruskal_p_value"]].drop_duplicates()

        unique_kruskal.loc[:, "kruskal_p_value"] = statsmodels.stats.multitest.multipletests(
            unique_kruskal["kruskal_p_value"],
            alpha=self.alpha,
            method=self.correction_method,
            maxiter=self.max_iter,
            is_sorted=False,
            returnsorted=False
        )[1]
        unique_kruskal = unique_kruskal.set_index(self.inner_group_column)

        df_.loc[:, "kruskal_p_value"] = df_[self.inner_group_column].apply(
            lambda x: unique_kruskal.loc[x, "kruskal_p_value"]
        )

        passed_kruscal = df_.loc[:, "kruskal_p_value"] <= self.alpha

        if passed_kruscal.sum() == 0:
            df_.loc[:, "Rejected"] = False
            df_.loc[:, "adj_p-val"] = np.nan
            df_.loc[:, "p_value"] = np.nan
            return df_

        df_.loc[passed_kruscal, "adj_p-val"] = statsmodels.stats.multitest.multipletests(
            df_.loc[passed_kruscal, "p_value"],
            alpha=self.alpha,
            method=self.correction_method,
            maxiter=self.max_iter,
            is_sorted=False,
            returnsorted=False
        )[1]
        df_.loc[df_.loc[:, "adj_p-val"].isna(), "p_value"] = np.nan
        df_.loc[:, "Rejected"] = np.logical_and(
            df_.loc[:, "adj_p-val"] <= self.alpha,
            passed_kruscal
        )
        return df_

    def get_global_correction_two_groups(self, df_):
        df_ = df_.copy()

        df_.loc[:, "adj_p-val"] = statsmodels.stats.multitest.multipletests(
            df_.loc[:, "p_value"],
            alpha=self.alpha,
            method=self.correction_method,
            maxiter=self.max_iter,
            is_sorted=False,
            returnsorted=False
        )[1]
        df_.loc[:, "Rejected"] = df_.loc[:, "adj_p-val"] <= self.alpha
        return df_

    def plot_barplot_rejected_per_group(self):
        fig, ax = plt.subplots(1, figsize=(10, 10), dpi=300)
        df_to_plot = self.results.query("Rejected == True").groupby("group_pairs").count()[["Rejected"]]
        if not df_to_plot.empty:
            df_to_plot.plot(
                kind="bar", ax=ax)
            counts = self.results.query("Rejected == True").groupby("group_pairs").count()["Rejected"]
            max_counts = counts.max()
            for i, count in enumerate(counts):
                ax.text(i, count + (max_counts * 0.03), str(count), ha='center',
                        bbox=dict(boxstyle="round", fc="w", ec="k"))
            ax.spines[['right', 'top', 'bottom']].set_visible(False)
            return fig
        else:
            logger.warning("No group pairs were rejected")
            return None

    def plot_coverage(self, gi):
        signal_shape = len(self.signal_col_index)
        if self.results is None:
            logger.warning("Analysis has not been run yet. Running analysis first")
            self.run()

        df = self.df.loc[
            self.df[self.inner_group_column] == gi,
            [self.outer_group_column, self.inner_group_column] + self.signal_col_index
        ]
        df = df[~df[self.outer_group_column].isin(["NA", "NaN", "None"])].set_index(self.outer_group_column).drop(
            "genomic_interval", axis=1)
        colors_std, colors_means = {}, {}
        outer_groups = self.df[self.outer_group_column].unique()
        palette = list(seaborn.color_palette("Paired", len(outer_groups) * 2))
        for group in sorted(outer_groups):
            colors_std[group] = palette.pop()
            colors_means[group] = palette.pop()
        means = df.groupby(self.outer_group_column).mean()
        stds = df.groupby(self.outer_group_column).std()

        fig, ax = plt.subplots(figsize=(20, 10))
        fig.suptitle(gi, fontsize=20)

        for i in stds.index.to_list():
            ax.plot(means.loc[i], color=colors_means[i], label=i)
            ax.plot(means.loc[i] + stds.loc[i], linestyle='dashed', color="lightgray")
            ax.plot(means.loc[i] - stds.loc[i], linestyle='dashed', color="lightgray")
            ax.fill_between(np.arange(signal_shape), means.loc[i] - stds.loc[i], means.loc[i] + stds.loc[i],
                            color=colors_std[i],
                            alpha=0.3)
        ax.spines[['right', 'top', 'bottom', "left"]].set_visible(False)
        ax.set_xticks(list(range(0, signal_shape + 1, signal_shape // 4)),
                      list(range(0, signal_shape + 1, signal_shape // 4)))
        ax.axvline(signal_shape // 2, ls="-.", color="black",
                   label="GI center")
        fig.legend()
        ax.grid(linestyle='dashed')
        ax.set_ylabel("normalized coverage")
        ax.set_xlabel("position aruond GI center")
        return fig

    def get_output_folder(self):
        return f"{self.time_stamp}__{self.run_id}"

    def save_results_per_group_pair(self, any_rejected: bool):

        if self.results is None:
            logger.warning("Analysis has not been run yet. Running analysis first")
            self.run()
        run_folder = self.get_output_folder()
        res_by_group_pair = self.results.groupby("group_pairs")
        for group_pair in res_by_group_pair.groups:
            output_path_group_pair = self.output_path / run_folder / group_pair
            output_path_group_pair.mkdir(exist_ok=True, parents=True)
            group = res_by_group_pair.get_group(group_pair)
            (
                group
                    .sort_values("adj_p-val", ascending=True)
                    .to_csv(output_path_group_pair / f"fextract_diff_signal_analysis_{group_pair}.csv")
            )
            if self.enrichment_results.get(group_pair, None) is not None:
                self.enrichment_results.get(group_pair).to_csv(output_path_group_pair / f"enrichment_{group_pair}.csv")

            if any_rejected and self.save_indivitual_plots:
                rejected_gi = group.query('Rejected == True')
                for genomic_interval in rejected_gi[self.inner_group_column]:
                    fig = self.plot_coverage(gi=genomic_interval)
                    fig.savefig(output_path_group_pair / f"coverage_{genomic_interval}.png")
                    plt.close(fig)
        return self

    def save_results(self):
        output_path = self.output_path / self.get_output_folder()
        output_path.mkdir(exist_ok=True, parents=True)

        rejected_tfs = self.results.query("Rejected == True").genomic_interval.unique()

        df_heatmap = (
            self.extract_column_for_plot(self.df.loc[self.df[self.inner_group_column].isin(rejected_tfs)])
                .pivot_table(index="sample_name",
                             columns=self.inner_group_column,
                             values=self.value_column,
                             fill_value=0)
        )

        def get_outliars_mask(column, bound="lower"):
            Q1 = column.quantile(0.25)
            Q3 = column.quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            if bound == "lower":
                return column < lower_bound
            else:
                return column > upper_bound

        def replace_outliars(x_):
            x = x_.copy()
            mask_lower = get_outliars_mask(x, bound="lower")
            mask_upper = get_outliars_mask(x, bound="upper")
            x[mask_lower] = x[~mask_lower].min()
            x[mask_upper] = x[~mask_upper].max()

            return x

        any_rejected = True if rejected_tfs.shape[0] > 0 else False

        clusters = {"groups": self.df.groupby("sample_name")[self.outer_group_column].first()}
        filter_rows = self.results[self.inner_group_column].isin(rejected_tfs)
        adjusted_pvals = (
            self.results.loc[filter_rows]
                .pivot_table(index="group_pairs",
                             columns=self.inner_group_column,
                             values="adj_p-val",
                             fill_value=np.nan).T
        )
        mean_1 = self.results.mean_group_1.values
        mean_2 = self.results.mean_group_2.values
        fc = mean_1 / mean_2
        log2_fc = np.sign(fc) * np.where(np.logical_and(fc > 0, fc < 1), np.log2(fc), np.abs(np.log2(np.abs(fc))))
        self.results.loc[:, "log2_fc"] = log2_fc

        log2_fc_vals = (
            self.results
                .assign(log2_fc=log2_fc)
                .pivot_table(index="group_pairs",
                             columns=self.inner_group_column,
                             values="log2_fc",
                             fill_value=np.nan)
        ).T

        if not df_heatmap.empty:
            df_heatmap = pd.DataFrame(
                data=MinMaxScaler().fit_transform(df_heatmap.apply(lambda x: replace_outliars(x), axis=0)),
                columns=df_heatmap.columns,
                index=df_heatmap.index
            )

            figsize_height = 10 + (df_heatmap.shape[1] // 20) + (clusters["groups"].unique().shape[0] * 2)
            figsize_width = 10 + ((df_heatmap.shape[1] // 20) * 3)
            self.plot_heatmap(
                df_heatmap,
                clusters=clusters,
                split_key="groups",
                adjusted_pval=adjusted_pvals,
                log2_fc=log2_fc_vals.loc[rejected_tfs, :],
                output_path=output_path,
                filename="fextract_diff_signals_heatmap.png",
                save=True,
                figsize=(
                    figsize_width,
                    figsize_height,
                ),
                cmap="RdYlBu_r"
            )
        fig = self.plot_barplot_rejected_per_group()
        if fig:
            fig.savefig(output_path / "fextract_diff_signals_per_group_bar_plot.png")
            plt.close(fig)
        self.save_results_per_group_pair(any_rejected=any_rejected)

        return self

    def _get_enrichment(self, genes: list[str]) -> pd.DataFrame | None:
        if self.results is not None:
            if not self.results.empty:
                try:
                    string_ids = stringdb.get_string_ids(genes).queryItem.to_list()
                    sleep(1)
                    return stringdb.get_enrichment(string_ids, )[
                        ["category", "preferredNames", "p_value", "fdr", "description"]]
                except ValueError as e:
                    logger.warning("There was an error with the string database server. "
                                   "Enrichment analysis could not be retrieved. the error was:\n"
                                   f"{e}")
                except KeyError as e:
                    logger.warning("Get enrichment failed due to a databse problem")
                    return None

            else:
                logger.warning("results are empty, enrichment can not be calculated")
                return None
        else:
            logger.warning("Analysis has not been run yet. Running analysis first")
            self.run()
            return self._get_enrichment(genes)


def parse_indices(ctx, param, value: str) -> tuple[int, int] | None:
    return tuple([int(i) for i in value.split(",")])[0:2] if value else None


def create_path(ctx, param, value):
    if not pathlib.Path(value).exists():
        pathlib.Path(value).mkdir(parents=True, exist_ok=True)
    return value


@click.group(name="post_extraction_analysis_commands")
@click.option("--path_to_res_summary",
              required=True,
              type=click.Path(exists=True,
                              path_type=pathlib.Path,
                              dir_okay=True,
                              file_okay=False,
                              readable=True,
                              resolve_path=False,
                              allow_dash=True),
              help="Path to the directory containing the fextract in batch results."
              )
@click.option("--signal_length", type=int, default=4000, show_default=True)
@click.option("--center_signal_indices", type=str, default="1800,2200", callback=parse_indices, show_default=True)
@click.option("--flanking_signal_indices", type=str, default="1000,3000", callback=parse_indices, show_default=True)
@click.option("--normalize", is_flag=True, default=False, show_default=True)
@click.option("--output_path",
              type=click.Path(
                  exists=False, path_type=pathlib.Path, dir_okay=True, file_okay=False,
                  readable=True, resolve_path=False, allow_dash=True, writable=True
              ),
              default=pathlib.Path.cwd(), show_default=True, callback=create_path)
@click.option("--overwrite", is_flag=True, default=False, show_default=True)
@click.option("--path_to_sample_sheet", type=click.Path(exists=True, path_type=pathlib.Path, dir_okay=False,
                                                        file_okay=True, readable=True, resolve_path=False,
                                                        allow_dash=True), default=None, show_default=True)
@click.option("--outer_group_column", type=str, default="group", show_default=True)
@click.option("--inner_group_column", type=str, default="genomic_interval", show_default=True)
@click.option("--value_column", type=str, default="amplitude", show_default=True)
@click.option("--correction_method", type=str, default="hs", show_default=True)
@click.option("--max_iter", type=int, default=1, show_default=True)
@click.option("--alpha", type=float, default=0.05, show_default=True)
@click.option("--remove_outliars", is_flag=True, default=False, show_default=True)
@click.option("--commit_hash", type=str, default="Not Provided", show_default=True)
@click.option("--save_indivitual_plots", is_flag=True, default=False, show_default=True)
@click.pass_context
def cli(ctx,
        path_to_res_summary: pathlib.Path,
        signal_length: int,
        center_signal_indices: tuple,
        flanking_signal_indices: tuple,
        normalize: bool,
        output_path: pathlib.Path,
        outer_group_column: str,
        inner_group_column: str,
        value_column: str,
        correction_method: str,
        max_iter: int,
        alpha: float,
        overwrite: bool,
        path_to_sample_sheet: pathlib.Path,
        remove_outliars: bool,
        commit_hash: str,
        save_indivitual_plots: bool):
    ctx.obj = dict()
    ctx.obj["path_to_res_summary"] = path_to_res_summary
    ctx.obj["signal_length"] = signal_length
    ctx.obj["center_signal_indices"] = center_signal_indices
    ctx.obj["flanking_signal_indices"] = flanking_signal_indices
    ctx.obj["normalize"] = normalize
    ctx.obj["output_path"] = output_path
    ctx.obj["outer_group_column"] = outer_group_column
    ctx.obj["inner_group_column"] = inner_group_column
    ctx.obj["value_column"] = value_column
    ctx.obj["correction_method"] = correction_method
    ctx.obj["max_iter"] = max_iter
    ctx.obj["alpha"] = alpha
    ctx.obj["overwrite"] = overwrite
    ctx.obj["path_to_sample_sheet"] = path_to_sample_sheet
    ctx.obj["remove_outliers"] = remove_outliars
    ctx.obj["commit_hash"] = commit_hash
    ctx.obj["save_indivitual_plots"] = save_indivitual_plots


@cli.command()
@click.pass_context
def get_differentially_active_genomic_intervals(ctx):
    starting_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    output_path = ctx.obj["output_path"] / "fextract_diff_signal_results"
    output_path.mkdir(exist_ok=True)
    accessibility_extraction_config = dict(
        start=ctx.obj["center_signal_indices"][0],
        end=ctx.obj["center_signal_indices"][1]
    )
    df_results = ResultsLoader(ctx.obj["path_to_res_summary"],
                               accessibility_extraction_config,
                               ctx.obj["signal_length"],
                               ctx.obj["flanking_signal_indices"],
                               ctx.obj["normalize"], path_to_sample_sheet=ctx.obj["path_to_sample_sheet"],
                               grouping_column=ctx.obj["outer_group_column"]).load()

    diff_signal_analysis = DiffSignalAnalysis(df=df_results,
                                              outer_group_column=ctx.obj["outer_group_column"],
                                              inner_group_column=ctx.obj["inner_group_column"],
                                              value_column=ctx.obj["value_column"],
                                              correction_method=ctx.obj["correction_method"],
                                              max_iter=ctx.obj["max_iter"],
                                              alpha=ctx.obj["alpha"],
                                              signal_col_index=[str(i) for i in range(ctx.obj["signal_length"])],
                                              output_path=output_path,
                                              overwrite=ctx.obj["overwrite"],
                                              remove_outliers=ctx.obj["remove_outliers"],
                                              save_indivitual_plots=ctx.obj["save_indivitual_plots"])

    diff_signal_analysis.run()
    diff_signal_analysis.save_results()

    metadata = dict(
        paccakge_repo="LBFextract",
        lbf_version=f"{lbfextract.__version__}",
        commithash=ctx.obj["commit_hash"],
        starting_time=starting_time,
        finisheing_time=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        run_id=diff_signal_analysis.run_id,
        path_to_res_summary=ctx.obj["path_to_res_summary"],
        signal_length=ctx.obj["signal_length"],
        center_signal_indices=ctx.obj["center_signal_indices"],
        flanking_signal_indices=ctx.obj["flanking_signal_indices"],
        normalize=ctx.obj["normalize"],
        output_path=output_path,
        outer_group_column=ctx.obj["outer_group_column"],
        inner_group_column=ctx.obj["inner_group_column"],
        value_column=ctx.obj["value_column"],
        correction_method=ctx.obj["correction_method"],
        max_iter=ctx.obj["max_iter"],
        alpha=ctx.obj["alpha"],
        overwrite=ctx.obj["overwrite"],
        remove_outliers=ctx.obj["remove_outliers"],
        save_indivitual_plots=ctx.obj["save_indivitual_plots"]
    )

    results = dict(
        results=diff_signal_analysis.results,
        enrichment_results=diff_signal_analysis.enrichment_results,
        data_for_plot=diff_signal_analysis.data_for_plot

    )
    import dill
    with open(output_path / diff_signal_analysis.get_output_folder() / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=4, sort_keys=True, default=str)
    with open(output_path / diff_signal_analysis.get_output_folder() / "results.pkl", "wb") as f:
        dill.dump(results, f)
    logger.info("finished")


@cli.command()
@click.pass_context
def generate_sample_sheet(ctx):
    output_path = ctx.obj["output_path"] / "fextract_diff_signal_results"
    output_path.mkdir(exist_ok=True)
    df_results = ResultsLoader(ctx.obj["path_to_res_summary"], ctx.obj["signal_length"],
                               ctx.obj["center_signal_indices"], ctx.obj["flanking_signal_indices"],
                               ctx.obj["normalize"], path_to_sample_sheet=ctx.obj["path_to_sample_sheet"])
    sample_sheet = df_results.generate_sample_sheet()
    sample_sheet.to_csv(output_path / "sample_sheet.csv", index=True)


if __name__ == "__main__":
    cli(obj={})
