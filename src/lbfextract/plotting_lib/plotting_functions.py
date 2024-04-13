import itertools
import logging
from typing import Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)


def plot_signal(summary: np.array,
                apply_savgol: bool = True,
                savgol_window_length: int = 33,
                savgol_polyorder: int = 3,
                ax: matplotlib.pyplot.Axes = None,
                fig: matplotlib.pyplot.Figure = None,
                color: Union[int, str] = None,
                label: str = None,
                label_center_line: str = None,
                title=None,
                title_font_size=20,
                general_font_size=15,
                plot_center_line_label=True,
                line_type="-"):
    """
    function to plot the signal summarised according to a specific function

    :param summary: signal per position
    :param apply_savgol: whether the smoothing need to be applied
    :param savgol_window_length: length of the window used to smooth
    :param savgol_polyorder: order of polynomial used by the savgol filter
    :param ax: matplotlib axis of the plot
    :param fig: matplotlib Figure
    :param color:
    :param label:
    :return:
    """

    if ax is None:
        fig, ax = plt.subplots(1, figsize=(10, 10))

    if title:
        fig.suptitle(title, fontsize=title_font_size)

    std = None

    if len(summary.shape) > 1:
        summary = np.nanmean(summary, axis=0)
        std = summary.std(axis=0)

    if apply_savgol:
        summary = savgol_filter(summary, savgol_window_length, savgol_polyorder)
        std = savgol_filter(std, savgol_window_length,
                            savgol_polyorder) if std is not None else None

    ax.plot(summary, label=f"{label}" if label else "signal summary", c=color, linestyle=line_type)
    if std is not None:
        color = "lightgray"
        ax.fill_between(np.arange(summary.shape[0]), summary - std, summary + std, alpha=0.5, color=color)

    label_center_line = label_center_line if label_center_line else "TFBS center" if plot_center_line_label else None

    ax.set_xticks([0, summary.shape[0] // 2, summary.shape[0]],
                  labels=[-(summary.shape[0] // 2), 0, summary.shape[0] // 2])
    ax.axvline(summary.shape[0] // 2, ls="-.", color="red",
               label=label_center_line)
    ax.legend()
    ax.set_ylabel(label, fontsize=general_font_size)
    ax.set_xlabel("positions", fontsize=general_font_size)
    ax.spines[['right', 'top', 'bottom']].set_visible(False)
    return fig, ax

def plot_signal_batch(df: pd.DataFrame,
                      apply_savgol: bool = False,
                      savgol_window_length: int = 11,
                      savgol_polyorder: int = 3,
                      signal="coverage",
                      title="coverage at TFBS",
                      figsize=(20, 10),
                      ax=None,
                      fig=None,
                      color=None,
                      label=None,
                      alpha=1,
                      linewidth=1,
                      flanking=1500,
                      xlabel=None,
                      top=5,
                      bottom=5,
                      window_center=50,
                      mask_rows=None,
                      ):
    """
    function to plot the signal summarised according to a specific function
    Args:
        df: array describing the signal
        apply_savgol: wether or not to apply savgol_filter
        savgol_window_length: window to smoth over
        savgol_polyorder: degree of polynomial to use

    Returns:
        matplolib Figure object representing the plot

    """

    df = df.copy()

    if df.shape[0] > 20:
        if top + bottom > df.shape[0]:
            logger.error(f"top + bottom is too large, it should be less than {df.shape[0]}")
            top = 10

        if flanking * 2 > df.shape[1]:
            logger.error(f"flanking is too large, it should be less than {df.shape[1] // 2}")
            flanking = (df.shape[1] // 5) * 2

        if window_center * 2 > df.shape[1]:
            logger.error(f"window_center is too large, it should be less than {df.shape[1] // 2}")

        mask = np.logical_or(np.arange(df.shape[1]) <= flanking,
                             np.arange(df.shape[1]) > df.shape[1] - flanking)
        center = df.shape[1] // 2

        df["amplitude"] = df.iloc[:, center - window_center: center + window_center].mean(axis=1) / df.iloc[:,
                                                                                                    mask].mean(axis=1)

        if mask_rows is None:
            mask_rows = np.logical_or(np.arange(df.shape[0]) <= bottom, np.arange(df.shape[0]) > df.shape[0] - top)
            df = df.sort_values(by="amplitude", ascending=False).iloc[mask_rows, :-2]
        else:
            df = df.iloc[mask_rows, :-2]

    plus_minus = u"\u00B1"
    if not fig:
        fig, ax = plt.subplots(figsize=figsize)
        fig.suptitle(title, fontsize=20)
        ax.set_xlabel(xlabel if xlabel else "signal per interval", fontsize=15)
        ax.set_ylabel(signal, fontsize=15)
    if apply_savgol:
        index = df.index.copy()
        df = savgol_filter(df, savgol_window_length, savgol_polyorder)
        df = pd.DataFrame(df, index=index)
    ax.plot(df.values.flatten(), label=label if label else '', c=color, alpha=alpha, linewidth=linewidth)
    xtick_labels = ([str(int(-df.shape[1] / 2))]
                    + list(itertools.chain(*[[i[0], i[1]] for i in itertools.zip_longest(df.index.to_list(), [
                plus_minus + str(int(df.shape[1] / 2))] * df.shape[0])]))
                    )
    xtick_labels[-1] = str(int(+df.shape[1] / 2))
    ax.set_xticks([int(i) for i in range(0, df.shape[0] * df.shape[1] + 1, int(df.shape[1] / 2))],
                  labels=xtick_labels
                  )
    colors = seaborn.color_palette(n_colors=df.shape[0])

    for count, i in enumerate(
            [i for i in range(df.shape[1] // 2, df.shape[0] * df.shape[1], df.shape[1])]):
        ax.axvline(i, ls="-.", color=colors[count], alpha=0.5, linewidth=linewidth)
    for count, i in enumerate([i for i in range(0, df.shape[0] * df.shape[1] + 1, df.shape[1])]):
        ax.axvline(i, ls=":", color="gray", alpha=0.5)
    ax.spines[['right', 'top', 'bottom']].set_visible(False)
    ax.annotate('', xy=(0, -0.1), xycoords='axes fraction', xytext=(1, -0.1))

    return fig, ax


def plot_heatmap_signal_batch(fig=None, ax=None, array=None,
                              title=None,
                              title_font_size=20,
                              general_font_size=15):
    if not ax:
        fig, ax = plt.subplots(figsize=(20, 10))

    if title:
        ax.set_title(title, fontsize=title_font_size)

    seaborn.heatmap(array, ax=ax, cmap="viridis", cbar_kws={'label': 'coverage'}, rasterized=True)
    ax.set_xticks([0, array.shape[1] // 2, array.shape[1]],
                  labels=[-(array.shape[1] // 2), "TFBS center", array.shape[1] // 2], rotation=0)
    ax.set_yticks(ax.get_yticks(), ax.get_yticklabels(), rotation=0)
    ax.axvline(array.shape[1] // 2, ls="-.", color="white",
               label="TFBS center")
    ax.set_xlabel("position", fontsize=general_font_size)

    return fig, ax


def plot_heatmap_kde_amplitude(array=None,
                               fig=None,
                               ax=None,
                               title=None,
                               title_font_size=20,
                               general_font_size=15,
                               tf_to_annotate=("CTCF",),
                               ylabel=None,
                               flanking=1500,
                               annotation_center_line="center",
                               window_center=50,
                               top=5,
                               bottom=5,
                               ):
    array = array.copy()
    if not ax:
        fig, ax = plt.subplots(2, 2, figsize=(20, 10))

    if title:
        fig.suptitle(title, fontsize=title_font_size)

    pca = PCA(2)
    trf_array = pd.DataFrame(pca.fit_transform(array))
    trf_array.columns = ["PC1", "PC2"]

    if flanking * 2 > array.shape[1]:
        logger.error(f"flanking is too large, it should be less than {array.shape[1] // 2}")
        flanking = (array.shape[1] // 5) * 2

    if window_center * 2 > array.shape[1]:
        logger.error(f"window_center is too large, it should be less than {array.shape[1] // 2}")

    mask = np.logical_or(np.arange(array.shape[1]) <= flanking,
                         np.arange(array.shape[1]) > array.shape[1] - flanking)
    center = array.shape[1] // 2

    trf_array["amplitude"] = (
            array.iloc[:, center - window_center: center + window_center].mean(axis=1) / array.iloc[:, mask].mean(
        axis=1)
    ).values
    trf_array.index = array.index
    if trf_array.shape[0] > 2:
        # with less than 2 rows, seaborn.kdeplot fails
        seaborn.kdeplot(
            data=trf_array, x="PC1", y="PC2", ax=ax[0, 1], color="lightgray")
    seaborn.scatterplot(data=trf_array, x="PC1", y="PC2", ax=ax[0, 1], sizes="amplitude")

    array["amplitude"] = trf_array.amplitude.copy()
    row_slice = slice(None, None, None)
    if array.shape[0] > 20:

        if top + bottom > array.shape[0]:
            logger.error(f"top + bottom is too large, it should be less than {array.shape[0]}")
            top = 10

        row_slice = np.logical_or(np.arange(array.shape[0]) > (array.shape[0] - top),
                                  np.arange(array.shape[0]) <= bottom)
    df = array.sort_values(by="amplitude", ascending=False)[row_slice].copy()
    tf_to_annotate = tf_to_annotate if tf_to_annotate is not None else df.index.to_list()
    marker_size = 20 if array.shape[0] < 20 else 10 if array.shape[0] < 50 else 5
    for tf in tf_to_annotate:
        if tf in trf_array.index:
            x, y, _ = trf_array.loc[tf].to_list()
            ax[0, 1].plot(x, y, marker="o", markersize=marker_size // 2, markeredgecolor="red", markerfacecolor="red")
            ax[0, 1].annotate(tf, xy=(x, y))

    for (point, amplitude) in zip(np.arange(df.shape[0]), df.amplitude.to_list()):
        ax[0, 0].axhline(point, color="lightgray")
        ax[0, 0].plot(amplitude, point, marker="o", markersize=marker_size, markeredgecolor="blue",
                      markerfacecolor="lightblue")
    ax[0, 0].set_yticks(np.arange(df.shape[0]), labels=df.index.to_list(), rotation=0)
    ax[0, 0].set_xlabel("amplitude", fontsize=general_font_size)

    df.drop(columns=["amplitude"], inplace=True)

    fig, _ = create_heatmap_bed_per_position(df=df, fig=fig, ax=ax[1, 0], title="", ylabel=None,
                                             general_font_size=general_font_size,
                                             title_font_size=title_font_size,
                                             annotation_center_line=annotation_center_line)

    fig, _ = create_masked_corr_matrix_plot(df, fig=fig, ax=ax[1, 1], title=None, title_font_size=title_font_size,
                                            general_font_size=general_font_size,
                                            axlabel=ylabel)
    for i in range(2):
        for l in range(2):
            if i == 0 and l == 0:
                ax[i, l].spines[['right', 'top']].set_visible(False)

            if i == 0 and l == 1:
                ax[i, l].spines[['right', 'top']].set_visible(False)
    fig.tight_layout()

    return fig, ax


def correlation_map_plot(df):
    corr_mat = df.T.corr().stack().reset_index(name="correlation")
    fig = seaborn.relplot(
        data=corr_mat,
        x="level_0", y="level_1", hue="correlation", size="correlation",
        palette="vlag", hue_norm=(-1, 1), edgecolor=".7",
        height=10, sizes=(50, 250), size_norm=(-.2, .8),
    )
    fig.set(xlabel="", ylabel="", aspect="equal")
    fig.despine(left=True, bottom=True)
    fig.ax.margins(.02)
    for label in fig.ax.get_xticklabels():
        label.set_rotation(90)
    for artist in fig.legend.legendHandles:
        artist.set_edgecolor(".7")
    return fig


def create_masked_corr_matrix_plot(df, fig=None, ax=None, title=None, title_font_size=20, general_font_size=15,
                                   axlabel=None, label_color_bar=""):
    df = df.T.copy().corr()
    if not ax or not fig:
        fig, ax = plt.subplots(1, 1, figsize=(20, 10))

    mask = np.triu(np.ones_like(df, dtype=bool))
    annot = True
    if df.shape[0] > 20 or df.shape[1] > 20:
        annot = False
    seaborn.heatmap(df, mask=mask, vmin=-1, vmax=1, annot=annot, cmap='viridis', ax=ax,
                    cbar_kws={'label': label_color_bar, 'orientation': 'horizontal'}, annot_kws={"size": 8})
    if title:
        ax.set_title(title, fontdict={'fontsize': title_font_size}, pad=16)
    if df.shape[0] <= 20 and df.shape[1] <= 20:
        ax.set_xticks(np.arange(df.shape[0]) + 0.5, labels=df.index.to_list(), rotation=90, fontsize=8)
        ax.set_yticks(np.arange(df.shape[0]) + 0.5, labels=df.index.to_list(), rotation=0, fontsize=8)
    else:
        ax.set_yticklabels([])
        ax.set_xticklabels([])
    if axlabel:
        ax.set_ylabel(axlabel, fontsize=general_font_size)
        ax.set_xlabel(axlabel, fontsize=general_font_size)
    return fig, ax


def create_heatmap_bed_per_position(df, fig=None, ax=None, title="", ylabel=None, general_font_size=15,
                                    title_font_size=20, annotation_center_line: str = None,
                                    label_color_bar=""):
    if not ax or not fig:
        fig, ax = plt.subplots(2, 2, figsize=(20, 10))

    seaborn.heatmap(df, ax=ax, cmap="viridis", cbar_kws={'label': label_color_bar, 'orientation': 'horizontal'},
                    rasterized=True, cbar=True)
    if df.shape[0] > 20:
        ax.set_yticks([])
    else:
        ax.set_yticks(np.arange(df.shape[0]) + 0.5, labels=df.index.to_list(), rotation=0, fontsize=8)
    if title:
        ax.set_title(title, fontdict={'fontsize': title_font_size}, pad=16)
    ax.set_xticks([0, df.shape[1] // 2, df.shape[1]],
                  labels=[-(df.shape[1] // 2), annotation_center_line, df.shape[1] // 2], rotation=0,
                  fontsize=general_font_size)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=general_font_size / 2)
    ax.axvline(df.shape[1] // 2, ls="-.", color="white",
               label=annotation_center_line)
    ax.set_xlabel("position", fontsize=general_font_size)
    ax.set_ylabel(ylabel=ylabel, fontsize=general_font_size)
    return fig, ax


def plot_fragment_length_distribution(array: np.array, start_pos: int, end_pos: int):
    fig = plt.figure(constrained_layout=True, figsize=(15, 5))
    gs = fig.add_gridspec(4, 7)
    ax_heatmap = fig.add_subplot(gs[0:3, 0:])
    axs_distr_plots = [fig.add_subplot(gs[3, i]) for i in range(7)]

    df = pd.DataFrame(array, index=range(start_pos, end_pos))
    seaborn.heatmap(df, ax=ax_heatmap, rasterized=True)

    labels = [int(i.get_text()) - int(df.shape[1] / 2) for i in ax_heatmap.get_xticklabels()]
    ax_heatmap.set_xticklabels(labels)
    min_, max_ = [], []
    for count, i in enumerate(np.array_split(array, 7, axis=1)):
        mean_ax_1 = i.mean(axis=1)
        axs_distr_plots[count].plot(range(start_pos, end_pos), mean_ax_1)
        axs_distr_plots[count].grid()
        max_ += [mean_ax_1.max()]
        min_ += [mean_ax_1.min()]
    for count, i in enumerate(np.array_split(array, 7, axis=1)):
        axs_distr_plots[count].set_ylim((min(min_), max(max_)))
    return fig


def plot_ellipses(ax, weights, means, covars, covariance_type=None):
    """taken from https://scikit-learn.org/stable/auto_examples/mixture/plot_concentration_prior.html#sphx-glr-auto-examples-mixture-plot-concentration-prior-py"""
    for n in range(means.shape[0]):
        if covariance_type == 'full':
            cov_n = covars[n]
        elif covariance_type == 'tied':
            cov_n = covars
        elif covariance_type == 'diag':
            cov_n = np.eye(2) * covars[n]
        elif covariance_type == 'spherical':
            cov_n = np.eye(2) * covars[n]
        else:
            raise ValueError("Invalid covariance_type, must be one of 'spherical', 'tied', 'diag', 'full'")

        eig_vals, eig_vecs = np.linalg.eigh(cov_n)
        unit_eig_vec = eig_vecs[0] / np.linalg.norm(eig_vecs[0])
        angle = np.arctan2(unit_eig_vec[1], unit_eig_vec[0])
        # Ellipse needs degrees
        angle = 180 * angle / np.pi
        # eigenvector normalization
        eig_vals = 2 * np.sqrt(2) * np.sqrt(eig_vals)
        ell = matplotlib.patches.Ellipse(
            means[n], eig_vals[0], eig_vals[1], angle=180 + angle, edgecolor="black"
        )
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(weights[n])
        ell.set_facecolor("#56B4E9")
        ax.add_artist(ell)


