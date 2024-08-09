import itertools
import logging
import pathlib
from typing import Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn
from PyComplexHeatmap import ClusterMapPlotter, HeatmapAnnotation, anno_simple, anno_scatterplot
from matplotlib.axes import Axes
from matplotlib.figure import Figure
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
                title: str = None,
                title_font_size: int = 20,
                general_font_size: int = 15,
                plot_center_line_label: bool = True,
                line_type: str = "-") -> tuple[Figure, Axes]:
    """
    Plot a signal with optional Savitzky-Golay smoothing.

    This function plots the given signal, optionally applying a Savitzky-Golay
    filter to smooth the data. Various customization options for the plot are
    available, including axis, figure, color, labels, and title settings.

    :param summary: A numpy array containing the signal data to be plotted.
    :param apply_savgol: A boolean flag indicating whether to apply Savitzky-Golay
                         smoothing to the signal. Default is True.
    :param savgol_window_length: The length of the window used for Savitzky-Golay
                                 smoothing. Must be an odd integer. Default is 33.
    :param savgol_polyorder: The order of the polynomial used for Savitzky-Golay
                             smoothing. Default is 3.
    :param ax: A matplotlib Axes object to plot on. If None, a new axis will be created.
    :param fig: A matplotlib Figure object to plot on. If None, a new figure will be created.
    :param color: The color of the plot line. Can be a string (e.g., 'r' or 'blue') or an
                  integer. If None, the default color will be used.
    :param label: A label for the plot line. This will be used in the legend if provided.
    :param label_center_line: A label for the center line, if it is plotted.
    :param title: The title of the plot. If None, no title will be displayed.
    :param title_font_size: The font size of the title. Default is 20.
    :param general_font_size: The general font size for the plot text. Default is 15.
    :param plot_center_line_label: A boolean flag indicating whether to plot the center line label.
                                   Default is True.
    :param line_type: The type of line to plot (e.g., '-', '--', '-.', ':'). Default is '-'.

    :return: A tuple containing the matplotlib Figure and Axes objects used for the plot.
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
    plus_minus = u"\u00B1"
    if std is not None:
        color = "lightgray"
        ax.fill_between(np.arange(summary.shape[0]), summary - std, summary + std, alpha=0.5, color=color,
                        label=plus_minus + "std")

    label_center_line = label_center_line if label_center_line else "TFBS center" if plot_center_line_label else None

    ax.set_xticks([0, summary.shape[0] // 2, summary.shape[0]],
                  labels=["-" + str(summary.shape[0] // 2), 0, "+" + str(summary.shape[0] // 2)])
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
                      signal: str = "coverage",
                      title: str = "coverage at TFBS",
                      figsize: tuple[int, int] = (20, 10),
                      ax: Axes = None,
                      fig: Figure = None,
                      color: str = None,
                      label: str = None,
                      alpha: float = 1,
                      linewidth: int = 1,
                      flanking: int = 1500,
                      xlabel: str = None,
                      top: int = 5,
                      bottom: int = 5,
                      window_center: int = 50,
                      mask_rows: list[bool] = None,
                      max_signals: int = 10,
                      sort_values: bool = True
                      ) -> tuple[Figure, Axes]:
    """
    Plots the signal for the top and bottom BED files based on their peak or deep in the central part of the signal one
    next to the others optionally using Savitzky-Golay filter and various customization options.

    :param df: The DataFrame containing the signal data to plot.
    :param apply_savgol: Whether to apply the Savitzky-Golay filter to smooth the signal. Default is False.
    :param savgol_window_length: The window length parameter for the Savitzky-Golay filter. Default is 11.
    :param savgol_polyorder: The polynomial order parameter for the Savitzky-Golay filter. Default is 3.
    :param signal: The signal label for the y-axis. Default is "coverage".
    :param title: The title for the plot. Default is "coverage at TFBS".
    :param figsize: The size of the figure. Default is (20, 10).
    :param ax: The axes object to use for plotting. If not provided, new axes will be created.
    :param fig: The figure object to use for plotting. If not provided, a new figure will be created.
    :param color: The color of the plot line. Default is None.
    :param label: The label for the plot line. Default is None.
    :param alpha: The transparency level of the plot line. Default is 1.
    :param linewidth: The width of the plot line. Default is 1.
    :param flanking: The number of flanking base pairs to consider in the analysis. Default is 1500.
    :param xlabel: The label for the x-axis. Default is None.
    :param top: The number of top rows to consider for analysis. Default is 5.
    :param bottom: The number of bottom rows to consider for analysis. Default is 5.
    :param window_center: The number of base pairs to consider around the center for amplitude calculation. Default is
                          50.
    :param mask_rows: The mask rows to apply to the DataFrame. Default is None.
    :param max_signals: It describes the maximum number of BED files signals to be plotted.
    :param sort_values: It describes whether the profiles should be sorted based on their amplitude
    :return: A tuple containing the figure and axes objects (fig, ax).
    """

    df = df.copy()
    n_intervals = df.shape[0]
    bps = df.shape[1]
    n = max_signals

    if (top + bottom > n >= n_intervals) or (n >= top + bottom >= n_intervals):
        top, bottom = n_intervals // 2, n_intervals // 2

    if (top + bottom > n_intervals > n) or (n_intervals > top + bottom > n):
        top, bottom = n // 2, n // 2

    two_fifths = (bps // 5) * 2

    if flanking >= two_fifths:
        logger.warning(f"flanking is too large, {two_fifths} will be used instead")
        flanking = two_fifths

    if window_center > two_fifths // 4:
        logger.warning(f"window_center is too large, {two_fifths // 4} will be used instead")
        window_center = two_fifths // 4

    mask_indices = np.arange(bps)
    mask = np.logical_or(
        mask_indices <= flanking,
        mask_indices > bps - flanking
    )
    center = bps // 2
    df["amplitude"] = (
            df.iloc[:, center - window_center: center + window_center]
            .mean(axis=1) / df.iloc[:, mask].mean(axis=1)
    )
    indices_intervals = np.arange(n_intervals)
    if mask_rows is None:
        mask_rows = np.logical_or(
            indices_intervals < bottom,
            indices_intervals >= n_intervals - top
        )
    df["orig_row_number"] = np.arange(df.shape[0])
    df = df.sort_values(by="amplitude", ascending=False)
    df = df.iloc[mask_rows, :]
    df = df.sort_values(by="orig_row_number", ascending=True).iloc[:, :-3] if not sort_values else df.iloc[:, :-3]
    n_intervals = df.shape[0]
    bps = df.shape[1]
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

    ax.plot(
        df.values.flatten(),
        label=label if label else '',
        c=color,
        alpha=alpha,
        linewidth=linewidth
    )

    index_list = df.index.to_list()
    left_xtick = [str(int((-bps - 1) / 2))]
    middle_xtick = [plus_minus + str(int((bps + 1) / 2))]
    right_xtick = list(itertools.chain(
        *[[i[0], i[1]] for i in itertools.zip_longest(index_list, middle_xtick * n_intervals)]
    ))
    right_xtick[-1] = "+" + str(int((+bps + 1) / 2))
    xtick_labels = (left_xtick + right_xtick)
    xticks = [int(i) for i in range(0, n_intervals * bps + 1, int(bps / 2))]

    ax.set_xticks(xticks, labels=xtick_labels)
    colors = seaborn.color_palette(n_colors=n_intervals)

    for count, i in enumerate([i for i in range(bps // 2, n_intervals * bps, bps)]):
        ax.axvline(i, ls="-.", color=colors[count], alpha=0.5, linewidth=linewidth)

    for count, i in enumerate([i for i in range(0, n_intervals * bps + 1, bps)]):
        ax.axvline(i, ls=":", color="gray", alpha=0.5)

    ax.spines[['right', 'top', 'bottom']].set_visible(False)
    ax.annotate('', xy=(0, -0.1), xycoords='axes fraction', xytext=(1, -0.1))

    return fig, ax


def plot_heatmap_signal_batch(fig=None, ax=None,
                              array=None,
                              title=None,
                              title_font_size=20,
                              general_font_size=15) -> tuple[Figure, Axes]:
    """
    Plots a heatmap of the given array using seaborn's heatmap function.

    :param fig: The figure object to use for plotting. If not provided, a new figure will be created.
    :param ax: The axes object to use for plotting. If not provided, new axes will be created.
    :param array: The 2D array to plot as a heatmap.
    :param title: The title for the heatmap. If not provided, no title will be set.
    :param title_font_size: The font size for the title. Default is 20.
    :param general_font_size: The font size for the axis labels. Default is 15.
    :return: A tuple containing the figure and axes objects (fig, ax).
    """

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
                               ) -> tuple[Figure, Axes]:
    """
    Plots various analyses of the given array, including a PCA plot, signal plot, signal per position plot,
    and correlation plot.

    :param array: The 2D array to plot and analyze.
    :param fig: The figure object to use for plotting. If not provided, a new figure will be created.
    :param ax: The axes object to use for plotting. If not provided, new axes will be created.
    :param title: The title for the plots. If not provided, no title will be set.
    :param title_font_size: The font size for the title. Default is 20.
    :param general_font_size: The font size for the axis labels. Default is 15.
    :param tf_to_annotate: Transcription factors to annotate on the PCA plot. Default is ("CTCF",).
    :param ylabel: The label for the y-axis of the correlation plot.
    :param flanking: The number of flanking base pairs to consider in the analysis. Default is 1500.
    :param annotation_center_line: The center line annotation for the heatmap. Default is "center".
    :param window_center: The number of base pairs to consider around the center for amplitude calculation.
                          Default is 50.
    :param top: The number of top rows to consider for analysis. Default is 5.
    :param bottom: The number of bottom rows to consider for analysis. Default is 5.
    :return: A tuple containing the figure and axes objects (fig, ax).
    """
    array = array.copy()
    if not ax:
        fig, ax = plt.subplots(2, 2, figsize=(20, 10))

    if title:
        fig.suptitle(title, fontsize=title_font_size)

    pca = PCA(2)
    trf_array = pd.DataFrame(pca.fit_transform(array))
    trf_array.columns = ["PC1", "PC2"]

    bps = array.shape[1]
    n_intervals = array.shape[0]
    n = 10

    if (top + bottom > n >= n_intervals) or (n >= top + bottom >= n_intervals):
        top, bottom = n_intervals // 2, n_intervals // 2

    if (top + bottom > n_intervals > n) or (n_intervals > top + bottom > n):
        top, bottom = n // 2, n // 2

    two_fifths = (bps // 5) * 2

    if flanking > two_fifths:
        logger.warning(f"flanking is too large, {two_fifths} will be used instead")
        flanking = two_fifths

    if window_center > two_fifths / 4:
        logger.warning(f"window_center is too large, {two_fifths / 4} will be used instead")
        window_center = two_fifths / 4

    mask = np.logical_or(np.arange(bps) <= flanking,
                         np.arange(bps) > bps - flanking)
    center = bps // 2

    trf_array["amplitude"] = (
            array.iloc[:, center - window_center: center + window_center]
            .mean(axis=1) / array.iloc[:, mask].mean(axis=1)
    ).values
    trf_array.index = array.index
    if trf_array.shape[0] > 2:
        # with less than 2 rows, seaborn.kdeplot fails
        seaborn.kdeplot(
            data=trf_array, x="PC1", y="PC2", ax=ax[0, 1], color="lightgray")
    seaborn.scatterplot(data=trf_array, x="PC1", y="PC2", ax=ax[0, 1], sizes="amplitude")
    ax[0, 1].set_title("PCA plot", fontsize=general_font_size)

    array["amplitude"] = trf_array.amplitude.copy()
    indices_intervals = np.arange(n_intervals)

    mask_rows = np.logical_or(
        indices_intervals < bottom,
        indices_intervals >= n_intervals - top
    )
    df = array.sort_values(by="amplitude", ascending=False)[mask_rows].copy()
    tf_to_annotate = tf_to_annotate if tf_to_annotate is not None else trf_array.index.to_list()
    marker_size = 20 if n_intervals < 20 else 10 if n_intervals < 50 else 5
    for tf in tf_to_annotate:
        if tf in trf_array.index:
            x, y, _ = trf_array.loc[tf].to_list()
            ax[0, 1].plot(x, y, marker="o", markersize=marker_size // 2,
                          markeredgecolor="red", markerfacecolor="red")
            ax[0, 1].annotate(tf, xy=(x, y))

    for (point, amplitude) in zip(np.arange(df.shape[0]), df.amplitude.to_list()):
        ax[0, 0].axhline(point, color="lightgray")
        ax[0, 0].plot(amplitude, point, marker="o", markersize=marker_size, markeredgecolor="blue",
                      markerfacecolor="lightblue")
    ax[0, 0].set_yticks(np.arange(df.shape[0]), labels=df.index.to_list(), rotation=0)
    ax[0, 0].set_xlabel("Signal", fontsize=general_font_size)
    ax[0, 0].set_title("Signal plot", fontsize=general_font_size)
    df.drop(columns=["amplitude"], inplace=True)

    fig, _ = create_heatmap_bed_per_position(df=df, fig=fig, ax=ax[1, 0], title="", ylabel=None,
                                             general_font_size=general_font_size,
                                             title_font_size=title_font_size,
                                             annotation_center_line=annotation_center_line)
    ax[1, 0].collections[0].colorbar.set_label('Signal scale', fontsize=general_font_size)
    ax[1, 0].set_title("Signal per position plot", fontsize=general_font_size)

    fig, _ = create_masked_corr_matrix_plot(df, fig=fig, ax=ax[1, 1], title=None, title_font_size=title_font_size,
                                            general_font_size=general_font_size,
                                            axlabel=ylabel,
                                            )
    ax[1, 1].collections[0].colorbar.set_label('Correlation scale', fontsize=general_font_size)
    ax[1, 1].set_title("Correlation plot", fontsize=general_font_size)
    for row in range(2):
        for col in range(2):
            if row == 0 and col == 0:
                ax[row, col].spines[['right', 'top']].set_visible(False)

            if row == 0 and col == 1:
                ax[row, col].spines[['right', 'top']].set_visible(False)
    fig.tight_layout()

    return fig, ax


def correlation_map_plot(df: pd.DataFrame) -> Figure:
    """
    This function creates a plot of the correlation matrix between the provided signals.

    :params df: pandas DataFrame containing the signal per BED file to be plotted
    :return: A matplotlib.Figure object containing the plot
    """

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


def create_masked_corr_matrix_plot(df: pd.DataFrame,
                                   fig: Figure = None,
                                   ax: Axes = None,
                                   title: str = None,
                                   title_font_size: int = 20,
                                   general_font_size: int = 15,
                                   axlabel: str = None,
                                   label_color_bar: str = "") -> tuple[Figure, Axes]:
    """
    Creates a masked correlation matrix plot.

    :param df: The DataFrame to compute and plot the correlation matrix from.
    :param fig: The figure object to use for plotting. If not provided, a new figure will be created.
    :param ax: The axes object to use for plotting. If not provided, new axes will be created.
    :param title: The title for the plot. If not provided, no title will be set.
    :param title_font_size: The font size for the title. Default is 20.
    :param general_font_size: The font size for the axis labels. Default is 15.
    :param axlabel: The label for the x and y axes. Default is None.
    :param label_color_bar: The label for the color bar. Default is an empty string.

    :return: A tuple containing the figure and axes objects (fig, ax).
    """

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


def create_heatmap_bed_per_position(df,
                                    fig: Figure = None,
                                    ax: Axes = None,
                                    title: str = "",
                                    ylabel: str = None,
                                    general_font_size: int = 15,
                                    title_font_size: int = 20,
                                    annotation_center_line: str = None,
                                    label_color_bar: str = "") -> tuple[Figure, Axes]:
    """
    Creates a heatmap of the given DataFrame containing the signals for each BED file with optional annotations and
    customizations.

    :param df: The DataFrame to plot as a heatmap.
    :param fig: The figure object to use for plotting. If not provided, a new figure will be created.
    :param ax: The axes object to use for plotting. If not provided, new axes will be created.
    :param title: The title for the plot. Default is an empty string.
    :param ylabel: The label for the y-axis. Default is None.
    :param general_font_size: The font size for the axis labels. Default is 15.
    :param title_font_size: The font size for the title. Default is 20.
    :param annotation_center_line: The label for the vertical center line annotation. Default is None.
    :param label_color_bar: The label for the color bar. Default is an empty string.

    :return: A tuple containing the figure and axes objects (fig, ax).
    """

    if not ax or not fig:
        fig, ax = plt.subplots(2, 2, figsize=(20, 10))

    seaborn.heatmap(df, ax=ax, cmap="viridis", cbar_kws={'label': label_color_bar,
                                                         'orientation': 'horizontal',
                                                         },
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


def plot_fragment_length_distribution(array: np.array,
                                      start_pos: int,
                                      end_pos: int,
                                      title: str = "FLD per position") -> Figure:
    """
    Plots the fragment length distribution per position.

    :param array: The array containing the fragment length data.
    :param start_pos: The starting position for the distribution plot.
    :param end_pos: The ending position for the distribution plot.
    :param title: The title for the plot. Default is "FLD per position".

    :return: The figure object containing the plots.
    """
    general_font_size = 15
    small_font_size = 12
    fig = plt.figure(constrained_layout=True, figsize=(15, 10))
    fig.suptitle(title, ha="center", fontsize=30)
    gs = fig.add_gridspec(10, 7)
    ax_heatmap = fig.add_subplot(gs[3:10, 0:])
    axs_distr_plots = [fig.add_subplot(gs[1:3, i]) for i in range(7)]

    df = pd.DataFrame(array, index=range(start_pos, end_pos))
    seaborn.heatmap(df, ax=ax_heatmap, rasterized=True, cbar_kws={
        "orientation": "horizontal",
    })
    ax_heatmap.set_yticks([0, array.shape[0] // 2, array.shape[0]], [start_pos, (start_pos + end_pos) // 2, end_pos],
                          fontsize=general_font_size)
    ax_heatmap.set_xticks([0, array.shape[-1] // 2, array.shape[-1]], [-array.shape[-1] // 2, 0, +array.shape[-1] // 2],
                          fontsize=general_font_size)
    ax_heatmap.set_xlabel("Position around genomic interval center", fontsize=20)
    ax_heatmap.set_ylabel("Fragment Length Distribution", fontsize=20)
    ax_heatmap.collections[0].colorbar.set_label('Density Scale', fontsize=small_font_size)
    min_, max_ = [], []
    for count, i in enumerate(np.array_split(array, 7, axis=1)):
        mean_ax_1 = i.mean(axis=1)
        axs_distr_plots[count].plot(range(start_pos, end_pos), mean_ax_1)

        axs_distr_plots[count].grid()
        max_ += [mean_ax_1.max()]
        min_ += [mean_ax_1.min()]
    pos = 0
    for count, i in enumerate(np.array_split(array, 7, axis=1)):
        if count == 3:
            axs_distr_plots[count].set_xlabel("fragment length", fontsize=small_font_size)
        if count == 0:
            axs_distr_plots[count].set_ylabel("density", fontsize=small_font_size)

        next_pos = pos + i.shape[1]
        axs_distr_plots[count].set_ylim((min(min_), max(max_)))
        axs_distr_plots[count].set_xticks([start_pos, end_pos],
                                          [start_pos, end_pos],
                                          fontsize=small_font_size)
        y_ticks = axs_distr_plots[count].get_yticks()
        y_ticks_labels = axs_distr_plots[count].get_yticklabels()
        axs_distr_plots[count].set_yticks([y_ticks[0], y_ticks[-2]],
                                          [y_ticks_labels[0], y_ticks_labels[-2]],
                                          fontsize=small_font_size)
        axs_distr_plots[count].set_title(f"\u03BCFLD: [{pos - 2000}, {next_pos - 2000})", fontsize=small_font_size)
        pos = next_pos
    return fig


def plot_heatmap_diff_active_gi(data,
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
    data_for_plot = dict(
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

    return fig, data_for_plot
