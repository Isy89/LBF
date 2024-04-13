import polars as pl
import numpy as np
from lbfextract.transcription_factor_analysis.schemas import AccessibilityConfig


def get_chromatin_accessibility_coverage(df_: pl.DataFrame, metadata: AccessibilityConfig) -> np.ndarray:
    """
    Function to extract the chromatin accessibility at TFBSs for given TFs. Since the accessibility is inversely 
    correlated with coverage. The sign of the amplitude is swapped and measured to be the min if there is a dip or
    the max, in case there is a peak in the central part of the signal. the central part is defined in the metadata
    """
    central_part_mean = df_[:, metadata.start:metadata.end].mean(axis=1)

    values_relative_to_1 = df_[:, metadata.start:metadata.end] - 1

    return np.where(central_part_mean > 0,
                    values_relative_to_1.max(axis=1) * -1,
                    values_relative_to_1.min(axis=1) * -1)


def get_chromatin_accessibility_entropy(df_: pl.DataFrame, metadata: AccessibilityConfig) -> np.ndarray:
    """
    Function to extract the chromatin accessibility at TFBSs for given TFs from entropy derived signals.
    """
    values_relative_to_1 = df_[:, metadata.start:metadata.end] - 1

    return values_relative_to_1.max(axis=1)
