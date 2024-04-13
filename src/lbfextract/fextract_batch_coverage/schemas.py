from typing import List

import voluptuous
from voluptuous import Schema

from lbfextract.fextract.schemas import Config


class PlotConfig(Config):
    title_font_size: int = None
    general_font_size: int = None
    tf_to_annotate: List[str] = None
    ylabel: str = None
    flanking: int = None
    annotation_center_line: str = None
    apply_savgol: bool = None
    savgol_window_length: int = None
    savgol_polyorder: int = None
    color: str = None
    xlabel: str = None
    top: int = None
    bottom: int = None
    window_center: int = None
    schema = Schema(
        {
            "title_font_size": voluptuous.Coerce(int, msg="title_font_size should be a integer"),
            "general_font_size": voluptuous.Coerce(int, msg="general_font_size should be a integer"),
            "tf_to_annotate": List[str],
            "ylabel": str,
            "flanking": voluptuous.Coerce(int, msg="flanking should be a integer"),
            "annotation_center_line": str,
            "apply_savgol": bool,
            "savgol_window_length": voluptuous.Coerce(int, msg="savgol_window_length should be a integer"),
            "savgol_polyorder": voluptuous.Coerce(int, msg="savgol_polyorder should be a integer"),
            "color": str,
            "xlabel": str,
            "window_center": voluptuous.Coerce(int, msg="window_center should be a integer"),
            "top": voluptuous.Coerce(int, msg="top should be a integer"),
            "bottom": voluptuous.Coerce(int, msg="bottom should be a integer"),

        }
    )
