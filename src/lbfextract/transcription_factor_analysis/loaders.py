import hashlib
import logging
import pathlib
from collections import defaultdict
from copy import copy
from typing import Callable

import numpy as np
import pandas as pd
import polars as pl

from lbfextract.fextract.schemas import Config
from lbfextract.transcription_factor_analysis.accessibility_extraction import get_chromatin_accessibility_coverage, \
    get_chromatin_accessibility_entropy
from lbfextract.transcription_factor_analysis.schemas import AccessibilityConfig

logger = logging.getLogger(__name__)


class ResultsLoader:
    signal_types = defaultdict(dict, {
        "coverage": {"fun": get_chromatin_accessibility_coverage,
                     "validator": AccessibilityConfig},
        "entropy": {"fun": get_chromatin_accessibility_entropy,
                    "validator": AccessibilityConfig}
    })

    @classmethod
    def register_signal_type(cls, signal_name, fun: Callable, validator: Config):
        cls.signal_types[signal_name]["fun"] = fun
        cls.signal_types[signal_name]["validator"] = validator

    def __init__(self, path_to_res_summary: pathlib.Path,
                 accessibility_extraction_config: dict,
                 signal_length: int = 4000,
                 flanking_signal_indices: tuple = (1000, 3000),
                 normalize: bool = False,
                 path_to_sample_sheet: pathlib.Path = None,
                 grouping_column: str = None,
                 signal_type: str = "coverage"):

        self.signal_type = self.check_signal_type_compatibility(signal_type)
        self.path_to_res_summary = path_to_res_summary
        self.signal_length = signal_length
        self.accessibility_extraction_config = accessibility_extraction_config or dict(start=1800, end=2200)
        self.flanking_signal_indices = flanking_signal_indices
        self.normalize = normalize
        self.grouping_column = grouping_column
        self.sample_sheet = self.load_sample_sheet(
            path_to_sample_sheet) if path_to_sample_sheet else self.generate_sample_sheet()

    def check_signal_type_compatibility(self, signal_type: str) -> str:
        signal_types_keys = list(self.signal_types.keys())
        check = True if signal_type in signal_types_keys else False
        if not check:
            raise ValueError(
                f"signal: {self.signal_type} is not compatible, possible signals: {' '.join(list(self.signal_types.keys()))}")
        return signal_type

    def load_sample_sheet(self, path_to_sample_sheet: pathlib.Path) -> pd.DataFrame:
        sample_sheet = pd.read_csv(path_to_sample_sheet, sep=",", index_col=0)
        not_float_cols = (sample_sheet.dtypes.apply(lambda x: x.name) == "object").to_list()
        sample_sheet.loc[:, not_float_cols] = sample_sheet.loc[:, not_float_cols].copy().applymap(
            lambda x: "NA" if pd.isna(x) else x)
        return sample_sheet

    @staticmethod
    def _hash_path_sample_results(sample_path: str | pathlib.Path):
        if isinstance(sample_path, pathlib.Path):
            sample_path = str(sample_path)
        input_bytes = sample_path.encode('utf-8')

        md5_hash = hashlib.md5()
        md5_hash.update(input_bytes)
        hash_result = md5_hash.hexdigest()

        return hash_result

    def generate_sample_sheet(self):
        paths_to_sample_result = list(self.path_to_res_summary.glob("**/*csv"))
        sample_names = [i.parent.parent.stem for i in paths_to_sample_result]
        bed_file_metadata = [i.parent.stem for i in paths_to_sample_result]
        index_df = range(len(paths_to_sample_result))
        sample_sheet_df = pd.DataFrame(
            columns=["group", "tumor_fraction", "cov", "signal_type"],
            index=index_df)
        sample_sheet_df["sample_name"] = pd.Series(sample_names, index=index_df)
        sample_sheet_df["bed_file_metadata"] = pd.Series(bed_file_metadata, index=index_df)
        sample_sheet_df["path_to_res_summary"] = pd.Series(self.path_to_res_summary, index=index_df)
        sample_sheet_df["paths_to_sample_result"] = paths_to_sample_result
        sample_sheet_df.index = sample_sheet_df.paths_to_sample_result.apply(
            lambda x: self._hash_path_sample_results(x)
        ).to_list()
        return sample_sheet_df

    def get_res_df_polars(self) -> pl.DataFrame:
        dfs = []
        for count, i in enumerate(self.path_to_res_summary.glob("**/*csv")):
            path_exists = i.exists()
            path_hash = self._hash_path_sample_results(str(i))
            path_hash_exists = path_hash in self.sample_sheet.index
            if not path_exists or not path_hash_exists:
                msg = f"path {i} {'exists' if path_exists else 'does not exist'}"
                msg1 = f"path hash {path_hash} {'is in index' if path_hash_exists else 'is not in index'}"
                logger.warning(f"{msg + ' but ' + msg1 if path_exists else msg1 + ' but ' + msg}")
                continue
            parsed_name = self.sample_sheet.loc[path_hash]
            index_transcription_factor_df = ["genomic_interval"] + [str(i) for i in list(range(self.signal_length))]
            gi_df = pl.scan_csv(i, new_columns=index_transcription_factor_df).collect()
            parsed_sample_name_df = pl.DataFrame([
                copy(parsed_name.to_list()) for _ in range(gi_df.shape[0])],
                schema=parsed_name.index.to_list()
            )
            dfs.append(pl.concat([parsed_sample_name_df, gi_df], how="horizontal"))
        return pl.concat(dfs, how="vertical")

    def get_accessibility(self, df_: pl.DataFrame) -> np.ndarray:
        metadata = self.signal_types[self.signal_type]["validator"](self.accessibility_extraction_config)
        return self.signal_types[self.signal_type]["fun"](df_, metadata)

    def normalize_df(self, df: pl.DataFrame) -> pl.DataFrame:
        signal_col_index = [str(i) for i in list(range(self.signal_length))]
        signal_col_index_l_flank = [str(i) for i in list(range(self.flanking_signal_indices[0]))]
        signal_col_index_r_flank = [
            str(i) for i
            in list(range(self.flanking_signal_indices[1], self.signal_length))
        ]
        df[signal_col_index] = (
                df[signal_col_index] /
                (
                        0.5 * (df[signal_col_index_l_flank].mean(axis=1) + df[signal_col_index_r_flank].mean(axis=1))
                )
        )
        return df

    def load(self) -> pl.DataFrame:
        res_polar_df = self.normalize_df(self.get_res_df_polars()) if self.normalize else self.get_res_df_polars()
        signal_col_index = [str(i) for i in list(range(self.signal_length))]
        res_polar_df = res_polar_df.with_columns(
            pl.Series(self.get_accessibility(res_polar_df[signal_col_index])).alias("amplitude"))
        return res_polar_df
