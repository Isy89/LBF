from __future__ import annotations

import pathlib
from typing import Optional

import voluptuous
from voluptuous import ALLOW_EXTRA, Schema


class LbfextractInvalidConfigError(ValueError):
    def __init__(self, class_name: str, x: dict, schema: voluptuous.Schema):
        self.class_name = class_name
        self.x = x
        self.schema = schema
        self.input_as_string = "{\n  " + "\n  ".join([f"{key} : {val}" for key, val in x.items()]) + "\n}"
        msg = f"\n\nInvalid input was found in {self.class_name}\n" \
              f"The input provided was:\n {self.input_as_string}\n" \
              f"The schema that need to be followed is:\n" \
              f"{self.schema}\n\n" \
              f"If you are using LBFextract from the command line check " \
              f"that the provided arguments match the defined schema.\n" \
              f"If you are using LBFextract from the App object, " \
              f"check that the argument passed to the specific config ({self.class_name}) respect the schema."

        super().__init__(msg)


class Config:
    schema = Schema({}, extra=ALLOW_EXTRA)

    def __init__(self, config_dict: Optional[dict] = None):
        if config_dict is None:
            config_dict = {}
        self.__dict__.update(self._validate_input(config_dict))

    def _validate_input(self, x: dict):
        try:
            return self.schema(x)
        except voluptuous.error.Invalid as e:
            raise LbfextractInvalidConfigError(self.__class__.__name__, x, self.schema) from e

    def to_dict(self):
        return self.__dict__

    def __str__(self):
        return "\n".join([f"object of class: {self.__class__.__name__}"] +
                         ["  Attributes:"] +
                         [f"    {key} = {val}" for key, val in self.__dict__.items()])

    def __repr__(self):
        return self.__str__()


class ReadFetcherConfig(Config):
    window: int = None
    flanking_region_window: int = None
    extra_bases: int = None
    n_binding_sites: int = None
    cores: int = None
    f: int = None
    F: int = None

    schema = Schema({
        "window": voluptuous.Coerce(int, msg="window should be a integer"),
        "flanking_region_window": voluptuous.Coerce(int, msg="flanking_region_window should be a integer"),
        "extra_bases": voluptuous.Coerce(int, msg="extra_bases should be a integer"),
        "n_binding_sites": voluptuous.Coerce(int, msg="n_binding_sites should be a integer"),
        "cores": voluptuous.Coerce(int, msg="cores should be a integer"),
        "f": voluptuous.Coerce(int,
                               msg="f should be an integer representing the samtools flag to be used to include reads."
                               ),
        "F": voluptuous.Coerce(int,
                               msg="F should be an integer representing the samtools flag to be used to exclude reads."
                               ),
    }, )


class SingleSignalTransformerConfig(Config):
    signal_transformer = None
    n = None
    window_size = None
    flip_based_on_strand = None
    gc_correction = None
    tag = None
    read_start = None
    read_end = None
    peaks = None
    max_fragment_length = None
    min_fragment_length = None
    possible_signal_transformers = {"coverage",
                                    "coverage_dyads",
                                    "middle_point_coverage",
                                    "middle_n_points_coverage",
                                    "sliding_window_coverage",
                                    "peter_ulz_coverage",
                                    "wps_coverage"}

    schema = Schema(
        {
            "n": voluptuous.Coerce(int, msg="n should be a integer"),
            "window_size": voluptuous.Coerce(int, msg="flanking_window should be a integer"),
            "signal_transformer": voluptuous.validators.In(
                possible_signal_transformers,
                msg="signal_transformer should be one of the following: " + ", ".join(
                    possible_signal_transformers)),
            "flip_based_on_strand": voluptuous.Coerce(bool, msg="flip_based_on_strand should be a boolean"),
            "gc_correction": voluptuous.Coerce(bool, msg="whether gc correction should be performed or not"),
            "tag": voluptuous.Coerce(str,
                                     msg="the bam file tag to be used to extract the gc coefficient from each read"),
            "read_start": voluptuous.Coerce(int, msg="the start of the region to used of a read"),
            "read_end": voluptuous.Coerce(int, msg="the end of the region to used of a read"),
            "peaks": voluptuous.Coerce(list, msg="peaks should be a boolean"),
            "max_fragment_length": voluptuous.Coerce(int, msg="max_fragment_length should be a integer"),
            "min_fragment_length": voluptuous.Coerce(int, msg="min_fragment_length should be a integer"),

        }
    )


class SignalSummarizer(Config):
    bedfile: pathlib.Path = None
    summarization_method: str = None
    schema = Schema(
        {
            "bed_file": voluptuous.Coerce(pathlib.Path, msg="bedfile should be a pathlib.Path"),
            "summarization_method": voluptuous.validators.In(["mean", "median", "max", "min", "skip"]),
        }
    )


class AppExtraConfig(Config):
    ctx: dict = None
    cores: int = None
