import voluptuous
from voluptuous import Schema
from lbfextract.fextract.schemas import Config


class SingleSignalTransformerConfig(Config):
    min_fragment_length = None
    max_fragment_length = None
    signal_transformer = None
    n = None
    subsample = None
    possible_signal_transformers = {"entropy", "fld", "fld_middle", "fld_middle_n", "fld_dyad", "fld_peter_ulz"}
    n_bins_pos = None
    n_bins_len = None
    flip_based_on_strand = None
    gc_correction = None
    tag = None
    w = None
    read_start = None
    read_end = None
    peaks = None
    schema = Schema(
        {
            "flip_based_on_strand":  voluptuous.Coerce(bool, msg="flip_based_on_strand should be a boolean"),
            "min_fragment_length": voluptuous.Coerce(int, msg="n should be a integer"),
            "max_fragment_length": voluptuous.Coerce(int, msg="n should be a integer"),
            "n": voluptuous.All(voluptuous.Coerce(int, msg="n should be a integer"),
                                voluptuous.Range(min=1, msg="n should be greater than 1")),
            "w": voluptuous.Any(None,
                                voluptuous.All(voluptuous.Coerce(int, msg="n should be a integer"),
                                voluptuous.Range(min=1, msg="n should be greater than 1"))),
            "subsample": bool,
            "signal_transformer": voluptuous.validators.In(
                possible_signal_transformers,
                msg="signal_transformer should be one of the following: " + ", ".join(
                    possible_signal_transformers)),
            "n_bins_len": voluptuous.Any(
                None,
                voluptuous.All(
                    voluptuous.Coerce(int, msg="n should be a integer"),
                    voluptuous.Range(min=1, msg="n should be greater than 1")
                )
            ),
            "n_bins_pos": voluptuous.Any(
                None,
                voluptuous.All(
                    voluptuous.Coerce(int, msg="n should be a integer"),
                    voluptuous.Range(min=1, msg="n should be greater than 1")
                )
            ),
            "gc_correction": voluptuous.Coerce(bool, msg="gc_correction should be a boolean"),
            "tag": voluptuous.Coerce(str,
                                     msg="tag should be a string"),
            "read_start": voluptuous.Coerce(int, msg="the start of the region to used of a read"),
            "read_end": voluptuous.Coerce(int, msg="the end of the region to used of a read"),
            "peaks": voluptuous.Coerce(list, msg="peacks should be a boolean"),
        }
    )


