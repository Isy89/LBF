import voluptuous
from voluptuous import Schema
from lbfextract.fextract.schemas import Config


class SignalSummarizer(Config):
    smooth = None
    sigma = None
    flanking = None

    schema = Schema(
        {
            "smooth": bool,
            "flanking": voluptuous.All(voluptuous.Coerce(int, msg="n should be a integer"),
                                       voluptuous.Range(min=0, msg="n should be positive")),
            "sigma": voluptuous.All(voluptuous.Coerce(float, msg="n should be a float"),
                                    voluptuous.Range(min=0, msg="n should be positive"))
        }
    )
