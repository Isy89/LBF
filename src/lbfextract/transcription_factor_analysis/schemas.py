import voluptuous
from lbfextract.fextract.schemas import Config

class AccessibilityConfig(Config):
    start: int
    end: int
    schema = voluptuous.Schema({
        voluptuous.Required("start"): int,
        voluptuous.Required("end"): int

    })