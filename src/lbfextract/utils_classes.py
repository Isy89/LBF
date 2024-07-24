from __future__ import annotations

import functools
import logging
from abc import ABC
from collections import defaultdict
from timeit import default_timer
from typing import Optional

import numpy as np
from pympler import asizeof


class Signal(ABC):
    def __init__(self, array: np.array, metadata: Optional[dict], tags: Optional[tuple]):
        self.array = array
        self.metadata = metadata
        self.tags = tags


class TimerAndMemoryProfiler(object):
    def __init__(self, title: str = '', timer_logger: logging.Logger = None, debug=False):
        self.title = title
        self.timer = default_timer
        self.logger = timer_logger
        self.debug = debug
        self.elapsed_secs = None
        self.memory_used_by_object = None
        self.f = None

    def __enter__(self):
        if self.debug:
            self.start = self.timer()
        return self

    def __exit__(self, *args):
        if self.debug:
            end = self.timer()
            self.elapsed_secs = end - self.start
            if self.logger:
                self.logger.info(self.__repr__())
            else:
                print(self.__repr__())

    def __repr__(self):
        return f"{self.__class__.__name__} obj:\n" \
               f"{'Profiling ' + self.f.__name__ if self.f else self.title}\n" \
               f"o result memory used: {self.memory_used_by_object}\n" \
               f"o elapsed time in secs: {self.elapsed_secs}\n"

    def __str__(self):
        return self.__repr__()

    def __call__(self, f):
        self.f = f

        @functools.wraps(f)
        def decorated(*args, **kwds):
            with self:
                result = f(*args, **kwds)
                self.memory_used_by_object = f"{asizeof.asizeof(result) / 1024 / 1024} MB"
                return result

        return decorated


class Tracer(object):
    functions_calls = defaultdict(lambda: 0)

    def __init__(self, logger: logging.Logger = None, debug=False):
        self.logger = logger
        self.debug = debug
        self.f = None

    def __repr__(self):
        return "\n".join([f"{k}: {v}" for k, v in self.functions_calls.items()])

    def __str__(self):
        return self.__repr__()

    def __call__(self, f):
        self.f = f
        if not self.debug:
            @functools.wraps(f)
            def decorated(*args, **kwds):
                return f(*args, **kwds)
        else:
            @functools.wraps(f)
            def decorated(*args, **kwds):
                self.functions_calls[self.f.__module__ + "." + self.f.__name__] += 1
                return f(*args, **kwds)

        return decorated


class ExternalModuleFilter(logging.Filter):
    def filter(self, record):
        return record.name.startswith('lbfextract')
