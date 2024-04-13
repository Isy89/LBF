import numpy as np
import pandas as pd


def adapt_indices(start, end, length_region):
    start = start if 0 <= start <= length_region else 0
    end = end if 0 <= end <= length_region else length_region if end > length_region else 0
    return np.arange(start, end, dtype=np.int64)


class TfbsEntropy:
    def __init__(self, min_fragment_length=100, max_fragment_length=400):
        self.min_fragment_length = min_fragment_length
        self.max_fragment_length = max_fragment_length

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        relative_fragment_length_range = self.max_fragment_length - self.min_fragment_length
        tensor = np.zeros((relative_fragment_length_range, region_length))
        for read in x.reads_per_interval:
            relative_fragment_length = read.template_length - self.min_fragment_length

            if 0 <= relative_fragment_length < relative_fragment_length_range:
                relative_start = read.pos - start
                indices = adapt_indices(relative_start, relative_start + read.template_length, region_length)
                if indices is not None:
                    tensor[relative_fragment_length, indices] += 1
        return tensor
