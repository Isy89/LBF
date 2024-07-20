import bottleneck
import numpy as np
import pandas as pd

from lbfextract.utils import adapt_indices


class TFBSCoverage:
    def __init__(self, gc_correction: bool = False, tag: str = None):
        self.gc_correction = gc_correction
        self.tag = tag

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        tensor = np.zeros(region_length)
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_start = read.pos - start

            index = adapt_indices(relative_start, relative_start + read.template_length, region_length)
            if index is not None:
                tensor[index] += 1 * gc_coef
        return tensor


class TFBSCoverageAroundDyads:
    def __init__(self, n=1, gc_correction: bool = False, tag: str = None, peaks: list = None):
        self.n = n
        self.gc_correction = gc_correction
        self.tag = tag
        self.peaks = peaks

    def get_relative_start_end(self, read, start) -> list:
        relative_start = read.pos - start
        f = np.abs(read.template_length)
        s = []
        nucleosome_length = self.peaks[0]
        n_of_possible_nucleosomes = f // nucleosome_length
        remainder = f % nucleosome_length
        p_same = ((nucleosome_length - remainder) / nucleosome_length)
        p_next = 1 - p_same
        n_of_possible_nucleosomes = np.random.choice(
            [n_of_possible_nucleosomes,
             n_of_possible_nucleosomes + 1],
            p=[p_same, p_next]
        )
        expanded_fragment_length = n_of_possible_nucleosomes * nucleosome_length if n_of_possible_nucleosomes > 0 else nucleosome_length
        middle_point = f // 2
        relative_middle_point = relative_start + middle_point
        relative_start = relative_middle_point - (expanded_fragment_length // 2)
        if n_of_possible_nucleosomes > 0:
            m = expanded_fragment_length // (n_of_possible_nucleosomes * 2)
        else:
            m = expanded_fragment_length // 2
        for i in range(1, n_of_possible_nucleosomes * 2, 2):
            n_s = relative_start + (m * i) - self.n
            n_e = relative_start + (m * i) + self.n
            s.append((n_s, n_e))
        return s

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        tensor = np.zeros(region_length)
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_starts_ends_list = self.get_relative_start_end(read, start)
            for n_s, n_e in relative_starts_ends_list:
                index_col = adapt_indices(n_s, n_e, region_length)
                if index_col is not None:
                    tensor[index_col] += 1 * gc_coef
        return tensor


class TFBSMiddlePointCoverage:
    def __init__(self, gc_correction: bool = False, tag: str = None):
        self.gc_correction = gc_correction
        self.tag = tag

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        tensor = np.zeros(region_length)
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_start = read.pos - start
            middle_point = relative_start + (read.template_length // 2)
            indices = adapt_indices(middle_point, middle_point + 1, region_length)
            if indices is not None:
                tensor[indices] += 1 * gc_coef
        return tensor


class TFBSNmiddlePointCoverage:

    def __init__(self, n=1, gc_correction: bool = False, tag: str = None):
        self.n = n
        self.gc_correction = gc_correction
        self.tag = tag

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        tensor = np.zeros(region_length)
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_start = read.pos - start
            middle_point = relative_start + (read.template_length // 2)
            indices = adapt_indices(middle_point - self.n, middle_point + self.n, region_length)
            if indices is not None:
                tensor[indices] += 1 * gc_coef
        return tensor


class TFBSSlidingWindowCoverage:

    def __init__(self, window_size: int, gc_correction: bool = False, tag: str = None):
        self.window_size = window_size
        self.gc_correction = gc_correction
        self.tag = tag

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        tensor = np.zeros(region_length)
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_start = read.pos - start
            indices = adapt_indices(relative_start, relative_start + read.template_length, region_length)
            if indices is not None:
                tensor[indices] += 1 * gc_coef
        tensor = bottleneck.move_mean(tensor, window=self.window_size, min_count=1)
        return tensor


class FragmentLengthDistribution:

    def __init__(self, min_fragment_length=100, max_fragment_length=400, gc_correction: bool = False, tag: str = None):
        self.min_fragment_length = min_fragment_length
        self.max_fragment_length = max_fragment_length
        self.gc_correction = gc_correction
        self.tag = tag

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        relative_fragment_length = self.max_fragment_length - self.min_fragment_length
        tensor = np.zeros((relative_fragment_length + 1, region_length))
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            fragment_length = read.template_length

            if self.min_fragment_length <= fragment_length <= self.max_fragment_length:
                relative_start = read.pos - start
                indices = adapt_indices(relative_start, relative_start + read.template_length,
                                        region_length)
                if indices is not None:
                    tensor[relative_fragment_length, indices] += 1 * gc_coef
        return tensor


class PeterUlzCoverage:

    def __init__(self, gc_correction: bool,
                 tag: str,
                 read_start: int = 53,
                 read_end: int = 113):
        self.gc_correction = gc_correction
        self.tag = tag
        self.read_start = read_start
        self.read_end = read_end

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        tensor = np.zeros(region_length)
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_start = read.pos - start

            relative_end = relative_start + read.template_length
            first_read_dyad_pos = [relative_start + self.read_start, relative_start + self.read_end]
            second_read_dyad_pos = [relative_end - self.read_end, relative_end - self.read_start]
            for relative_start, relative_end in [first_read_dyad_pos, second_read_dyad_pos]:
                indices = adapt_indices(relative_start, relative_end, region_length)
                if indices is not None:
                    tensor[indices] += 1 * gc_coef
        return tensor


class WPSCoverage(TFBSCoverage):
    def __init__(self, gc_correction: bool = False, tag: str = None, window_size: int = None,
                 min_fragment_length: int = None, max_fragment_length: int = None):
        super().__init__(gc_correction, tag)
        self.window_size = window_size
        self.min_fragment_length = min_fragment_length
        self.max_fragment_length = max_fragment_length

    def get_minus_one_indices(self, relative_start, relative_end, region_length):
        starting_start = relative_start - self.window_size // 2
        starting_end = relative_start + self.window_size // 2
        ending_start = relative_end - self.window_size // 2
        ending_end = relative_end + self.window_size // 2
        indices_start = adapt_indices(starting_start, starting_end, region_length)
        indices_end = adapt_indices(ending_start, ending_end, region_length)
        return indices_start, indices_end

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        regions_minus_one = np.zeros(region_length)
        regions_plus_one = np.zeros(region_length)
        for read in x.reads_per_interval:
            if not self.min_fragment_length < read.template_length <= self.max_fragment_length:
                continue
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_start = read.pos - start
            relative_end = relative_start + read.template_length
            indices = adapt_indices(relative_start + (self.window_size // 2), relative_end - (self.window_size//2), region_length)
            indices_start, indices_end = self.get_minus_one_indices(relative_start, relative_end, region_length)
            if indices_start is not None:
                regions_minus_one[indices_start] += 1 * gc_coef
            if indices_end is not None:
                regions_minus_one[indices_end] += 1 * gc_coef
            if indices is not None:
                regions_plus_one[indices] += 1 * gc_coef
        wps_g = (regions_plus_one - regions_minus_one) 
        return wps_g - bottleneck.move_median(wps_g, window=100, min_count=1)

