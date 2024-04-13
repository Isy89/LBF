import random

import numpy as np
import pandas as pd
import pysam

from lbfextract.fextract.signal_transformer import adapt_indices


class TfbsFragmentLengthDistribution:
    def __init__(self,
                 min_fragment_length: int = 100,
                 max_fragment_length: int = 400,
                 gc_correction: bool = False,
                 tag: str = None):
        self.min_fragment_length = min_fragment_length
        self.max_fragment_length = max_fragment_length
        self.gc_correction = gc_correction
        self.tag = tag

    def get_relative_start_end(self,
                               read: pysam.AlignedSegment,
                               start: int):
        relative_start = read.reference_start - start
        relative_end = relative_start + read.template_length
        return relative_start, relative_end

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        relative_fragment_length_range = self.max_fragment_length - self.min_fragment_length
        tensor = np.zeros((relative_fragment_length_range, region_length))
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_fragment_length = read.template_length - self.min_fragment_length
            if 0 <= relative_fragment_length < relative_fragment_length_range:
                relative_start, relative_end = self.get_relative_start_end(read, start)
                index_col = adapt_indices(relative_start, relative_end, region_length)
                if index_col is not None:
                    tensor[relative_fragment_length, index_col] += 1 * gc_coef
        return tensor


class TfbsFragmentLengthDistributionMiddlePoint(TfbsFragmentLengthDistribution):

    def get_relative_start_end(self, read: pysam.AlignedSegment, start: int):
        relative_start = read.reference_start - start
        relative_end = relative_start + read.template_length
        middle_point = (relative_start + relative_end) // 2
        return middle_point, middle_point + 1


class TfbsFragmentLengthDistributionMiddleNPoints(TfbsFragmentLengthDistribution):

    def __init__(self, min_fragment_length=100, max_fragment_length=400, gc_correction: bool = False, tag: str = None,
                 n=5):
        super().__init__(min_fragment_length, max_fragment_length, gc_correction, tag)
        self.n = n

    def get_relative_start_end(self, read: pysam.AlignedSegment, start: int):
        relative_start = read.reference_start - start
        relative_end = relative_start + read.template_length
        middle_point = (relative_start + relative_end) // 2
        return middle_point - self.n, middle_point + self.n + 1


class TfbsFragmentLengthDistributionDyad(TfbsFragmentLengthDistributionMiddleNPoints):

    def __init__(self, min_fragment_length=100, max_fragment_length=400, gc_correction: bool = False, tag: str = None,
                 n=5, peaks: list = None):
        super().__init__(min_fragment_length, max_fragment_length, gc_correction, tag, n)
        self.peaks = peaks

    def get_relative_start_end(self, read: pysam.AlignedSegment, start: int) -> list:
        relative_start = read.reference_start - start
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
        for i in range(1, n_of_possible_nucleosomes *2, 2):
            n_s = relative_start + (m * i) - self.n
            n_e = relative_start + (m * i) + self.n
            s.append((n_s, n_e))
        return s

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        relative_fragment_length_range = self.max_fragment_length - self.min_fragment_length
        tensor = np.zeros((relative_fragment_length_range, region_length))
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_fragment_length = read.template_length - self.min_fragment_length
            if 0 <= relative_fragment_length < relative_fragment_length_range:
                relative_starts_ends_list = self.get_relative_start_end(read, start)
                for n_s, n_e in relative_starts_ends_list:
                    index_col = adapt_indices(n_s, n_e, region_length)
                    if index_col is not None:
                        tensor[relative_fragment_length, index_col] += 1 * gc_coef
        return tensor


class PeterUlzFragmentLengthDistribution(TfbsFragmentLengthDistribution):

    def __init__(self,
                 min_fragment_length: int,
                 max_fragment_length: int,
                 gc_correction: bool,
                 tag: str,
                 read_start: int = 53,
                 read_end: int = 113
                 ):
        super().__init__(min_fragment_length, max_fragment_length, gc_correction, tag)
        self.read_start = read_start
        self.read_end = read_end

    def __call__(self, x: pd.Series):
        start = x.Start
        end = x.End
        region_length = end - start
        relative_fragment_length_range = self.max_fragment_length - self.min_fragment_length
        tensor = np.zeros((relative_fragment_length_range, region_length))
        for read in x.reads_per_interval:
            gc_coef = read.get_tag(self.tag) if self.gc_correction and read.has_tag(self.tag) else 1
            relative_start = read.pos - start
            relative_fragment_length = read.template_length - self.min_fragment_length
            if 0 <= relative_fragment_length < relative_fragment_length_range:
                relative_end = relative_start + read.template_length
                first_read_dyad_pos = [relative_start + self.read_start, relative_start + self.read_end]
                second_read_dyad_pos = [relative_end - self.read_end, relative_end - self.read_start]
                for relative_start, relative_end in [first_read_dyad_pos, second_read_dyad_pos]:
                    indices = adapt_indices(relative_start, relative_end, region_length)
                    if indices is not None:
                        tensor[relative_fragment_length, indices] += 1 * gc_coef
        return tensor
