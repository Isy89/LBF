import pathlib
from collections import defaultdict
import pandas as pd
import numpy as np

path_here = pathlib.Path(__file__).parent


class DummyDatasetGenerator:
    def __init__(self, n_tf, groups: int = 2, dimensions: int = 4000, n_samples: int = 100, n_diff_active: float = 0.2,
                 random_seed: int = 42, noise_level=50, sign=1):
        self.random_state = np.random.RandomState(seed=random_seed)
        self.n_tf = n_tf
        self.groups = groups
        self.dimensions = dimensions
        self.n_samples = n_samples
        self.n_diff_active = int(self.n_tf * n_diff_active)
        self.tfs = self.pick_tf_names(self.load_tfs())
        self.noise_level = noise_level
        self.sign = sign

    def check_n_tf(self):
        if self.n_tf > len(self.tfs):
            raise ValueError(f"Number of TFs must be less than or equal to {len(self.tfs)}")

    def pick_tf_names(self, tfs: list[str]) -> np.array:
        return self.random_state.choice(tfs, self.n_tf, replace=False)

    def load_tfs(self) -> np.array:
        list_tfs = pd.read_csv(path_here / "tfs.csv", index_col=0).iloc[:, 0].to_list()
        self.random_state.shuffle(list_tfs)
        return list_tfs

    def nucleosome_signal_sim(self, damping_factor: float = 1, peak: float = 1, noise_level: float = 0.1) -> np.array:
        signal = np.hstack([
            (np.sin(np.linspace(20, 85, self.dimensions)) * np.exp(
                -damping_factor * (np.arange(self.dimensions) / self.dimensions)))[::-1],
            (np.sin(np.linspace(20, 85, self.dimensions)) * np.exp(
                -damping_factor * (np.arange(self.dimensions) / self.dimensions)))
        ])
        noise = self.random_state.normal(0, 0.1, self.dimensions * 2)
        return (signal * peak + (noise * noise_level)) + 1

    def generate_similar_tf_signals(self) -> pd.DataFrame:
        results = defaultdict(dict)
        for tf in range(self.n_tf - self.n_diff_active):
            damping_factor = self.random_state.randint(5, 15)
            peak = self.sign * self.random_state.poisson(lam=3, size=1)[0]
            for group in range(self.groups):
                for sample in range(self.n_samples):
                    results[f"group:{group}_sample:{sample}"][tf] = self.nucleosome_signal_sim(
                        damping_factor=damping_factor,
                        peak=peak,
                        noise_level=self.random_state.rand() * self.noise_level
                    )
        return pd.DataFrame(results)

    def generate_different_tf_signals(self) -> pd.DataFrame:
        signal_values = self.sign * self.random_state.poisson(lam=3, size=self.n_diff_active)
        amplitude_per_group_per_tf = {group: np.roll(signal_values, group) for group in range(self.groups)}
        results = defaultdict(dict)
        for tf in range(self.n_diff_active):
            damping_factor = self.random_state.randint(5, 15)
            for group in range(self.groups):
                for sample in range(self.n_samples):
                    results[f"group:{group}_sample:{sample}"][tf] = self.nucleosome_signal_sim(
                        damping_factor=damping_factor,
                        peak=amplitude_per_group_per_tf[group][tf],
                        noise_level=self.random_state.rand() * self.noise_level
                    )
        return pd.DataFrame(results)

    def create_dataset(self):
        tf_signals_similar = self.generate_similar_tf_signals()
        tf_signals_different = self.generate_different_tf_signals()
        df = pd.concat([tf_signals_similar, tf_signals_different])
        n_rows = df.shape[0]
        index = self.tfs[:n_rows]
        return df.set_index(pd.Index(index))

    @staticmethod
    def split_dataset_by_sample(dataset: pd.DataFrame):
        columns = dataset.columns
        return {column: pd.DataFrame(dataset[column].to_list(), index=dataset.index.copy()) for column in columns}

    def save_dataset_by_sample(self, dataset: pd.DataFrame, output_path: str):
        path = pathlib.Path(output_path)
        path.mkdir(parents=True, exist_ok=True)
        for sample, df in self.split_dataset_by_sample(dataset).items():
            df.to_csv(path / f"{sample}.csv")

    @staticmethod
    def save_dataset(dataset: pd.DataFrame, path: pathlib.Path) -> None:
        dataset.to_csv(path / "dataset.csv")
