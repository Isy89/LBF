import datetime
import pandas as pd


def generate_time_stamp() -> str:
    """
    Generate a time stamp.
    """
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def remove_outliers(data, threshold=1.5) -> pd.Series:
    """
    Remove outliers from a Pandas Series using the IQR method.

    :param pandas.Series data: The Pandas Series of data.
    :param float threshold: The threshold value to determine outliers (default: 1.5).
    :return: The data with outliers removed.
    """
    q1 = data.quantile(0.25)
    q3 = data.quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - threshold * iqr
    upper_bound = q3 + threshold * iqr
    filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
    return filtered_data
