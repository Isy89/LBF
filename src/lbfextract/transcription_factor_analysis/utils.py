import datetime


def generate_time_stamp():
    """
    Generate a time stamp.

    Returns:
        str: The time stamp.
    """
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def remove_outliers(data, threshold=1.5):
    """
    Remove outliers from a Pandas Series using the IQR method.

    Args:
        data (pandas.Series): The Pandas Series of data.
        threshold (float): The threshold value to determine outliers (default: 1.5).

    Returns:
        pandas.Series: The data with outliers removed.
    """
    q1 = data.quantile(0.25)
    q3 = data.quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - threshold * iqr
    upper_bound = q3 + threshold * iqr
    filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
    return filtered_data
