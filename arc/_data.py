"""Internal configuration for ARC's user data directory."""

from platformdirs import user_data_dir


DPATH = user_data_dir("arc")
