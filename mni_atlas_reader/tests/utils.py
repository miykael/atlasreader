"""
Utility functions for testing mni_atlas_reader.
"""

from os.path import join as pjoin
from pkg_resources import resource_filename


def get_test_data_path(fname=None):
    """
    Returns path to test data directory

    If `fname` is supplied, return path to `fname` in test data directory

    Parameters
    ----------
    fname : str, optional
        Filename of test data. Default: None

    Returns
    -------
    path : str
        Path to test data
    """
    path = resource_filename('mni_atlas_reader', 'tests/data')
    return pjoin(path, fname) if fname is not None else path
