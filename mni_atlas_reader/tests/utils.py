"""
Utility functions for testing mni_atlas_reader.
"""

from os.path import abspath, dirname, join, sep


def get_test_data_path():
    """
    Returns the path to test datasets, terminated with separator.

    Test-related data should be in kept in mni_atlas_reader/tests/data, and can
    be accessed by joining the path returned from this function with the
    desired data filename.

    Based on a function by Yaroslav Halchenko used in Neurosynth.
    """
    return abspath(join(dirname(__file__), 'data') + sep)
