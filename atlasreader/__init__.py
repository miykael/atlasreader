__all__ = ['__version__', 'create_output', 'get_statmap_info']

from atlasreader.info import __version__
from atlasreader.atlasreader import create_output, get_statmap_info

_LICENSE_MESSAGE = """\
The Python package you are importing, `atlasreader`, is licensed under the
BSD-3 license; however, the atlases it uses are separately licensed under more
restrictive frameworks.

By using `atlasreader`, you agree to abide by the license terms of the
individual atlases. Information on these terms can be found online at
https://github.com/miykael/atlasreader.
"""


def _first_import():
    import os.path as op
    from pkg_resources import resource_filename
    imported = resource_filename('atlasreader', 'data/.imported')
    if not op.exists(imported):
        print(_LICENSE_MESSAGE)
        with open(imported, 'w'):
            pass


_first_import()
