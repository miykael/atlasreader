__all__ = ['__version__', 'create_output', 'get_statmap_info']

from atlasreader.info import __version__
from atlasreader.atlasreader import create_output, get_statmap_info

_LICENSE_MESSAGE = """\
The Python package you are importing, AtlasReader, is licensed under the
BSD-3 license; however, the atlases it uses are separately licensed under more
restrictive frameworks.
By using AtlasReader, you agree to abide by the license terms of the
individual atlases. Information on these terms can be found online at:
https://github.com/miykael/atlasreader/tree/master/atlasreader/data
"""

print(_LICENSE_MESSAGE)
