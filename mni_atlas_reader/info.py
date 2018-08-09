__version__ = '0.0.1'

NAME = 'mni_atlas_reader'
MAINTAINER = 'Michael Notter'
EMAIL = 'michaelnotter@hotmail.com'
VERSION = __version__
LICENSE = 'MIT'
DESCRIPTION = ('A toolbox for generating cluster reports from statistical '
               'maps in MNI space')
LONG_DESCRIPTION = ('')
URL = 'http://github.com/miykael/{name}'.format(name=NAME)
DOWNLOAD_URL = ('https://github.com/miykael/{name}/archive/{ver}.tar.gz'
                .format(name=NAME, ver=__version__))

INSTALL_REQUIRES = [
    'matplotlib',
    'nibabel',
    'nilearn',
    'numpy',
    'pandas',
    'scipy',
    'scikit-learn'
]

TESTS_REQUIRE = [
    'pytest',
    'pytest-cov'
]

PACKAGE_DATA = {
    'mni_atlas_reader': ['data/*']
}
