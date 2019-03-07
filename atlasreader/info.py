__version__ = '0.1.1'

NAME = 'atlasreader'
MAINTAINER = 'Michael Notter'
EMAIL = 'michaelnotter@hotmail.com'
VERSION = __version__
LICENSE = 'MIT'
DESCRIPTION = ('A toolbox for generating cluster reports from statistical '
               'maps')
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
    'scikit-image',
    'scikit-learn'
]

TESTS_REQUIRE = [
    'pytest',
    'pytest-cov'
]

PACKAGE_DATA = {
    'atlasreader': [
        'data/*', 'data/atlases/*', 'data/templates/*'
    ],
    'atlasreader.tests': [
        'data/*'
    ]
}
