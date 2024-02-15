# Package name
NAME = "atlasreader"

# Version of the package
VERSION = "0.3.1"

# Short description of the package
DESCRIPTION = "A toolbox for generating cluster reports from statistical maps."

# Long description of the package
LONG_DESCRIPTION = """A toolbox for generating cluster reports from statistical maps."""

# Package maintainer
MAINTAINER = "Michael Notter"

# Maintainer email
EMAIL = "michaelnotter@hotmail.com"

# URL to the project homepage
URL = "https://github.com/miykael/atlasreader"

# URL to download the project
DOWNLOAD_URL = "https://github.com/miykael/atlasreader/releases"

# List of dependencies
INSTALL_REQUIRES = [
    "matplotlib>=3.7",
    "nibabel>=5.0",
    "nilearn>=0.10",
    "numpy>=1.22.0",
    "pandas>=2.0",
    "scikit-image>=0.21.0",
    "scikit-learn>=1.0",
    "scipy>=1.10"
]

# Package data to include
PACKAGE_DATA = {
    "atlasreader": ["data/atlas/*"],
}

# Dependencies for running tests
TESTS_REQUIRE = [
    "coverage",
    "pytest>=6.0.0",
    "pytest-cov",
    "nbval"
]

# License of the package
LICENSE = "BSD-3-Clause"
