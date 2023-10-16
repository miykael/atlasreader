import pytest
from pathlib import Path


@pytest.fixture
def data_dir():
    return Path(__file__).parent / 'data'


@pytest.fixture
def stat_img(data_dir):
    return data_dir / 'collection_658' / "image_10426.nii.gz"
