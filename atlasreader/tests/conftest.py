import pytest
from pathlib import Path

@pytest.fixture
def stat_img():
    return Path(__file__).parent / 'data' / 'collection_658' / "image_10426.nii.gz"