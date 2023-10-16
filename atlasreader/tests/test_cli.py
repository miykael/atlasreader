import subprocess


def test_cli(tmp_path, stat_img):
    output_dir = tmp_path / "mni_test"
    output_dir.mkdir()

    ret = subprocess.run(
        [
            "atlasreader",
            "--atlas",
            "harvard_oxford",
            "aal",
            "--threshold",
            "7.0",
            "--probability",
            "5",
            "--mindist",
            "20",
            "--outdir",
            output_dir,
            stat_img,
            "20",
        ]
    )
    assert ret.returncode == 0

    # test if output exists and if the key .csv and .png files were created
    assert output_dir.exists()
    assert len([x for x in output_dir.iterdir()]) > 0

    stat_img_name = stat_img.stem[:11]
    for ending in ["_clusters.csv", "_peaks.csv", ".png"]:
        assert (output_dir / f"{stat_img_name}{ending}").exists()
