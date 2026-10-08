# Copyright 2026 Mark Diekhans
"""Tests of the FLAIR module wrappers."""
import pytest

from flair_validate import REPO_ROOT, FlairValidateError, flair_modules
from flair_validate.flair_modules import FlairConfig

FLAIR_DIR = REPO_ROOT.parent / "flair-dev"
needs_flair = pytest.mark.skipif(not (FLAIR_DIR / "bin/flair").exists(),
                                 reason=f"no FLAIR checkout at {FLAIR_DIR}")


@pytest.fixture
def config():
    return FlairConfig(FLAIR_DIR)


def test_scalar_option():
    assert flair_modules.option_args({"genome": "g.fa"}) == ["--genome", "g.fa"]


def test_path_option(tmp_path):
    genome = tmp_path / "g.fa"
    assert flair_modules.option_args({"genome": genome}) == ["--genome", str(genome)]


def test_flag_option():
    assert flair_modules.option_args({"no_stringent": True}) == ["--no_stringent"]


def test_omitted_options():
    assert flair_modules.option_args({"junction_tab": None, "no_stringent": False}) == []


def test_list_option():
    assert flair_modules.option_args({"reads": ["a.fq", "b.fq"]}) == ["--reads", "a.fq", "b.fq"]


def test_options_are_ordered():
    args = flair_modules.option_args({"genome": "g", "annot": "a"})
    assert args == ["--annot", "a", "--genome", "g"]


@needs_flair
def test_config(config):
    assert config.flair_dir == FLAIR_DIR
    assert config.prog == FLAIR_DIR / "bin/flair"
    assert config.src_dir == FLAIR_DIR / "src"
    assert repr(config) == f"FlairConfig('{FLAIR_DIR}')"


def test_config_not_a_checkout(tmp_path):
    with pytest.raises(FlairValidateError, match="no flair program at"):
        FlairConfig(tmp_path)


def test_config_without_src(tmp_path):
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin/flair").write_text("")
    with pytest.raises(FlairValidateError, match="no src directory at"):
        FlairConfig(tmp_path)


@needs_flair
def test_command(config, tmp_path):
    cmd = flair_modules.flair_command(config, tmp_path, "transcriptome", {"genome": "g.fa"})
    assert cmd[0] == "/usr/bin/env"
    assert f"--chdir={tmp_path}" in cmd
    assert f"PYTHONPATH={config.src_dir}" in cmd
    timed = cmd.index("/usr/bin/time")
    assert cmd[timed + 1:timed + 3] == ["--verbose", "--output=time"]
    assert cmd[-4:] == [str(config.prog), "transcriptome", "--genome", "g.fa"]


@needs_flair
def test_run(config, tmp_path):
    "flair runs, and the three record files are written"
    stage_dir = tmp_path / "stage"
    flair_modules.transcriptome(config, stage_dir, genome_aligned_bam="x.bam", genome="g.fa",
                                sample_name="s", help=True)
    assert "flair transcriptome" in (stage_dir / flair_modules.CMD_FILE).read_text()
    assert "usage: flair transcriptome" in (stage_dir / flair_modules.LOG_FILE).read_text()
    time_text = (stage_dir / flair_modules.TIME_FILE).read_text()
    assert "Maximum resident set size" in time_text
    assert "Exit status: 0" in time_text


@needs_flair
def test_run_failure(config, tmp_path):
    "a failed run says which stage failed and where its output is"
    stage_dir = tmp_path / "stage"
    with pytest.raises(FlairValidateError, match="flair transcriptome failed in"):
        flair_modules.transcriptome(config, stage_dir, genome_aligned_bam="x.bam", genome="g.fa",
                                    sample_name="s", no_such_option=True)
    assert "unrecognized arguments" in (stage_dir / flair_modules.LOG_FILE).read_text()
