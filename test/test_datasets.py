# Copyright 2026 Mark Diekhans
"""Tests of the dataset catalog loader."""
import copy
import json

import pytest

from flair_validate import FlairValidateDataError
from flair_validate.datasets import (Access, FileRole, Library, Platform, Validation,
                                     load_datasets)

WTC11 = {
    "name": "wtc11-chr22-pb",
    "description": "WTC11 PacBio cDNA from LRGASP, alignments subset to chr22.",
    "access": "public",
    "organism": "human",
    "assembly": "GRCh38",
    "annotation": "gencode-v38",
    "platform": "pacbio",
    "library": "cdna",
    "regions": ["chr22"],
    "validations": ["transcriptome", "splice-junctions", "regression"],
    "source": "ENCFF370NFS",
    "files": {
        "alignments": "data/wtc11-chr22/WTC11.ENCFF370NFS.chr22.genomealigned.bam",
        "genome-fasta": "data/wtc11-chr22/GRCh38.chr22.genome.fa",
        "reads": ["data/wtc11-chr22/rep1.fastq.gz", "data/wtc11-chr22/rep2.fastq.gz"],
    },
}

SMARCA4 = {
    "name": "smarca4-pb",
    "description": "SMARCA4 knockdown PacBio Kinnex, Brooks lab, unpublished.",
    "access": "restricted",
    "organism": "human",
    "assembly": "GRCh38",
    "platform": "pacbio",
    "validations": ["performance"],
    "files": {"alignments": "data/smarca4/B1A_kd_induced_rep1.bam"},
}


def write_catalog(tmp_path, *datasets):
    "write a catalog holding the given dataset objects, return its path"
    metadata_dir = tmp_path / "metadata"
    metadata_dir.mkdir(exist_ok=True)
    datasets_json = metadata_dir / "datasets.json"
    datasets_json.write_text(json.dumps({"datasets": list(datasets)}))
    return datasets_json


def load_bad(tmp_path, *datasets):
    "load a catalog expected to be rejected, return the message"
    with pytest.raises(FlairValidateDataError) as excinfo:
        load_datasets(write_catalog(tmp_path, *datasets))
    return str(excinfo.value)


def edited(dataset, **changes):
    "copy of a dataset with fields changed; a None value drops the field"
    entry = copy.deepcopy(dataset)
    for field, value in changes.items():
        if value is None:
            del entry[field]
        else:
            entry[field] = value
    return entry


def test_full_entry(tmp_path):
    datasets = load_datasets(write_catalog(tmp_path, WTC11, SMARCA4))
    assert len(datasets) == 2
    wtc11 = datasets.by_name("wtc11-chr22-pb")
    assert wtc11.access is Access.public
    assert wtc11.platform is Platform.pacbio
    assert wtc11.library is Library.cdna
    assert wtc11.annotation == "gencode-v38"
    assert wtc11.regions == ("chr22",)
    assert wtc11.validations == frozenset((Validation.transcriptome, Validation.splice_junctions,
                                           Validation.regression))
    assert not wtc11.restricted
    assert not wtc11.whole_genome


def test_optional_fields_omitted(tmp_path):
    datasets = load_datasets(write_catalog(tmp_path, SMARCA4))
    smarca4 = datasets.by_name("smarca4-pb")
    assert smarca4.annotation is None
    assert smarca4.library is None
    assert smarca4.source is None
    assert smarca4.regions == ()
    assert smarca4.whole_genome
    assert smarca4.restricted


def test_paths_resolve_against_root(tmp_path):
    datasets = load_datasets(write_catalog(tmp_path, WTC11))
    wtc11 = datasets.by_name("wtc11-chr22-pb")
    assert wtc11.path(FileRole.genome_fasta) == tmp_path / "data/wtc11-chr22/GRCh38.chr22.genome.fa"
    assert len(wtc11.paths(FileRole.reads)) == 2
    assert wtc11.paths(FileRole.cage_peaks) == ()


def test_path_rejects_multi_valued_role(tmp_path):
    datasets = load_datasets(write_catalog(tmp_path, WTC11))
    with pytest.raises(FlairValidateDataError, match="has 2 files with role 'reads'"):
        datasets.by_name("wtc11-chr22-pb").path(FileRole.reads)


def test_missing_reports_absent_files(tmp_path):
    datasets = load_datasets(write_catalog(tmp_path, SMARCA4))
    smarca4 = datasets.by_name("smarca4-pb")
    assert smarca4.missing() == (tmp_path / "data/smarca4/B1A_kd_induced_rep1.bam",)
    (tmp_path / "data/smarca4").mkdir(parents=True)
    (tmp_path / "data/smarca4/B1A_kd_induced_rep1.bam").touch()
    assert smarca4.missing() == ()


def test_selection(tmp_path):
    datasets = load_datasets(write_catalog(tmp_path, WTC11, SMARCA4))
    assert [d.name for d in datasets.public()] == ["wtc11-chr22-pb"]
    assert [d.name for d in datasets.for_validation(Validation.performance)] == ["smarca4-pb"]
    assert datasets.for_validation(Validation.quantification) == ()


def test_unknown_name(tmp_path):
    datasets = load_datasets(write_catalog(tmp_path, WTC11))
    with pytest.raises(FlairValidateDataError, match="no dataset named 'nope'"):
        datasets.by_name("nope")


def test_catalog_not_found(tmp_path):
    with pytest.raises(FlairValidateDataError, match="can't read dataset catalog"):
        load_datasets(tmp_path / "metadata/datasets.json")


def test_invalid_json(tmp_path):
    datasets_json = write_catalog(tmp_path, WTC11)
    datasets_json.write_text("{not json")
    with pytest.raises(FlairValidateDataError, match="is not valid JSON"):
        load_datasets(datasets_json)


def test_missing_field(tmp_path):
    assert "missing required field(s) 'assembly'" in load_bad(tmp_path, edited(WTC11, assembly=None))


def test_unknown_field(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, genome="GRCh38"))
    assert "unknown field(s) 'genome'" in msg
    assert "permitted fields are" in msg


def test_bad_enum_value(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, platform="nanopore"))
    assert "field 'platform' has 'nanopore'" in msg
    assert "'ont'" in msg


def test_bad_validation_value(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, validations=["transcriptome", "speed"]))
    assert "field 'validations' has 'speed'" in msg


def test_empty_validations(tmp_path):
    assert "must name at least one" in load_bad(tmp_path, edited(WTC11, validations=[]))


def test_bad_name(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, name="WTC11 chr22"))
    assert "name 'WTC11 chr22' must match" in msg


def test_duplicate_name(tmp_path):
    msg = load_bad(tmp_path, WTC11, edited(WTC11, description="a second copy"))
    assert "'wtc11-chr22-pb' used more than once" in msg


def test_unknown_file_role(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, files={"isoforms": "data/x.bed"}))
    assert "unknown role 'isoforms'" in msg


def test_empty_files(tmp_path):
    assert "must name at least one file" in load_bad(tmp_path, edited(WTC11, files={}))


def test_absolute_path_rejected(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, files={"alignments": "/private/groups/brookslab/x.bam"}))
    assert "must be relative and start with 'data/'" in msg


def test_path_outside_data_rejected(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, files={"alignments": "results/x.bam"}))
    assert "must be relative and start with 'data/'" in msg


def test_parent_path_rejected(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, files={"alignments": "data/../../x.bam"}))
    assert "must be relative and start with 'data/'" in msg


def test_source_path_rejected(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, source="/private/groups/brookslab/cafelton"))
    assert "never a file path" in msg


def test_wrong_type(tmp_path):
    assert "must be a string, not int" in load_bad(tmp_path, edited(WTC11, organism=1))


def test_regions_not_a_list(tmp_path):
    msg = load_bad(tmp_path, edited(WTC11, regions="chr22"))
    assert "must be a list of strings, not str" in msg
