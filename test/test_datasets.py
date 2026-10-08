# Copyright 2026 Mark Diekhans
"""Tests of the dataset catalog loader."""
import copy
import json

import pytest

from flair_validate import FlairValidateDataError
from flair_validate.datasets import (Access, CagePeaks, Library, Platform, PolyaPeaks,
                                     Reference, Sample, SpliceJunctions, load_datasets)

REF = {
    "type": "reference",
    "name": "grch38-chr22-gencode-v38",
    "description": "GRCh38 chr22 with the chr22 subset of GENCODE v38.",
    "access": "public",
    "assembly": "GRCh38",
    "annotation": "gencode-v38",
    "regions": ["chr22"],
    "genome_fasta": "data/wtc11-chr22/GRCh38.chr22.genome.fa",
    "annotation_gtf": "data/wtc11-chr22/gencode.v38.annotation.chr22.gtf",
}

SAMPLE = {
    "type": "sample",
    "name": "wtc11-chr22-pb",
    "description": "WTC11 PacBio cDNA from LRGASP, alignments subset to chr22.",
    "access": "public",
    "platform": "pacbio",
    "library": "cdna",
    "reference": "grch38-chr22-gencode-v38",
    "regions": ["chr22"],
    "source": "ENCFF370NFS",
    "alignments": "data/wtc11-chr22/WTC11.chr22.bam",
    "alignments_bed": "data/wtc11-chr22/WTC11.chr22.bed",
}

READS = {
    "type": "sample",
    "name": "wtc11-lrgasp-ont",
    "description": "Full size WTC11 ONT cDNA, three replicates, unaligned.",
    "access": "public",
    "platform": "ont",
    "reads": ["data/wtc11-lrgasp/rep1.fastq.gz", "data/wtc11-lrgasp/rep2.fastq.gz"],
}

JUNCTIONS = {
    "type": "splice_junctions",
    "name": "wtc11-chr22-junctions",
    "description": "WTC11 short read junctions from STAR.",
    "access": "public",
    "reference": "grch38-chr22-gencode-v38",
    "junctions_tab": "data/wtc11-chr22/WTC11_all.SJ.out.tab",
}

CAGE = {
    "type": "cage_peaks",
    "name": "a549-cage",
    "description": "A549 CAGE peaks, TSS evidence.",
    "access": "restricted",
    "reference": "grch38-chr22-gencode-v38",
    "cage_peaks": "data/a549-drna/A549_CAGE.bed",
}

POLYA = {
    "type": "polya_peaks",
    "name": "wtc11-polya",
    "description": "WTC11 poly(A) peaks, TES evidence.",
    "access": "restricted",
    "reference": "grch38-chr22-gencode-v38",
    "polya_peaks": "data/a549-drna/WTC11_all_polyApeaks_fixed.bed",
}


def write_catalog(tmp_path, *entries):
    "write a catalog holding the given entries, return its path"
    metadata_dir = tmp_path / "metadata"
    metadata_dir.mkdir(exist_ok=True)
    datasets_json = metadata_dir / "datasets.json"
    datasets_json.write_text(json.dumps({"datasets": list(entries)}))
    return datasets_json


def load_good(tmp_path, *entries):
    "load a catalog expected to be accepted"
    return load_datasets(write_catalog(tmp_path, *entries))


def load_bad(tmp_path, *entries):
    "load a catalog expected to be rejected, return the message"
    with pytest.raises(FlairValidateDataError) as excinfo:
        load_datasets(write_catalog(tmp_path, *entries))
    return str(excinfo.value)


def edited(entry, **changes):
    "copy of an entry with fields changed; a None value drops the field"
    copied = copy.deepcopy(entry)
    for field, value in changes.items():
        if value is None:
            del copied[field]
        else:
            copied[field] = value
    return copied


def test_reference(tmp_path):
    ref = load_good(tmp_path, REF).by_name("grch38-chr22-gencode-v38")
    assert isinstance(ref, Reference)
    assert ref.access is Access.public
    assert not ref.restricted
    assert ref.assembly == "GRCh38"
    assert ref.annotation == "gencode-v38"
    assert ref.regions == ("chr22",)
    assert not ref.whole_genome
    assert ref.genome_fasta.path == tmp_path / "data/wtc11-chr22/GRCh38.chr22.genome.fa"
    assert ref.annotation_gtf.path == tmp_path / "data/wtc11-chr22/gencode.v38.annotation.chr22.gtf"


def test_reference_minimal(tmp_path):
    entry = edited(REF, annotation=None, annotation_gtf=None, regions=None)
    ref = load_good(tmp_path, entry).by_name(REF["name"])
    assert ref.annotation is None
    assert ref.annotation_gtf is None
    assert ref.whole_genome


def test_sample(tmp_path):
    sample = load_good(tmp_path, REF, SAMPLE).by_name("wtc11-chr22-pb")
    assert isinstance(sample, Sample)
    assert sample.platform is Platform.pacbio
    assert sample.library is Library.cdna
    assert sample.source == "ENCFF370NFS"
    assert sample.aligned
    assert len(sample.alignments) == 1
    assert sample.alignments[0].index_path == tmp_path / "data/wtc11-chr22/WTC11.chr22.bam.bai"
    assert sample.reference.name == "grch38-chr22-gencode-v38"
    assert sample.reference.assembly == "GRCh38"


def test_sample_reads_only(tmp_path):
    sample = load_good(tmp_path, READS).by_name("wtc11-lrgasp-ont")
    assert len(sample.reads) == 2
    assert not sample.aligned
    assert sample.library is None
    assert sample.reference is None
    assert sample.reads[0].compressed


def test_evidence(tmp_path):
    catalog = load_good(tmp_path, REF, JUNCTIONS, CAGE, POLYA)
    junctions = catalog.by_name("wtc11-chr22-junctions")
    assert isinstance(junctions, SpliceJunctions)
    assert junctions.introns_bed is None
    assert junctions.junctions_tab.path.name == "WTC11_all.SJ.out.tab"
    assert isinstance(catalog.by_name("a549-cage"), CagePeaks)
    assert isinstance(catalog.by_name("wtc11-polya"), PolyaPeaks)
    assert catalog.by_name("wtc11-polya").reference.name == REF["name"]


def test_of_type_and_public(tmp_path):
    catalog = load_good(tmp_path, REF, SAMPLE, JUNCTIONS, CAGE)
    assert [d.name for d in catalog.of_type(Sample)] == ["wtc11-chr22-pb"]
    assert [d.name for d in catalog.of_type(Reference)] == ["grch38-chr22-gencode-v38"]
    assert [d.name for d in catalog.public()] == ["grch38-chr22-gencode-v38", "wtc11-chr22-pb",
                                                  "wtc11-chr22-junctions"]
    assert len(catalog) == 4


def test_missing_files(tmp_path):
    bam = tmp_path / "data/wtc11-chr22/WTC11.chr22.bam"
    bam.parent.mkdir(parents=True)
    bam.write_text("")
    sample = load_good(tmp_path, REF, SAMPLE).by_name("wtc11-chr22-pb")
    assert len(sample.data_files()) == 2
    assert sample.missing() == (tmp_path / "data/wtc11-chr22/WTC11.chr22.bed",)


def test_unknown_name(tmp_path):
    with pytest.raises(FlairValidateDataError, match="no dataset named 'nope'"):
        load_good(tmp_path, REF).by_name("nope")


def test_missing_type(tmp_path):
    assert "missing required field 'type'" in load_bad(tmp_path, edited(REF, type=None))


def test_unknown_type(tmp_path):
    msg = load_bad(tmp_path, edited(REF, type="genome"))
    assert "unknown type 'genome'" in msg
    assert "'reference'" in msg


def test_missing_field(tmp_path):
    assert "missing required field(s) 'assembly'" in load_bad(tmp_path, edited(REF, assembly=None))


def test_unknown_field(tmp_path):
    msg = load_bad(tmp_path, edited(REF, platform="pacbio"))
    assert "unknown field(s) 'platform'" in msg


def test_bad_name(tmp_path):
    assert "name 'WTC11 chr22' must match" in load_bad(tmp_path, edited(REF, name="WTC11 chr22"))


def test_duplicate_name_across_types(tmp_path):
    msg = load_bad(tmp_path, REF, edited(SAMPLE, name=REF["name"]))
    assert "'grch38-chr22-gencode-v38' used more than once" in msg


def test_bad_access(tmp_path):
    assert "must be one of 'public', 'restricted'" in load_bad(tmp_path, edited(REF, access="open"))


def test_source_is_a_path(tmp_path):
    msg = load_bad(tmp_path, edited(REF, source="/private/groups/brookslab/thing.fa"))
    assert "never a file path" in msg


def test_path_outside_data(tmp_path):
    entry = edited(REF, genome_fasta="/private/groups/brookslab/GRCh38.fa")
    assert "must be relative and start with 'data/'" in load_bad(tmp_path, entry)


def test_path_with_dotdot(tmp_path):
    entry = edited(REF, genome_fasta="data/../../GRCh38.fa")
    assert "must be relative and start with 'data/'" in load_bad(tmp_path, entry)


def test_bad_suffix(tmp_path):
    entry = edited(REF, genome_fasta="data/wtc11-chr22/GRCh38.chr22.genome.2bit")
    msg = load_bad(tmp_path, entry)
    assert "file 'GRCh38.chr22.genome.2bit' for 'genome_fasta'" in msg
    assert "'.fa'" in msg


def test_path_not_a_string(tmp_path):
    assert "must be a string, not int" in load_bad(tmp_path, edited(REF, genome_fasta=3))


def test_unresolved_reference(tmp_path):
    msg = load_bad(tmp_path, edited(SAMPLE, reference="nope"))
    assert "names reference 'nope', which is not in the catalog" in msg


def test_reference_is_not_a_reference(tmp_path):
    msg = load_bad(tmp_path, REF, READS, edited(SAMPLE, reference="wtc11-lrgasp-ont"))
    assert "which is a sample entry, not a reference" in msg


def test_sample_without_data(tmp_path):
    entry = edited(SAMPLE, alignments=None, alignments_bed=None)
    assert "needs 'reads', 'alignments' or both" in load_bad(tmp_path, REF, entry)


def test_aligned_sample_without_reference(tmp_path):
    entry = edited(SAMPLE, reference=None)
    assert "'reference' is required for a sample with alignments" in load_bad(tmp_path, entry)


def test_junctions_without_a_form(tmp_path):
    entry = edited(JUNCTIONS, junctions_tab=None)
    assert "needs 'junctions_tab', 'introns_bed' or both" in load_bad(tmp_path, REF, entry)


def test_introns_bed(tmp_path):
    entry = edited(JUNCTIONS, junctions_tab=None,
                   introns_bed="data/wtc11-chr22/wtc11.chr22.introns.bed")
    junctions = load_good(tmp_path, REF, entry).by_name(JUNCTIONS["name"])
    assert junctions.junctions_tab is None
    assert junctions.introns_bed.path.name == "wtc11.chr22.introns.bed"


def test_not_json(tmp_path):
    metadata_dir = tmp_path / "metadata"
    metadata_dir.mkdir()
    datasets_json = metadata_dir / "datasets.json"
    datasets_json.write_text("{nope")
    with pytest.raises(FlairValidateDataError, match="is not valid JSON"):
        load_datasets(datasets_json)


def test_unreadable(tmp_path):
    with pytest.raises(FlairValidateDataError, match="can't read dataset catalog"):
        load_datasets(tmp_path / "metadata/datasets.json")


def test_real_catalog():
    "the committed catalog must load and every reference must resolve"
    catalog = load_datasets()
    assert len(catalog) > 0
    for dataset in catalog:
        if getattr(dataset, "reference_name", None) is not None:
            assert dataset.reference.name == dataset.reference_name
