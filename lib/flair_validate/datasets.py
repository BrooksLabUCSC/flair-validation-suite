# Copyright 2026 Mark Diekhans
"""
Catalog of the data used to validate FLAIR.

The catalog is ``metadata/datasets.json``, written by hand.  It is a JSON object
with a single field ``datasets``, holding a list of entries::

    {"datasets": [
        {"type": "reference", "name": "grch38-gencode-v48", ...},
        {"type": "sample", "name": "a549-drna", ...},
        {"type": "cage_peaks", "name": "a549-cage", ...}
    ]}

Each entry's ``type`` selects the class that parses it.  The types are:

    ``reference``             Reference
    ``sample``                Sample
    ``splice_junctions``      SpliceJunctions
    ``cage_peaks``            CagePeaks
    ``polya_peaks``           PolyaPeaks
    ``truth_transcriptome``   TruthTranscriptome

Each class documents its own fields; the file kinds they hold are documented in
``flair_validate.datafiles``.  Read both with ``python3 -m pydoc``.

The catalog holds the parts of a run, never a combination of them.  Which
reference, samples and evidence make up a run is a decision, and decisions are
code: an experiment names them by symbol.  So there is no entry that ties a
reference to a sample, and no field saying what a combination can validate.

Names are unique across every type, since an experiment names an entry without
saying what kind it is.
"""
import json
import re
from collections import Counter
from enum import StrEnum
from pathlib import Path

from flair_validate import REPO_ROOT, FlairValidateDataError
from flair_validate.datafiles import (AlignmentsBed, AnnotationGtf, CagePeaksBed, DataFile,
                                      GenomeBam, GenomeFasta, IntronsBed, JunctionsTab,
                                      PolyaPeaksBed, Reads, TruthTranscriptomeGtf)
from flair_validate.jsonparse import (check_keys, get_enum, get_opt_enum, get_opt_str, get_str,
                                      get_str_list, quoted, type_name)

DATASETS_JSON = REPO_ROOT / "metadata/datasets.json"

NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

COMMON_REQUIRED = ("type", "name", "description", "access")
COMMON_OPTIONAL = ("source",)


class Access(StrEnum):
    """Whether the files may be released publicly."""
    public = "public"
    restricted = "restricted"


class Platform(StrEnum):
    """Sequencing platform producing the reads."""
    pacbio = "pacbio"
    ont = "ont"
    illumina = "illumina"
    simulated = "simulated"


class Library(StrEnum):
    """Library preparation: complementary DNA or direct RNA."""
    cdna = "cdna"
    drna = "drna"


class Dataset:
    """Base class of every catalog entry.  Fields every entry has:

    ``type``
        Which kind of entry this is, selecting the class that parses it.
    ``name``
        Symbolic name, matching ``[A-Za-z0-9][A-Za-z0-9._-]*``, unique across the
        whole catalog.  An experiment names an entry by this, and it names output
        directories, so it holds no whitespace and nothing a shell or a path
        reserves.
    ``description``
        English description: what the data is and where it came from.
    ``access``
        ``public`` or ``restricted``.  Restricted data is Brooks lab data that is
        not publicly releasable.
    ``source``
        Optional public provenance: an accession or a URL.  Never a path into
        restricted storage, since this file is committed.
    """

    dataset_type = None
    required = ()
    optional = ()
    regions = ()

    def __init__(self, entry, root_dir, where):
        self.name = _parse_name(entry, where)
        self.description = get_str(entry, "description", where)
        self.access = get_enum(Access, entry, "access", where)
        self.source = _parse_source(entry, where)
        self.root_dir = Path(root_dir)

    def __repr__(self):
        return f"{type(self).__name__}('{self.name}')"

    @classmethod
    def required_fields(cls):
        "names of the fields this type must have"
        return COMMON_REQUIRED + cls.required

    @classmethod
    def optional_fields(cls):
        "names of the fields this type may have"
        return COMMON_OPTIONAL + cls.optional

    @property
    def restricted(self):
        "is access to this entry's data restricted?"
        return self.access is Access.restricted

    @property
    def whole_genome(self):
        "does this entry cover the whole genome?"
        return len(self.regions) == 0

    def data_files(self):
        "every data file this entry names"
        return _collect_files(vars(self).values())

    def missing(self):
        "locations of files that do not exist, which is normal for restricted data"
        return tuple(f.path for f in self.data_files() if not f.exists)


class Reference(Dataset):
    """A genome and the annotation of it that a run is done against.  Written once
    and named by every sample and evidence set aligned to it, which is what makes a
    change of annotation a change of one entry.  Added fields:

    ``assembly``
        Genome assembly, such as ``GRCh38`` or ``sacCer3``.
    ``genome_fasta``
        Genome sequence.
    ``annotation``
        Optional annotation release, such as ``gencode-v38``.
    ``annotation_gtf``
        Optional reference annotation.
    ``regions``
        Optional list of regions the reference is cut down to, such as
        ``["chr22"]``.  Absent means the whole genome.
    """

    dataset_type = "reference"
    required = ("assembly", "genome_fasta")
    optional = ("annotation", "annotation_gtf", "regions")

    def __init__(self, entry, root_dir, where):
        super().__init__(entry, root_dir, where)
        self.assembly = get_str(entry, "assembly", where)
        self.annotation = get_opt_str(entry, "annotation", where)
        self.genome_fasta = GenomeFasta.get(entry, root_dir, where)
        self.annotation_gtf = AnnotationGtf.get(entry, root_dir, where, required=False)
        self.regions = get_str_list(entry, "regions", where, required=False)


class ReferencedDataset(Dataset):
    """Base of the entries whose files carry coordinates on a reference, so they
    name one.  Added field:

    ``reference``
        Name of the ``reference`` entry the files are on.  The loader replaces it
        with the entry itself, available as ``reference``; the name as written stays
        in ``reference_name``.
    """

    def __init__(self, entry, root_dir, where):
        super().__init__(entry, root_dir, where)
        self.reference_name = get_opt_str(entry, "reference", where)
        self.reference = None


class Sample(ReferencedDataset):
    """The sequence data of one biological sample, as reads, as alignments, or both.
    One sample, not a group of them: a run combines samples by naming several.
    Added fields:

    ``platform``
        ``pacbio``, ``ont``, ``illumina`` or ``simulated``.
    ``library``
        Optional ``cdna`` or ``drna``.
    ``reference``
        Required when the sample has alignments, since those carry coordinates.
        Absent is allowed for reads alone, which do not.
    ``reads``
        Unaligned reads, one file or several.
    ``alignments``
        Reads aligned to the genome, one file or several.
    ``alignments_bed``
        The same alignments as BED12.
    ``regions``
        Optional list of regions the data is cut down to.

    At least one of ``reads`` and ``alignments`` is required.
    """

    dataset_type = "sample"
    required = ("platform",)
    optional = ("library", "reference", "reads", "alignments", "alignments_bed", "regions")

    def __init__(self, entry, root_dir, where):
        super().__init__(entry, root_dir, where)
        self.platform = get_enum(Platform, entry, "platform", where)
        self.library = get_opt_enum(Library, entry, "library", where)
        self.reads = Reads.get(entry, root_dir, where, required=False)
        self.alignments = GenomeBam.get(entry, root_dir, where, required=False)
        self.alignments_bed = AlignmentsBed.get(entry, root_dir, where, required=False)
        self.regions = get_str_list(entry, "regions", where, required=False)
        self._check_has_data(where)
        self._check_reference(where)

    @property
    def aligned(self):
        "does this sample come with alignments?"
        return len(self.alignments) > 0

    def _check_has_data(self, where):
        if not (self.reads or self.alignments):
            raise FlairValidateDataError(f"{where}: a sample needs 'reads', 'alignments' or both; "
                                         f"it has neither")

    def _check_reference(self, where):
        if (self.reference_name is None) and (self.alignments or self.alignments_bed):
            raise FlairValidateDataError(f"{where}: field 'reference' is required for a sample "
                                         f"with alignments, since they carry coordinates on one; "
                                         f"name the reference entry they were aligned to")


class SpliceJunctions(ReferencedDataset):
    """Splice junctions, evidence for which junctions are real.  Added fields:

    ``reference``
        Name of the reference the junctions are on.
    ``junctions_tab``
        STAR ``SJ.out.tab``, passed to FLAIR as ``--junction_tab``.
    ``introns_bed``
        intron-prospector introns with read support in the score column, passed as
        ``--junction_bed``.

    The two files are alternative ways to say the same thing, from different tools.
    At least one is required; both are allowed, for the same junctions in the two
    forms.
    """

    dataset_type = "splice_junctions"
    required = ("reference",)
    optional = ("junctions_tab", "introns_bed")

    def __init__(self, entry, root_dir, where):
        super().__init__(entry, root_dir, where)
        self.junctions_tab = JunctionsTab.get(entry, root_dir, where, required=False)
        self.introns_bed = IntronsBed.get(entry, root_dir, where, required=False)
        if (self.junctions_tab is None) and (self.introns_bed is None):
            raise FlairValidateDataError(f"{where}: needs 'junctions_tab', 'introns_bed' or both; "
                                         f"it has neither")


class CagePeaks(ReferencedDataset):
    """CAGE peaks, evidence for transcription start sites.  Added fields:
    ``reference`` and ``cage_peaks``."""

    dataset_type = "cage_peaks"
    required = ("reference", "cage_peaks")

    def __init__(self, entry, root_dir, where):
        super().__init__(entry, root_dir, where)
        self.cage_peaks = CagePeaksBed.get(entry, root_dir, where)


class PolyaPeaks(ReferencedDataset):
    """Poly(A) or QuantSeq peaks, evidence for transcript end sites.  Added fields:
    ``reference`` and ``polya_peaks``."""

    dataset_type = "polya_peaks"
    required = ("reference", "polya_peaks")

    def __init__(self, entry, root_dir, where):
        super().__init__(entry, root_dir, where)
        self.polya_peaks = PolyaPeaksBed.get(entry, root_dir, where)


class TruthTranscriptome(ReferencedDataset):
    """The transcriptome a simulated sample was generated from, the one case where
    accuracy is measured against truth rather than against an annotation.  Added
    fields: ``reference`` and ``truth_transcriptome``."""

    dataset_type = "truth_transcriptome"
    required = ("reference", "truth_transcriptome")

    def __init__(self, entry, root_dir, where):
        super().__init__(entry, root_dir, where)
        self.truth_transcriptome = TruthTranscriptomeGtf.get(entry, root_dir, where)


DATASET_TYPES = (Reference, Sample, SpliceJunctions, CagePeaks, PolyaPeaks, TruthTranscriptome)
_TYPE_CLASSES = {cls.dataset_type: cls for cls in DATASET_TYPES}


class Datasets:
    """The catalog: every entry from one metadata/datasets.json."""

    def __init__(self, datasets, datasets_json):
        self.datasets = tuple(datasets)
        self.datasets_json = datasets_json
        self._by_name = {d.name: d for d in self.datasets}

    def __iter__(self):
        return iter(self.datasets)

    def __len__(self):
        return len(self.datasets)

    def find(self, name):
        "the entry with this name, None when there is none"
        return self._by_name.get(name)

    def by_name(self, name):
        "the entry with this name, whatever its type"
        dataset = self._by_name.get(name)
        if dataset is None:
            raise FlairValidateDataError(f"no dataset named '{name}' in {self.datasets_json}; "
                                         f"known names are {quoted(self._by_name)}")
        return dataset

    def of_type(self, dataset_cls):
        "entries of one type, such as of_type(Sample)"
        return tuple(d for d in self.datasets if isinstance(d, dataset_cls))

    def public(self):
        "entries that are not access restricted"
        return tuple(d for d in self.datasets if not d.restricted)


def load_datasets(datasets_json=DATASETS_JSON, root_dir=None):
    """Load and validate the catalog.  Paths in it are relative to root_dir, which
    defaults to the directory holding the metadata directory."""
    datasets_json = Path(datasets_json)
    if root_dir is None:
        root_dir = datasets_json.parent.parent
    entries = _get_dataset_list(_read_json(datasets_json), datasets_json)
    datasets = [_parse_dataset(e, i, datasets_json, Path(root_dir)) for i, e in enumerate(entries)]
    _check_unique_names(datasets, datasets_json)
    catalog = Datasets(datasets, datasets_json)
    _resolve_references(catalog, datasets_json)
    return catalog


def _read_json(datasets_json):
    try:
        with open(datasets_json) as fh:
            return json.load(fh)
    except OSError as ex:
        raise FlairValidateDataError(f"can't read dataset catalog {datasets_json}") from ex
    except json.JSONDecodeError as ex:
        raise FlairValidateDataError(f"dataset catalog {datasets_json} is not valid JSON") from ex


def _get_dataset_list(doc, datasets_json):
    check_keys(doc, ("datasets",), (), str(datasets_json))
    entries = doc["datasets"]
    if not isinstance(entries, list):
        raise FlairValidateDataError(f"{datasets_json}: field 'datasets' must be a list of "
                                     f"entries, not {type_name(entries)}")
    return entries


def _parse_dataset(entry, index, datasets_json, root_dir):
    where = f"{datasets_json}: dataset {index}"
    if not isinstance(entry, dict):
        raise FlairValidateDataError(f"{where} must be an object, not {type_name(entry)}")
    dataset_cls = _get_dataset_class(entry, where)
    where = f"{datasets_json}: {entry.get('type')} '{entry.get('name')}'"
    check_keys(entry, dataset_cls.required_fields(), dataset_cls.optional_fields(), where)
    return dataset_cls(entry, root_dir, where)


def _get_dataset_class(entry, where):
    if "type" not in entry:
        raise FlairValidateDataError(f"{where}: missing required field 'type'; it must be one of "
                                     f"{quoted(_TYPE_CLASSES)}")
    dataset_cls = _TYPE_CLASSES.get(get_str(entry, "type", where))
    if dataset_cls is None:
        raise FlairValidateDataError(f"{where}: unknown type '{entry['type']}'; types are "
                                     f"{quoted(_TYPE_CLASSES)}")
    return dataset_cls


def _check_unique_names(datasets, datasets_json):
    counts = Counter(d.name for d in datasets)
    dups = sorted(name for name, cnt in counts.items() if cnt > 1)
    if dups:
        raise FlairValidateDataError(f"{datasets_json}: name(s) {quoted(dups)} used more than "
                                     f"once; a name is unique across every type, since an "
                                     f"experiment names an entry without saying its type")


def _resolve_references(catalog, datasets_json):
    for dataset in catalog.of_type(ReferencedDataset):
        if dataset.reference_name is not None:
            dataset.reference = _lookup_reference(dataset, catalog, datasets_json)


def _lookup_reference(dataset, catalog, datasets_json):
    reference = catalog.find(dataset.reference_name)
    if reference is None:
        raise FlairValidateDataError(f"{datasets_json}: {dataset.dataset_type} '{dataset.name}' "
                                     f"names reference '{dataset.reference_name}', which is not "
                                     f"in the catalog")
    if not isinstance(reference, Reference):
        raise FlairValidateDataError(f"{datasets_json}: {dataset.dataset_type} '{dataset.name}' "
                                     f"names reference '{dataset.reference_name}', which is a "
                                     f"{reference.dataset_type} entry, not a reference")
    return reference


def _parse_name(entry, where):
    name = get_str(entry, "name", where)
    if not NAME_RE.match(name):
        raise FlairValidateDataError(f"{where}: name '{name}' must match {NAME_RE.pattern}; it "
                                     f"names directories and is typed on command lines, so use "
                                     f"letters, digits, '.', '_' and '-', starting with a letter "
                                     f"or digit")
    return name


def _parse_source(entry, where):
    source = get_opt_str(entry, "source", where)
    if (source is not None) and source.startswith("/"):
        raise FlairValidateDataError(f"{where}: field 'source' is '{source}'; this file is "
                                     f"committed, so source holds an accession or a URL, never a "
                                     f"file path, which may name restricted storage")
    return source


def _collect_files(values):
    "the DataFile objects among an entry's attribute values"
    files = []
    for value in values:
        if isinstance(value, DataFile):
            files.append(value)
        elif isinstance(value, tuple):
            files.extend(v for v in value if isinstance(v, DataFile))
    return tuple(files)
