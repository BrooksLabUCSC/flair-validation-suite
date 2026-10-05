# Copyright 2026 Mark Diekhans
"""
Catalog of the data sets used to validate FLAIR.

The catalog is ``metadata/datasets.json``, written by hand.  This module is the
only documentation of its format.

The file is a JSON object with a single field ``datasets``, holding a list of
data set objects.  A data set object has these fields:

``name``
    Symbolic name, matching ``[a-z0-9][a-z0-9-]*``, unique in the file.  Used to
    select a data set on a command line and to name output directories, so it
    carries no whitespace and nothing that FASTA or Newick reserves.  Required.
``description``
    English description: what the data is and where it came from.  Required.
``access``
    ``public`` or ``restricted``.  Restricted data is Brooks lab data that is
    not publicly releasable.  Required.
``organism``
    Common name, such as ``human`` or ``yeast``.  Required.
``assembly``
    Genome assembly, such as ``GRCh38`` or ``sacCer3``.  Required.
``annotation``
    Reference annotation, such as ``gencode-v38``.  Omitted when the data set
    has no reference annotation.
``platform``
    ``pacbio``, ``ont``, ``illumina`` or ``simulated``.  Required.
``library``
    ``cdna`` or ``drna``.  Omitted when it does not apply.
``regions``
    List of region strings limiting the data set, such as ``["chr22"]`` or
    ``["chr19:10900001-11100000"]``.  Omitted when the data set covers the whole
    genome.
``validations``
    Non-empty list naming the general functionality this data set can validate,
    from this closed list:

    ``transcriptome``
        Isoform and transcriptome assembly accuracy.
    ``splice-junctions``
        Splice junction accuracy and support.
    ``transcript-ends``
        TSS and TES accuracy; needs ``cage-peaks`` or ``polya-peaks``.
    ``quantification``
        Isoform and gene quantification accuracy.
    ``performance``
        Run time and peak memory.
    ``regression``
        Comparison of results between FLAIR versions.

    Required.
``source``
    Public provenance: an accession or a URL.  Never a path into restricted
    storage, since this file is committed.  Omitted when there is none.
``files``
    Non-empty object mapping a file role to one path or a list of paths, from
    this closed list of roles:

    ``reads``
        Long reads, FASTQ or FASTA.
    ``alignments``
        Reads aligned to the genome, BAM.
    ``alignments-bed``
        The same alignments as BED12, used as read evidence.
    ``genome-fasta``
        Genome sequence.
    ``annotation-gtf``
        Reference annotation.
    ``short-read-junctions``
        Short read splice junctions, such as a STAR ``SJ.out.tab``.
    ``cage-peaks``
        CAGE peaks, TSS evidence.
    ``polya-peaks``
        Poly(A) or QuantSeq peaks, TES evidence.
    ``truth-transcriptome``
        Known transcriptome, for simulated data.
    ``regions-tsv``
        Region list, columns ``chr``, ``start``, ``end``.

    Required.

Every path is relative and starts with ``data/``.  The ``data/`` directory is
not committed; it holds symlinks to the real files, which for restricted data
live outside the repository.  An absolute path, or one outside ``data/``, is
rejected, which is what keeps restricted locations out of git.

Data files are not required to exist when the catalog is loaded, as most users
have no access to the restricted ones.  Use ``Dataset.missing`` to find out.

Example of one entry::

    {"datasets": [
        {
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
                "alignments-bed": "data/wtc11-chr22/WTC11.ENCFF370NFS.chr22.genomealigned.bed",
                "genome-fasta": "data/wtc11-chr22/GRCh38.chr22.genome.fa",
                "annotation-gtf": "data/wtc11-chr22/gencode.v38.annotation.chr22.gtf"
            }
        }
    ]}
"""
import json
import re
from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from types import MappingProxyType

from flair_validate import FlairValidateDataError

REPO_ROOT = Path(__file__).parents[2]
DATASETS_JSON = REPO_ROOT / "metadata/datasets.json"

DATA_SUBDIR = "data"
NAME_RE = re.compile(r"^[a-z0-9][a-z0-9-]*$")


class Access(StrEnum):
    """Whether a data set may be released publicly."""
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


class Validation(StrEnum):
    """General functionality a data set can validate.  See the module
    documentation for what each one covers."""
    transcriptome = "transcriptome"
    splice_junctions = "splice-junctions"
    transcript_ends = "transcript-ends"
    quantification = "quantification"
    performance = "performance"
    regression = "regression"


class FileRole(StrEnum):
    """Role a file plays in a data set.  See the module documentation."""
    reads = "reads"
    alignments = "alignments"
    alignments_bed = "alignments-bed"
    genome_fasta = "genome-fasta"
    annotation_gtf = "annotation-gtf"
    short_read_junctions = "short-read-junctions"
    cage_peaks = "cage-peaks"
    polya_peaks = "polya-peaks"
    truth_transcriptome = "truth-transcriptome"
    regions_tsv = "regions-tsv"


@dataclass(frozen=True)
class Dataset:
    """One data set from the catalog, with its fields validated.  Paths are as
    written in the catalog, relative to root_dir; use path or paths to get
    locations to open."""
    name: str
    description: str
    access: Access
    organism: str
    assembly: str
    annotation: str | None
    platform: Platform
    library: Library | None
    regions: tuple[str, ...]
    validations: frozenset[Validation]
    source: str | None
    files: MappingProxyType
    root_dir: Path

    @property
    def restricted(self):
        "is access to this data set restricted?"
        return self.access is Access.restricted

    @property
    def whole_genome(self):
        "does this data set cover the whole genome?"
        return len(self.regions) == 0

    def suits(self, validation):
        "can this data set be used to validate the given functionality?"
        return validation in self.validations

    def paths(self, role):
        "locations of the files with the given role, empty if it has none"
        return tuple(self.root_dir / p for p in self.files.get(role, ()))

    def path(self, role):
        "location of the one file with the given role"
        paths = self.paths(role)
        if len(paths) != 1:
            raise FlairValidateDataError(f"dataset '{self.name}' has {len(paths)} files with role "
                                         f"'{role}', expected exactly one; use paths() for roles "
                                         f"that may have more than one file")
        return paths[0]

    def missing(self):
        "locations of files that do not exist, which is normal for restricted data"
        return tuple(p for role in self.files for p in self.paths(role) if not p.exists())


class Datasets:
    """The catalog: all data sets from one metadata/datasets.json."""

    def __init__(self, datasets, datasets_json):
        self.datasets = tuple(datasets)
        self.datasets_json = datasets_json
        self._by_name = {d.name: d for d in self.datasets}

    def __iter__(self):
        return iter(self.datasets)

    def __len__(self):
        return len(self.datasets)

    def by_name(self, name):
        "the data set with this name"
        dataset = self._by_name.get(name)
        if dataset is None:
            raise FlairValidateDataError(f"no dataset named '{name}' in {self.datasets_json}; "
                                         f"known datasets are {_quoted(self._by_name)}")
        return dataset

    def for_validation(self, validation):
        "data sets that can validate the given functionality"
        return tuple(d for d in self.datasets if d.suits(validation))

    def public(self):
        "data sets that are not access restricted"
        return tuple(d for d in self.datasets if not d.restricted)


def load_datasets(datasets_json=DATASETS_JSON, root_dir=None):
    """Load and validate the data set catalog.  Paths in it are relative to
    root_dir, which defaults to the directory holding the metadata directory."""
    datasets_json = Path(datasets_json)
    if root_dir is None:
        root_dir = datasets_json.parent.parent
    doc = _read_json(datasets_json)
    entries = _get_dataset_list(doc, datasets_json)
    datasets = [_parse_dataset(e, i, datasets_json, Path(root_dir)) for i, e in enumerate(entries)]
    _check_unique_names(datasets, datasets_json)
    return Datasets(datasets, datasets_json)


def _read_json(datasets_json):
    try:
        with open(datasets_json) as fh:
            return json.load(fh)
    except OSError as ex:
        raise FlairValidateDataError(f"can't read dataset catalog {datasets_json}") from ex
    except json.JSONDecodeError as ex:
        raise FlairValidateDataError(f"dataset catalog {datasets_json} is not valid JSON") from ex


def _get_dataset_list(doc, datasets_json):
    _check_keys(doc, ("datasets",), (), str(datasets_json))
    entries = doc["datasets"]
    if not isinstance(entries, list):
        raise FlairValidateDataError(f"{datasets_json}: field 'datasets' must be a list of dataset "
                                     f"objects, not {_type_name(entries)}")
    return entries


def _check_unique_names(datasets, datasets_json):
    counts = Counter(d.name for d in datasets)
    dups = sorted(name for name, cnt in counts.items() if cnt > 1)
    if dups:
        raise FlairValidateDataError(f"{datasets_json}: dataset name(s) {_quoted(dups)} used more "
                                     f"than once; each name must be unique")


_DATASET_REQUIRED = ("name", "description", "access", "organism", "assembly", "platform",
                     "validations", "files")
_DATASET_OPTIONAL = ("annotation", "library", "regions", "source")


def _parse_dataset(entry, index, datasets_json, root_dir):
    where = f"{datasets_json}: dataset {index}"
    if not isinstance(entry, dict):
        raise FlairValidateDataError(f"{where} must be an object, not {_type_name(entry)}")
    name = _parse_name(entry, where)
    where = f"{datasets_json}: dataset '{name}'"
    _check_keys(entry, _DATASET_REQUIRED, _DATASET_OPTIONAL, where)
    return _build_dataset(entry, name, where, root_dir)


def _build_dataset(entry, name, where, root_dir):
    return Dataset(name=name,
                   description=_get_str(entry, "description", where),
                   access=_get_enum(Access, entry, "access", where),
                   organism=_get_str(entry, "organism", where),
                   assembly=_get_str(entry, "assembly", where),
                   annotation=_get_opt_str(entry, "annotation", where),
                   platform=_get_enum(Platform, entry, "platform", where),
                   library=_get_opt_enum(Library, entry, "library", where),
                   regions=_get_str_list(entry, "regions", where, required=False),
                   validations=_parse_validations(entry, where),
                   source=_parse_source(entry, where),
                   files=_parse_files(entry["files"], where),
                   root_dir=root_dir)


def _parse_name(entry, where):
    name = _get_str(entry, "name", where)
    if not NAME_RE.match(name):
        raise FlairValidateDataError(f"{where}: name '{name}' must match {NAME_RE.pattern}; it "
                                     f"names directories and appears in sequence files, so it is "
                                     f"restricted to lower-case letters, digits and hyphens")
    return name


def _parse_validations(entry, where):
    validations = _get_enum_list(Validation, entry, "validations", where)
    if not validations:
        raise FlairValidateDataError(f"{where}: field 'validations' must name at least one of "
                                     f"{_quoted(Validation)}")
    return frozenset(validations)


def _parse_source(entry, where):
    source = _get_opt_str(entry, "source", where)
    if (source is not None) and source.startswith("/"):
        raise FlairValidateDataError(f"{where}: field 'source' is '{source}'; this file is "
                                     f"committed, so source holds an accession or a URL, never a "
                                     f"file path, which may name restricted storage")
    return source


def _parse_files(files, where):
    if not isinstance(files, dict):
        raise FlairValidateDataError(f"{where}: field 'files' must be an object mapping a role to "
                                     f"paths, not {_type_name(files)}")
    if len(files) == 0:
        raise FlairValidateDataError(f"{where}: field 'files' must name at least one file; roles "
                                     f"are {_quoted(FileRole)}")
    return MappingProxyType({_parse_role(r, where): _parse_paths(v, r, where)
                             for r, v in files.items()})


def _parse_role(role, where):
    try:
        return FileRole(role)
    except ValueError as ex:
        raise FlairValidateDataError(f"{where}: field 'files' has unknown role '{role}'; roles are "
                                     f"{_quoted(FileRole)}") from ex


def _parse_paths(value, role, where):
    values = value if isinstance(value, list) else [value]
    return tuple(_parse_path(v, role, where) for v in values)


def _parse_path(value, role, where):
    if not isinstance(value, str):
        raise FlairValidateDataError(f"{where}: path for role '{role}' must be a string, not "
                                     f"{_type_name(value)}")
    path = Path(value)
    _check_path_under_data(path, value, role, where)
    return path


def _check_path_under_data(path, value, role, where):
    if path.is_absolute() or (path.parts[0] != DATA_SUBDIR) or (".." in path.parts):
        raise FlairValidateDataError(f"{where}: path '{value}' for role '{role}' must be relative "
                                     f"and start with '{DATA_SUBDIR}/'; symlink the file under "
                                     f"{DATA_SUBDIR}/ rather than naming its location here, which "
                                     f"would commit the path of restricted data")


def _check_keys(obj, required, optional, where):
    missing = sorted(set(required) - set(obj))
    if missing:
        raise FlairValidateDataError(f"{where}: missing required field(s) {_quoted(missing)}")
    unknown = sorted(set(obj) - set(required) - set(optional))
    if unknown:
        raise FlairValidateDataError(f"{where}: unknown field(s) {_quoted(unknown)}; permitted "
                                     f"fields are {_quoted(tuple(required) + tuple(optional))}")


def _check_str(value, field, where):
    if not isinstance(value, str):
        raise FlairValidateDataError(f"{where}: field '{field}' must be a string, not "
                                     f"{_type_name(value)}")
    return value


def _get_str(obj, field, where):
    return _check_str(obj[field], field, where)


def _get_opt_str(obj, field, where):
    return _get_str(obj, field, where) if field in obj else None


def _get_str_list(obj, field, where, *, required=True):
    if field not in obj:
        if required:
            raise FlairValidateDataError(f"{where}: missing required field '{field}'")
        return ()
    values = obj[field]
    if not isinstance(values, list):
        raise FlairValidateDataError(f"{where}: field '{field}' must be a list of strings, not "
                                     f"{_type_name(values)}")
    return tuple(_check_str(v, field, where) for v in values)


def _check_enum(enum_cls, value, field, where):
    try:
        return enum_cls(value)
    except ValueError as ex:
        raise FlairValidateDataError(f"{where}: field '{field}' has '{value}', must be one of "
                                     f"{_quoted(enum_cls)}") from ex


def _get_enum(enum_cls, obj, field, where):
    return _check_enum(enum_cls, _get_str(obj, field, where), field, where)


def _get_opt_enum(enum_cls, obj, field, where):
    return _get_enum(enum_cls, obj, field, where) if field in obj else None


def _get_enum_list(enum_cls, obj, field, where):
    values = _get_str_list(obj, field, where)
    return tuple(_check_enum(enum_cls, v, field, where) for v in values)


def _type_name(value):
    return "null" if value is None else type(value).__name__


def _quoted(values):
    return ", ".join(f"'{v}'" for v in values)
