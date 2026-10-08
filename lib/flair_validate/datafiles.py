# Copyright 2026 Mark Diekhans
"""
The kinds of file a catalog entry can name.

One class per kind.  The class documents that kind, gives the field name it is
written under in ``metadata/datasets.json``, says which file name suffixes it
accepts, and owns whatever is specific to it, such as where its index file sits.

Field names are Python identifiers, so a kind is named the same way in the JSON
and in code.

Every path in the catalog is relative to the repository root and starts with
``data/``.  An absolute path, a path outside ``data/``, or one holding ``..`` is
rejected.  That is what keeps the locations of restricted data out of git:
``data/`` is not committed and holds symlinks to the real files.

Files are not required to exist.  Most users have no access to the restricted
ones.  Use ``exists``, or ``Dataset.missing``, to find out.

Read this with ``python3 -m pydoc flair_validate.datafiles``.
"""
from pathlib import Path

from flair_validate import FlairValidateDataError
from flair_validate.jsonparse import type_name

DATA_SUBDIR = "data"


class DataFile:
    """One file named in the catalog.  Subclasses set ``key``, the field name,
    ``suffixes``, the file name endings it accepts, and ``multiple``, whether the
    field may hold a list of paths rather than one."""

    key = None
    suffixes = ()
    multiple = False

    def __init__(self, rel_path, root_dir, where):
        self.rel_path = Path(rel_path)
        self.root_dir = Path(root_dir)
        self._check_under_data(rel_path, where)
        self._check_suffix(where)

    def __str__(self):
        return str(self.rel_path)

    def __repr__(self):
        return f"{type(self).__name__}('{self.rel_path}')"

    @property
    def path(self):
        "location of the file, for opening"
        return self.root_dir / self.rel_path

    @property
    def exists(self):
        "does the file exist?  false is normal for restricted data"
        return self.path.exists()

    @property
    def compressed(self):
        "is this file gzip or bgzip compressed?"
        return self.rel_path.name.endswith(".gz")

    @classmethod
    def get(cls, entry, root_dir, where, *, required=True):
        """build this kind of file from its field of a catalog entry, returning one
        object, or a tuple when the kind allows several.  Missing and not required
        gives None, or an empty tuple for a kind allowing several."""
        if cls.key not in entry:
            if required:
                raise FlairValidateDataError(f"{where}: missing required field '{cls.key}'")
            return () if cls.multiple else None
        return cls._parse(entry[cls.key], root_dir, where)

    @classmethod
    def _parse(cls, value, root_dir, where):
        if not cls.multiple:
            return cls(_check_path_str(value, cls.key, where), root_dir, where)
        values = value if isinstance(value, list) else [value]
        return tuple(cls(_check_path_str(v, cls.key, where), root_dir, where) for v in values)

    def _check_under_data(self, rel_path, where):
        if self.rel_path.is_absolute() or (self.rel_path.parts[0] != DATA_SUBDIR) or (".." in self.rel_path.parts):
            raise FlairValidateDataError(f"{where}: path '{rel_path}' for '{self.key}' must be "
                                         f"relative and start with '{DATA_SUBDIR}/'; symlink the "
                                         f"file under {DATA_SUBDIR}/ rather than naming its "
                                         f"location here, which would commit the path of "
                                         f"restricted data")

    def _check_suffix(self, where):
        if not any(self.rel_path.name.endswith(sfx) for sfx in self.suffixes):
            raise FlairValidateDataError(f"{where}: file '{self.rel_path.name}' for '{self.key}' "
                                         f"must end with one of "
                                         f"{', '.join(repr(s) for s in self.suffixes)}; name the "
                                         f"symlink under {DATA_SUBDIR}/ for what the file is, or "
                                         f"use the field for the kind it really is")


class GenomeFasta(DataFile):
    """Genome sequence, FASTA, optionally bgzip compressed.  A ``.fai`` sits beside
    it, and a ``.gzi`` as well when it is compressed."""

    key = "genome_fasta"
    suffixes = (".fa", ".fasta", ".fa.gz", ".fasta.gz")

    @property
    def fai_path(self):
        "location of the FASTA index"
        return Path(str(self.path) + ".fai")

    @property
    def gzi_path(self):
        "location of the bgzip index, None when the file is not compressed"
        return Path(str(self.path) + ".gzi") if self.compressed else None


class AnnotationGtf(DataFile):
    """Reference annotation, GTF.  The transcript models FLAIR is given and the ones
    an accuracy statistic compares against."""

    key = "annotation_gtf"
    suffixes = (".gtf", ".gtf.gz")


class Reads(DataFile):
    """Long reads, FASTQ or FASTA, unaligned.  A sample may have several files, one
    per replicate or per flow cell."""

    key = "reads"
    suffixes = (".fastq", ".fq", ".fasta", ".fa",
                ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz")
    multiple = True


class GenomeBam(DataFile):
    """Reads aligned to the genome, BAM.  A ``.bai`` sits beside it.  A sample may
    have several, one per replicate."""

    key = "alignments"
    suffixes = (".bam",)
    multiple = True

    @property
    def index_path(self):
        "location of the BAM index"
        return Path(str(self.path) + ".bai")


class AlignmentsBed(DataFile):
    """The same alignments as BED12, which is what the read recovery statistic reads
    rather than the BAM."""

    key = "alignments_bed"
    suffixes = (".bed", ".bed.gz")
    multiple = True


class JunctionsTab(DataFile):
    """Splice junctions as a STAR ``SJ.out.tab``, which STAR writes when short reads
    are aligned with it.  Passed to FLAIR as ``--junction_tab``."""

    key = "junctions_tab"
    suffixes = (".tab", ".tsv")


class IntronsBed(DataFile):
    """Introns called by intron-prospector from a BAM, BED, with the number of
    supporting reads in the score column.  Passed to FLAIR as ``--junction_bed``,
    which reads that score and drops anything under ``--junction_support``.

    An alternative to a STAR ``SJ.out.tab``, not a different kind of evidence: the
    two differ in the tool and the format, not in what they say.  intron-prospector
    takes any alignments, short read or long read, and junctions are called the same
    way from either.

    FLAIR's help for ``--junction_bed`` calls this long-read output.  That is the
    parenthetical of a help string, not a restriction."""

    key = "introns_bed"
    suffixes = (".bed", ".bed.gz")


class CagePeaksBed(DataFile):
    """CAGE peaks, BED.  Evidence for transcription start sites."""

    key = "cage_peaks"
    suffixes = (".bed", ".bed.gz")


class PolyaPeaksBed(DataFile):
    """Poly(A) or QuantSeq peaks, BED.  Evidence for transcript end sites."""

    key = "polya_peaks"
    suffixes = (".bed", ".bed.gz")


class TruthTranscriptomeGtf(DataFile):
    """The transcriptome a simulated data set was generated from, so the one case
    where accuracy can be measured against truth rather than against an annotation."""

    key = "truth_transcriptome"
    suffixes = (".gtf", ".gtf.gz")


class RegionsTsv(DataFile):
    """Regions a run is restricted to, columns ``chrom``, ``start``, ``end``."""

    key = "regions_tsv"
    suffixes = (".tsv", ".tab", ".bed")


DATA_FILE_TYPES = (GenomeFasta, AnnotationGtf, Reads, GenomeBam, AlignmentsBed,
                   JunctionsTab, IntronsBed, CagePeaksBed, PolyaPeaksBed,
                   TruthTranscriptomeGtf, RegionsTsv)


def _check_path_str(value, key, where):
    if not isinstance(value, str):
        raise FlairValidateDataError(f"{where}: path for '{key}' must be a string, not "
                                     f"{type_name(value)}")
    return value
