# Copyright 2026 Mark Diekhans
"""
Running the FLAIR modules.

One function per ``flair`` subcommand.  Each takes a FlairConfig saying which FLAIR
to run and the directory to run in, then that subcommand's required options by
name, then any other options as keywords.  A keyword becomes ``--name value``; True
becomes the bare flag, and None or False leaves the option off, so FLAIR's own
default applies::

    config = FlairConfig(flair_dir)
    flair_modules.transcriptome(config, stage_dir,
                                genome_aligned_bam=sample.alignments[0].path,
                                genome=reference.genome_fasta.path,
                                sample_name="wtc11-chr22-pb",
                                gtf=reference.annotation_gtf.path,
                                junction_tab=junctions.junctions_tab.path,
                                no_stringent=True)

Options are passed through rather than mirrored here, so this module does not go
stale when FLAIR gains an option.  FLAIR rejects what it does not know; run
``flair <subcommand> --help`` for the list.

Which FLAIR to run is the caller's to say, in the config.  Nothing here reads the
environment and nothing here has a default, so a run never depends on where it was
started from.  Selecting a version by git ref, rather than by a checkout the caller
already has, belongs here when it is built; see docs/design.md.

Each run writes three files into the directory it runs in:

``cmd``
    the command line as run, including the environment set for it, so the run can
    be repeated by hand.
``log``
    stdout and stderr together.
``time``
    output of ``/usr/bin/time --verbose``, which is where elapsed time and peak
    memory come from.
"""
from pathlib import Path

import pipettor

from flair_validate import FlairValidateError
from flair.pycbio.sys import fileOps

CMD_FILE = "cmd"
LOG_FILE = "log"
TIME_FILE = "time"


class FlairConfig:
    """Which FLAIR to run.  The caller builds one and passes it to every function
    here, so that what ran is a visible argument rather than something read from the
    environment at import time.

    ``flair_dir`` is a FLAIR checkout: the directory holding ``bin/flair`` and
    ``src``.  It is checked when the config is built, so a wrong path is reported
    before any work is done rather than at the first run."""

    def __init__(self, flair_dir):
        self.flair_dir = Path(flair_dir).absolute()
        self.prog = self.flair_dir / "bin/flair"
        self.src_dir = self.flair_dir / "src"
        self._check()

    def __repr__(self):
        return f"FlairConfig('{self.flair_dir}')"

    def _check(self):
        if not self.prog.exists():
            raise FlairValidateError(f"no flair program at {self.prog}; flair_dir must be a FLAIR "
                                     f"checkout, holding bin/flair and src")
        if not self.src_dir.is_dir():
            raise FlairValidateError(f"no src directory at {self.src_dir}; flair_dir must be a "
                                     f"FLAIR checkout, holding bin/flair and src")


def align(config, stage_dir, *, reads, **options):
    """Align reads to the genome.  Required: reads, one file or several.  See
    ``flair align --help`` for the rest."""
    return run_flair(config, stage_dir, "align", reads=reads, **options)


def transcriptome(config, stage_dir, *, genome_aligned_bam, genome, sample_name, **options):
    """Assemble a transcriptome from genome alignments.  Required:
    genome_aligned_bam, genome, sample_name.  The orthogonal evidence goes in as
    gtf, junction_tab or junction_bed.  See ``flair transcriptome --help``."""
    return run_flair(config, stage_dir, "transcriptome", genome_aligned_bam=genome_aligned_bam,
                     genome=genome, sample_name=sample_name, **options)


def quantify(config, stage_dir, *, manifest, genome, isoform_bed, **options):
    """Count reads per isoform.  Required: manifest, genome, isoform_bed.  See
    ``flair quantify --help``."""
    return run_flair(config, stage_dir, "quantify", manifest=manifest, genome=genome,
                     isoform_bed=isoform_bed, **options)


def combine(config, stage_dir, *, genome, **options):
    """Combine the transcriptomes of several samples into one.  Required: genome.
    See ``flair combine --help``."""
    return run_flair(config, stage_dir, "combine", genome=genome, **options)


def variantquant(config, stage_dir, **options):
    """Quantify isoforms by variant allele.  No required options.  See
    ``flair variantquant --help``."""
    return run_flair(config, stage_dir, "variantquant", **options)


def fusion(config, stage_dir, *, genome, gtf, genome_aligned_bam, sample_name, **options):
    """Call fusion transcripts.  Required: genome, gtf, genome_aligned_bam,
    sample_name.  See ``flair fusion --help``."""
    return run_flair(config, stage_dir, "fusion", genome=genome, gtf=gtf,
                     genome_aligned_bam=genome_aligned_bam, sample_name=sample_name, **options)


def diffexp(config, stage_dir, *, counts_matrix, output, **options):
    """Differential gene and isoform expression.  Required: counts_matrix, output.
    ``output`` is FLAIR's own output location, which is not stage_dir.  See
    ``flair diffexp --help``."""
    return run_flair(config, stage_dir, "diffexp", counts_matrix=counts_matrix, output=output, **options)


def diffsplice(config, stage_dir, *, isoform_bed, counts_matrix, output, **options):
    """Differential splicing events.  Required: isoform_bed, counts_matrix, output.
    ``output`` is FLAIR's own output location, which is not stage_dir.  See
    ``flair diffsplice --help``."""
    return run_flair(config, stage_dir, "diffsplice", isoform_bed=isoform_bed,
                     counts_matrix=counts_matrix, output=output, **options)


def alleles(config, stage_dir, *, tumor_bam, vcf, genome, **options):
    """Assign reads to alleles.  Required: tumor_bam, vcf, genome.  See
    ``flair alleles --help``."""
    return run_flair(config, stage_dir, "alleles", tumor_bam=tumor_bam, vcf=vcf, genome=genome, **options)


def isoalleles(config, stage_dir, *, allele_vcf, allele_read_map, genome, isoform_bed, iso_read_map,
               **options):
    """Isoforms per allele.  Required: allele_vcf, allele_read_map, genome,
    isoform_bed, iso_read_map.  See ``flair isoalleles --help``."""
    return run_flair(config, stage_dir, "isoalleles", allele_vcf=allele_vcf,
                     allele_read_map=allele_read_map, genome=genome, isoform_bed=isoform_bed,
                     iso_read_map=iso_read_map, **options)


def run_flair(config, stage_dir, subcommand, **options):
    """Run one flair subcommand in stage_dir, writing cmd, log and time there.  The
    subcommand functions above are the way to call this; use it directly only for a
    subcommand that has none."""
    stage_dir = Path(stage_dir).absolute()
    fileOps.ensureDir(stage_dir)
    cmd = flair_command(config, stage_dir, subcommand, options)
    _write_cmd(stage_dir, cmd)
    _run(stage_dir, cmd, subcommand)
    return stage_dir


def flair_command(config, stage_dir, subcommand, options):
    """the argv of a flair subcommand run: the environment it needs, the directory
    it runs in, and the timing wrapper, then flair itself.  It runs in stage_dir so
    that whatever the module writes to a relative path lands there."""
    return (["/usr/bin/env", f"--chdir={stage_dir}", f"PYTHONPATH={config.src_dir}",
             "/usr/bin/time", "--verbose", f"--output={TIME_FILE}",
             str(config.prog), subcommand]
            + option_args(options))


def option_args(options):
    "flair options as command line arguments, in a fixed order"
    args = []
    for name, value in sorted(options.items()):
        args.extend(_option_arg(name, value))
    return args


def _option_arg(name, value):
    flag = "--" + name
    if value is True:
        return [flag]
    if (value is None) or (value is False):
        return []
    if isinstance(value, (list, tuple)):
        return [flag] + [str(v) for v in value]
    return [flag, str(value)]


def _write_cmd(stage_dir, cmd):
    with fileOps.AtomicFileOpen(stage_dir / CMD_FILE) as fh:
        print(" ".join(cmd), file=fh)


def _run(stage_dir, cmd, subcommand):
    log_path = stage_dir / LOG_FILE
    try:
        with open(log_path, "w") as log_fh:
            pipettor.run([cmd], stdout=log_fh, stderr=log_fh)
    except pipettor.ProcessException as ex:
        raise FlairValidateError(f"flair {subcommand} failed in {stage_dir}; its output is in "
                                 f"{log_path}") from ex
