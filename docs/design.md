# Design

Draft.  Nothing here is built yet.  Open questions are at the end.

## Goal

- Run FLAIR on many combinations of input data, parameters and FLAIR version.
- Record accuracy statistics and run time and peak memory for every run.
- Compare any two runs on those statistics.
- Never run the same upstream work twice.
- First target is `flair transcriptome`.  Everything else follows the same pattern.

## Terminology

- **step**: a kind of work, such as `flair transcriptome`, `flair quantify`, SQANTI3.
  One Python function per step.
- **stage**: one run of one step with fixed inputs, parameters and FLAIR version.
  A node in the run tree.  Has a label and its own output directory.
- **run tree**: the one tree holding every stage ever run.  Not per experiment.
- **experiment**: a named question.  It declares the stages the question needs and
  the comparisons that answer it.  It owns no directories of its own beyond its
  report.
- **dataset**: a named set of input files, from `metadata/datasets.json`.
- **version**: a named FLAIR build.
- **metric**: one named number produced from a stage's output.
- **comparison**: the same metric from two stages, side by side.

## What varies

Three axes, chosen independently.

- **FLAIR version**: a git tag, branch or commit in a FLAIR checkout, a released
  version, or a working tree being developed.
- **Steps and their parameters**: which FLAIR modules run, in what order, with what
  options.  One branch point per parameter being tested.
- **Dataset**: which reads or alignments, which genome and annotation, and which
  orthogonal evidence (short-read junctions, CAGE peaks, poly(A) peaks).

Parameters sometimes have to differ between datasets to mean the same thing.  That
is a decision in the step function, not a table entry.  A step function that needs
evidence a dataset does not have raises an error naming the dataset and the role.

## The run tree

- One tree, under `results/runs/`, shared by every experiment.
- A stage's directory is its parent's directory plus its label, so the path reads as
  the history of how the stage was made:

```
results/runs/wtc11-chr22-pb/flair-2.0.0/transcriptome-default/sqanti/
results/runs/wtc11-chr22-pb/flair-2.0.0/transcriptome-no-sr-junctions/sqanti/
```

- The first two labels are the dataset and the FLAIR version.  Below that the labels
  are one per step run.
- **The path is the stage's identity.**  Declaring a stage at a path that already
  holds different work is an error naming both.  Nothing is keyed on a hash of the
  parameters, so a person reading a path knows what it is.
- Labels are given explicitly where the stage is declared, never derived from
  parameter values.  Derived labels stop being readable at about two parameters.
- Stages are declared by calling Python functions that build a branch:

```python
def transcriptome_default(tree, dataset, version):
    base = tree.root(dataset, version)
    return tree.transcriptome(base, "transcriptome-default")

def transcriptome_no_sr_junctions(tree, dataset, version):
    base = tree.root(dataset, version)
    return tree.transcriptome(base, "transcriptome-no-sr-junctions", junction_tab=None)
```

- Two experiments that both call `transcriptome_default` name the same path, so they
  get the same stage and it runs once.  Reuse across experiments is reuse of these
  functions, not a lookup.
- Declaring the same path twice with the same step, version, parameters and inputs is
  not an error.  It is how sharing works.

## Experiments

- An experiment is a question, written as a Python function.  It declares the stages
  it needs, then the comparisons that answer the question.

```python
def short_read_junctions(exp):
    "Do short-read junctions change transcriptome accuracy?"
    for dataset in exp.datasets_with(FileRole.short_read_junctions):
        with_sj = transcriptome_default(exp.tree, dataset, exp.version)
        without = transcriptome_no_sr_junctions(exp.tree, dataset, exp.version)
        exp.compare(with_sj.sqanti, without.sqanti, over=dataset.name)
```

- The experiment names what varies and what is held fixed.  Here the dataset and
  version are held fixed within each comparison and the junction input varies.
- Output is `results/experiments/<name>/`, holding `comparisons.tsv` and a report.
  No FLAIR output lives there; it points into `results/runs/`.
- Deleting an experiment does not delete any run.  Adding one usually adds no runs,
  because the branches it needs already exist.

## A stage directory

- The step's own output files, under the names FLAIR gives them.
- `cmd` the exact command line that was run.
- `log` stdout and stderr.
- `time` output of `/usr/bin/time --verbose`.
- `timing.tsv` elapsed seconds, CPU seconds and peak RSS, parsed from `time`.
- `stage.json` the stage identity: step, version with its git commit, parameters, and
  for each input file its path, size and modification time.
- `stats.tsv` for stages that produce statistics.  Long form, two columns, `metric`
  and `value`.  Long form so a new metric does not change a schema.

Re-running a stage whose `stage.json` does not match what is now declared is an
error.  It means a label was reused for different work.  The message says to pick a
new label or to remove the directory.

## Steps, in the order they are needed

- `align`: inserted only when a dataset has reads and no alignments.
- `transcriptome`: the first target.
- `sqanti`: SQANTI3 against the dataset's annotation; `*_classification.txt` and
  `*_junctions.txt` are the inputs to the metrics.
- `read-recovery`: how much of the read evidence the transcriptome represents.
  Carried over from `wtc11-chr22-eval/bin/evaluate_transcriptome_0226.py`, keeping
  only its read-relative metrics.  Needs the reads as BED12, which in the catalog
  today is only `wtc11-chr22-pb`.
- Later: `quantify`, `combine`, and the differential modules.

## Statistics and comparison

- Every statistics step writes `stats.tsv` in the same long form, whatever tool it
  wraps.  Parsing SQANTI output is one function; nothing downstream knows it was
  SQANTI.
- SQANTI3 measures agreement with the reference annotation, not with truth.  A
  transcript the annotation lacks is novel, which is either a discovery or an error,
  and SQANTI cannot tell them apart.  So an absolute SQANTI number is not a score.
  The difference between two runs that differ in one thing is the measurement, which
  is why comparison, not a threshold, is the unit of this suite.
- First cut of the SQANTI3 metrics: the structural classification relative to the
  reference, FSM, ISM, NIC and NNC, as counts and as percentages of the spliced
  isoforms.  More can be added without changing anything downstream.
- The two statistics answer different questions and both are needed.  SQANTI3 asks
  how the reported isoforms relate to the annotation.  `read-recovery` asks how much
  of the input the isoforms account for: what fraction of the genic regions the reads
  cover is covered by an isoform, and what fraction of the splice junction chains
  seen in reads appears in an isoform, exactly or as a subset.  A run can score well
  on one and badly on the other, and the pair is what says whether a parameter helped.
- `results/runs/stats.tsv` collects every stage: `stage`, `step`, `metric`, `value`.
- Comparisons are code.  A comparison names two stages and gets a row in
  `results/experiments/<name>/comparisons.tsv`: `metric`, `stage_a`, `stage_b`,
  `value_a`, `value_b`, `delta`.
- Because any stage can be compared with any other, a version comparison, a parameter
  comparison and a dataset comparison are the same operation.
- Speed and memory are metrics like any other, read from `timing.tsv`, so a parameter
  that costs an hour shows up next to the accuracy it bought.

## Snakemake driver

- The run tree is built in Python when the Snakefile is read.  By then every stage
  knows its inputs, its outputs and its command.
- One Snakemake rule per stage, generated in a loop with an explicit `name:`.  No
  wildcards and no path parsing, so `snakemake --list` prints the tree.
- `rule all` requires every leaf output, the collected `stats.tsv`, and every
  experiment report.
- Snakemake is the driver only.  Every rule body calls a function in
  `lib/flair_validate` or a program in `bin/`.  No analysis logic in the Snakefile.
- Execution is on one multi-core server.  `profiles/local/config.yaml` sets the core
  count and the per-rule thread counts; a cluster profile can be added beside it
  without touching a rule.  `wtc11-lrgasp-ont` is about 49 GB of reads and is what
  will force that.

## FLAIR versions

- There is no version catalog.  A version is a git ref, given on the command line.
- A ref is a tag, a branch, a commit, or the path of a working tree holding
  uncommitted work.
- `bin/flair-install` resolves the ref to a commit, checks it out, builds the
  environment, and records the resolved commit.
- Installs live under `build/flair/<label>/`, gitignored.
- The label is the ref as written when it is a tag or a branch.  A commit or a
  working tree needs a name given with it, since neither makes a readable directory
  and the label appears in every run path below it.
- `stage.json` records the resolved commit, and for a working tree a dirty marker.
- Running a stage sets `PATH` and `PYTHONPATH` to that install, the way
  `wtc11-chr22-eval/bin/run-transcriptome` does today.
- The suite's own code must not import from the FLAIR being evaluated.  Today
  `lib/flair_validate/__init__.py` imports `flair.pycbio` for `NoStackError`, which
  has to be severed before any stage runs.

## Code layout

```
bin/                       driver programs, no extension
lib/flair_validate/
  datasets.py              dataset catalog, exists
  versions.py              resolve a git ref to a commit and an install
  steps.py                 one function per step, builds the command line
  tree.py                  Stage, run tree, labels, paths, identity checking
  branches.py              named branch builders, shared between experiments
  runner.py                runs a command, writes cmd, log, time, timing.tsv, stage.json
  stats/sqanti.py          SQANTI output to stats.tsv
  stats/collect.py         run tree stats.tsv
  experiments/*.py         one module per experiment, each declaring stages and comparisons
metadata/datasets.json     dataset catalog
profiles/local/            snakemake profile for the multi-core server
Snakefile
results/runs/              every stage, gitignored
results/experiments/       comparisons and reports, gitignored
build/                     FLAIR installs, gitignored
```

## Code or data

- Data: facts about the world that the suite did not decide.  Which files make up a
  dataset.  That is `metadata/datasets.json`, and it is the only catalog.
- Code: everything the suite decides.  Which experiments exist, which parameters are
  worth varying, what a label means, how a dataset's quirks are handled.
- There is no configuration file naming runs, and no template expansion.

## Open questions

- `evaluate_transcriptome_0226.py` classifies isoforms as ISM by testing whether the
  junction chain string is a substring of a reference chain string.  That matches
  across chain boundaries and will call some isoforms ISM that are not.  Reimplement
  on junction tuples rather than carrying the bug forward.
- `read-recovery` needs reads as BED12 and only one dataset has that role.  Either
  derive it from the BAM in the step, or add the role to the datasets that need it.
