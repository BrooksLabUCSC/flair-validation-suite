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
  A node in the experiment tree.  Has a label and its own output directory.
- **experiment**: a named tree of stages, built by one Python function.
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

## The experiment tree

- An experiment is Python code that builds stages, not a table that is expanded.
- Sharing of upstream work is by construction: a stage object is built once and
  passed to several downstream calls.  There is no key matching or hashing behind
  the user's back.

```python
def build(exp, dataset, version):
    base = exp.root(dataset)
    tc = exp.transcriptome(base, "default", version)
    tc_nosj = exp.transcriptome(base, "no-sr-junctions", version, junction_tab=None)
    for stage in (tc, tc_nosj):
        exp.sqanti(stage, "sqanti")
```

- Each stage call takes an explicit label.  Labels are not derived from parameter
  values, because derived labels become unreadable as soon as there are more than
  two of them.
- A stage's directory is its parent's directory plus its label, so the path is the
  stage's name and reads as the history of how it was made:

```
results/transcriptome-params/wtc11-chr22-pb/flair-2.0.0/default/sqanti/
results/transcriptome-params/wtc11-chr22-pb/flair-2.0.0/no-sr-junctions/sqanti/
```

- Convention, not mechanism: the first two labels are the dataset and the version.
  An experiment comparing two versions in one tree puts the version lower down.
- Two children of one parent with the same label is an error.
- Two stages anywhere in one experiment with the same step, version, parameters and
  input files is an error naming both paths.  That is an accidental duplicate; the
  fix is to build the stage once and use it twice.
- `results/<experiment>/index.tsv` lists every stage: path, step, version, label,
  parameters.  This is the map from a path to what it means.

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

Re-running a stage whose `stage.json` does not match what the experiment now says is
an error.  It means a label was reused for different work.  The message says to pick
a new label or to remove the directory.

## Steps, in the order they are needed

- `align`: inserted only when a dataset has reads and no alignments.
- `transcriptome`: the first target.
- `sqanti`: SQANTI3 against the dataset's annotation; `*_classification.txt` and
  `*_junctions.txt` are the inputs to the metrics.
- `evaluate`: the comparison against read evidence carried over from
  `wtc11-chr22-eval/bin/evaluate_transcriptome_0226.py`.
- Later: `quantify`, `combine`, and the differential modules.

## Statistics and comparison

- Every statistics step writes `stats.tsv` in the same long form, whatever tool it
  wraps.  Parsing SQANTI output is one function; nothing downstream knows it was
  SQANTI.
- `results/<experiment>/stats.tsv` collects every stage: `stage`, `step`, `metric`,
  `value`.
- Comparisons are also code.  A comparison names two stages and gets
  `results/<experiment>/comparisons.tsv`: `metric`, `stage_a`, `stage_b`, `value_a`,
  `value_b`, `delta`.
- Because any stage can be compared with any other, a version comparison, a parameter
  comparison and a dataset comparison are the same operation.
- Speed and memory are metrics like any other, read from `timing.tsv`, so a parameter
  that costs an hour shows up next to the accuracy it bought.

## Snakemake driver

- The experiment tree is built in Python when the Snakefile is read.  By then every
  stage knows its inputs, its outputs and its command.
- One Snakemake rule per stage, generated in a loop with an explicit `name:`.  No
  wildcards and no path parsing, so `snakemake --list` prints the tree.
- `rule all` requires every leaf output plus the collected `stats.tsv`.
- Snakemake is the driver only.  Every rule body calls a function in
  `lib/flair_validate` or a program in `bin/`.  No analysis logic in the Snakefile.
- Big datasets need a cluster; one Snakemake job per stage, submitted through a
  profile.  `wtc11-lrgasp-ont` is about 49 GB of reads and sets the resource request.

## FLAIR versions

- A version is named, and resolves to a directory holding a FLAIR checkout plus the
  environment to run it.
- `bin/flair-install` materializes one: clone or fetch, check out the tag or commit,
  build the environment, record the resolved commit.
- Installs live under `build/flair/<version>/`, gitignored.
- A version may also point at an existing working tree, for testing uncommitted work.
  Its recorded commit is then the commit plus a dirty marker.
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
  versions.py              FLAIR version catalog and installs
  steps.py                 one function per step, builds the command line
  tree.py                  Stage, Experiment, labels, paths, duplicate detection
  runner.py                runs a command, writes cmd, log, time, timing.tsv, stage.json
  stats/sqanti.py          SQANTI output to stats.tsv
  stats/collect.py         per-experiment stats.tsv and comparisons.tsv
  experiments/*.py         one module per experiment, each defining build()
metadata/datasets.json     dataset catalog
metadata/flair-versions.json   version catalog
Snakefile
results/                   gitignored
build/                     gitignored
```

## Code or data

- Data: facts about the world that the suite did not decide.  Which files make up a
  dataset, where a FLAIR version comes from.  These are the two JSON catalogs.
- Code: everything the suite decides.  Which experiments exist, which parameters are
  worth varying, what a label means, how a dataset's quirks are handled.
- There is no configuration file naming runs, and no template expansion.

## Open questions

- Which cluster and which Snakemake profile.
- Whether `metadata/flair-versions.json` is worth having, or whether a version is
  just a git ref given on the command line and recorded in `stage.json`.
- Whether SQANTI3 is the accuracy tool, and which of its metrics matter.
  `evaluate_transcriptome_0226.py` covers some of the same ground.
- How to re-use a stage across experiments, not just within one.  Two experiments
  that start from the same dataset and version repeat the whole upstream branch.
- What `index.tsv` should hold so that a person can find a run months later.
