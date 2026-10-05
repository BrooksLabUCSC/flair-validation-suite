# flair-validation-suite

Pipeline validating FLAIR module accuracy and speed across data sets, parameters and FLAIR
versions.  Global ~/.claude/CLAUDE.md rules apply; only project-specific rules are below.

## Pipeline shape is undecided

No directory layout, config format or stage names have been chosen.  Ask before creating any of
them.  metadata/ exists but its contents are not defined.

## Restricted data

- Brooks lab data lives under /private/groups/brookslab and is not publicly releasable.
- Symlink it into a gitignored directory.  Never copy it into the repo.
- Each dataset declares whether has.  The public subset runs without it.
- Never run a recursive search (find, grep -r, ls -R, rg) in this tree or in
  /private/groups/brookslab without asking first.  data/ is symlinks into brookslab, so a
  walk of this tree is a walk of a shared networked file system, slow for everyone on it.
- Test or list specific known paths instead.  A path that is wrong or missing is something
  to report, not to go hunting for.

## Selecting a FLAIR version
- undecided

## Python

- flair.pycbio.sys.cli for argparse plus logging, 
- flair.pycbio.sys.fileOps for atomic and compression-aware  writes, 
- flair.pycbio.tsv for TSV I/O.
- snake_case for new code.  pycbio's own API is camelCase and stays camelCase at the call site.

## Shell

- `#!/bin/bash` then `set -beEu -o pipefail`.
- Driver scripts in bin/ have no extension.

## Timing

/usr/bin/time --verbose for wall clock and peak memory; austin for profiles.

## Scratch

/notes/ is scratch for doc, gitignored, and never staged.
/work/ is scratch for code and experiments, gitignored, and never staged.
