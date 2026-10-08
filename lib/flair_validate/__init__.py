# Copyright 2026 Mark Diekhans
"""FLAIR validation suite library."""
from pathlib import Path

from flair.pycbio import NoStackError

REPO_ROOT = Path(__file__).parents[2]


class FlairValidateError(Exception):
    """General error condition in the FLAIR validation suite."""
    pass


class FlairValidateDataError(FlairValidateError, NoStackError):
    """Error in metadata or data read by the validation suite.  Reported to the
    user without a stack trace by pycbio.sys.cli.ErrorHandler."""
    pass
