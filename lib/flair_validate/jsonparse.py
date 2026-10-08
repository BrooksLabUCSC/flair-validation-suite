# Copyright 2026 Mark Diekhans
"""
Checked access to fields of a JSON object read from a metadata file.

Every function takes ``where``, a string naming what is being parsed, and raises
``FlairValidateDataError`` with a message saying what is wrong and what is
permitted.  The caller never has to check a type itself.
"""
from flair_validate import FlairValidateDataError


def check_keys(obj, required, optional, where):
    "every required field present and no unknown field"
    missing = sorted(set(required) - set(obj))
    if missing:
        raise FlairValidateDataError(f"{where}: missing required field(s) {quoted(missing)}")
    unknown = sorted(set(obj) - set(required) - set(optional))
    if unknown:
        raise FlairValidateDataError(f"{where}: unknown field(s) {quoted(unknown)}; permitted "
                                     f"fields are {quoted(tuple(required) + tuple(optional))}")


def check_str(value, field, where):
    "value of field, which must be a string"
    if not isinstance(value, str):
        raise FlairValidateDataError(f"{where}: field '{field}' must be a string, not "
                                     f"{type_name(value)}")
    return value


def get_str(obj, field, where):
    "required string field"
    return check_str(obj[field], field, where)


def get_opt_str(obj, field, where):
    "string field, None when absent"
    return get_str(obj, field, where) if field in obj else None


def get_str_list(obj, field, where, *, required=True):
    "list of strings, empty when absent and not required"
    if field not in obj:
        if required:
            raise FlairValidateDataError(f"{where}: missing required field '{field}'")
        return ()
    values = obj[field]
    if not isinstance(values, list):
        raise FlairValidateDataError(f"{where}: field '{field}' must be a list of strings, not "
                                     f"{type_name(values)}")
    return tuple(check_str(v, field, where) for v in values)


def check_enum(enum_cls, value, field, where):
    "value of field, which must name a member of enum_cls"
    try:
        return enum_cls(value)
    except ValueError as ex:
        raise FlairValidateDataError(f"{where}: field '{field}' has '{value}', must be one of "
                                     f"{quoted(enum_cls)}") from ex


def get_enum(enum_cls, obj, field, where):
    "required field naming a member of enum_cls"
    return check_enum(enum_cls, get_str(obj, field, where), field, where)


def get_opt_enum(enum_cls, obj, field, where):
    "field naming a member of enum_cls, None when absent"
    return get_enum(enum_cls, obj, field, where) if field in obj else None


def type_name(value):
    "JSON name of a value's type, for an error message"
    return "null" if value is None else type(value).__name__


def quoted(values):
    "values as a quoted, comma separated list, for an error message"
    return ", ".join(f"'{v}'" for v in values)
