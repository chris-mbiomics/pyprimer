"""Conversion between JSON payloads and the library's dataclasses."""

from __future__ import annotations

import dataclasses
import enum
import types
import typing
from functools import cache
from typing import Any, TypeVar

T = TypeVar("T")


class SchemaError(ValueError):
    """Raised when a request payload does not match the expected structure."""


def build(cls: type[T], data: Any, path: str = "") -> T:
    """Construct dataclass ``cls`` from a JSON object, validating field types.

    Missing fields take their defaults; unknown fields are rejected.
    """
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise SchemaError(f"{path or cls.__name__}: expected an object")
    hints = _hints(typing.cast(Any, cls))
    names = {f.name for f in dataclasses.fields(cls)}  # type: ignore[arg-type]
    unknown = set(data) - names
    if unknown:
        raise SchemaError(f"{path or cls.__name__}: unknown field(s) {', '.join(sorted(unknown))}")
    kwargs = {name: _coerce(hints[name], value, f"{path}{name}") for name, value in data.items()}
    try:
        return cls(**kwargs)
    except (TypeError, ValueError) as exc:
        raise SchemaError(f"{path or cls.__name__}: {exc}") from exc


@cache
def _hints(cls: Any) -> dict[str, Any]:
    return typing.get_type_hints(cls)


def _coerce(hint: Any, value: Any, path: str) -> Any:
    origin = typing.get_origin(hint)
    if origin in (typing.Union, types.UnionType):
        options = typing.get_args(hint)
        if value is None:
            if type(None) in options:
                return None
            raise SchemaError(f"{path}: a value is required")
        last_error: SchemaError | None = None
        for option in options:
            if option is type(None):
                continue
            try:
                return _coerce(option, value, path)
            except SchemaError as exc:
                last_error = exc
        assert last_error is not None
        raise last_error
    if origin in (list, tuple):
        if not isinstance(value, list):
            raise SchemaError(f"{path}: expected a list")
        (item,) = typing.get_args(hint)[:1]
        items = [_coerce(item, v, f"{path}[{i}]") for i, v in enumerate(value)]
        return items if origin is list else tuple(items)
    if dataclasses.is_dataclass(hint):
        return build(typing.cast(type, hint), value, f"{path}.")
    if isinstance(hint, type) and issubclass(hint, enum.Enum):
        try:
            return hint(value)
        except ValueError as exc:
            raise SchemaError(f"{path}: invalid value {value!r}") from exc
    if hint is bool:
        if not isinstance(value, bool):
            raise SchemaError(f"{path}: expected true or false")
        return value
    if hint is int:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value != int(value):
            raise SchemaError(f"{path}: expected a whole number")
        return int(value)
    if hint is float:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise SchemaError(f"{path}: expected a number")
        return float(value)
    if hint is str:
        if not isinstance(value, str):
            raise SchemaError(f"{path}: expected text")
        return value
    raise SchemaError(f"{path}: unsupported field type")


def dump(obj: Any) -> Any:
    """Convert dataclasses (recursively) and enums to JSON-compatible values."""
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {f.name: dump(getattr(obj, f.name)) for f in dataclasses.fields(obj)}
    if isinstance(obj, enum.Enum):
        return obj.value
    if isinstance(obj, (list, tuple)):
        return [dump(v) for v in obj]
    if isinstance(obj, dict):
        return {k: dump(v) for k, v in obj.items()}
    return obj
