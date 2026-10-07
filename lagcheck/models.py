from dataclasses import dataclass, field
from typing import Any
import numpy as np


class InputError(ValueError):
    """Invalid or unsupported input; no scientific result was produced."""


class IntegrityError(InputError):
    """Bytes bound to a completed diagnostic or publication have changed."""


@dataclass
class SignalData:
    time: np.ndarray
    signal: np.ndarray
    metadata: dict = field(default_factory=dict)


@dataclass
class LaserData:
    lines: list[dict]
    metadata: dict = field(default_factory=dict)
    table: Any = None
