"""CLI output should remain UTF-8 under legacy console encodings."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    ("args", "stream_name", "expected_text", "expected_code"),
    [
        (["run", "--help"], "stdout", "µs s⁻¹", 0),
        (
            ["run", "unused.json", "--output", "unused-out", "--reference-rate-ppm", "µ"],
            "stderr",
            "µ",
            2,
        ),
    ],
)
def test_cli_reconfigures_standard_streams_to_utf8(
    args: list[str], stream_name: str, expected_text: str, expected_code: int
) -> None:
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "ascii:strict"
    env["PYTHONUTF8"] = "0"

    result = subprocess.run(
        [sys.executable, "-X", "utf8=0", "-B", "-m", "lagcheck", *args],
        cwd=Path(__file__).resolve().parents[1],
        env=env,
        capture_output=True,
    )

    output = getattr(result, stream_name).decode("utf-8")
    assert result.returncode == expected_code, output
    assert expected_text in output
