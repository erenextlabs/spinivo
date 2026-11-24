"""Read 1D proton FIDs from JEOL and Bruker acquisitions."""

from pathlib import Path

from .bruker import read_bruker
from .jeol import read_jeol
from .model import RawFID, UnsupportedInput

__all__ = ["RawFID", "UnsupportedInput", "read", "read_bruker", "read_jeol"]


def read(path):
    """Dispatch a JDF file or Bruker acquisition directory."""
    path = Path(path)
    if path.is_file() and path.suffix.lower() == ".jdf":
        return read_jeol(path)
    if path.is_dir() and any((path / name).is_file() for name in ("acqus", "fid", "ser")):
        return read_bruker(path)
    raise UnsupportedInput("not_found", path, "no JDF file or Bruker acquisition directory")
