# Python_Lumerical/__init__.py

from importlib.metadata import PackageNotFoundError, version

# --- Version ---
try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "0+unknown"  # fallback for source installs

# --- Import classes/functions from Constructor module ---
from Constructor.Constructor import (
    Constructor,
    loadingBar,
    Help,
    Logfile,
    Charge,
)

# --- Exported names ---
__all__ = [
    "__version__",
    "Constructor",
    "loadingBar",
    "Help",
    "Logfile",
    "Charge",
]