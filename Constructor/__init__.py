from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "0+unknown"

# Absolute imports
from Constructor.Constructor import Constructor, loadingBar, Logfile, Charge
from Constructor.Help import HelpDatabase, HelpSystem

__all__ = [
    "__version__",
    "Constructor",
    "loadingBar",
    "Logfile",
    "Charge",
    "HelpDatabase",
    "HelpSystem",
]
