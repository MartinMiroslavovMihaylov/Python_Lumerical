# Configuration file for the Sphinx documentation builder.

import os
import sys
import shutil
import subprocess

# --- Add main repo to sys.path ---
CODE_DIR = os.environ.get("CODE_DIR")
if not CODE_DIR or not os.path.isdir(CODE_DIR):
    CODE_DIR = os.path.abspath("../Python_Lumerical")  # fallback

sys.path.insert(0, CODE_DIR)

# --- Expose code examples under 'examples-src' ---
def _mount_examples():
    src = os.path.join(CODE_DIR, "Examples")
    dst = os.path.abspath(os.path.join(os.path.dirname(__file__), "examples-src"))

    if not os.path.isdir(src):
        print(f"[conf.py] No code examples at {src}; skipping.", file=sys.stderr)
        return

    try:
        # Remove existing folder, file, or symlink safely
        if os.path.exists(dst):
            if os.path.islink(dst):
                os.unlink(dst)  # remove symlink
                print(f"[conf.py] Removed symlink {dst}")
            elif os.path.isdir(dst):
                shutil.rmtree(dst)  # remove folder
                print(f"[conf.py] Removed folder {dst}")
            else:
                os.remove(dst)  # remove file if somehow exists
                print(f"[conf.py] Removed file {dst}")

        # Try creating a symlink first
        try:
            os.symlink(src, dst, target_is_directory=True)
            print(f"[conf.py] Linked examples-src -> {src}")
        except (OSError, NotImplementedError):
            # Windows fallback: use junction
            try:
                subprocess.run(
                    ["cmd", "/c", "mklink", "/J", dst, src],
                    check=True,
                    shell=True,
                )
                print(f"[conf.py] Junction examples-src -> {src}")
            except Exception:
                # Last resort: copy folder
                shutil.copytree(src, dst)
                print(f"[conf.py] Copied examples to '{dst}' (symlink/junction unavailable)")

    except Exception as e:
        print(f"[conf.py] Failed to mount examples: {e}", file=sys.stderr)

_mount_examples()

# --- Project information ---
project = "Python-Lumerical Constructor Script"
author = "Martin Mihaylov"
release = "01.01.2024"

# --- General configuration ---
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_rtd_theme",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_typehints = "description"

# --- Mock heavy imports ---
autodoc_mock_imports = [
    "numpy",
    "pandas",
    "matplotlib",
    "os",
    "sys",
    "matlab",
    "imt",
]

# --- HTML output options ---
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "titles_only": False,
    "includehidden": True,
    "sticky_navigation": True,
}
html_static_path = []
