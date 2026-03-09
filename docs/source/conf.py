import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))

project = "PolyGraph"
author = "Stephen Colodner"
copyright = "2026, Stephen Colodner"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",   # optional, harmless even with numpydoc
    "numpydoc",
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "furo"

autosummary_generate = True
autoclass_content = "class"
add_module_names = False

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False

master_doc = "API"