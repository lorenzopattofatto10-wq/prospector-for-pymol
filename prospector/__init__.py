"""Prospector for PyMOL - find, show and explain protein structures by name.

Install through PyMOL's Plugin Manager (Plugin > Plugin Manager > Install New
Plugin > choose prospector.zip).  Then type, in the PyMOL command line:

    find hemoglobin
    grab p53
    info
"""
__version__ = "2.9.3"
__license__ = "MIT"   # see LICENSE


def __init_plugin__(app=None):
    """Called by PyMOL when the plugin is loaded."""
    from .pymol_layer import commands
    commands.register()
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt("Prospector for PyMOL", commands.open_panel)
    except Exception:
        pass   # no Qt menu (e.g. command-line PyMOL): the commands still work
    print("Prospector %s loaded. New here? Type  start   |   All commands:  info"
          % commands.VERSION)
