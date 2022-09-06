#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**CalculIon.**

For :pep:`8` compliance, this namespace exposes a subset of the metadata
constants provided by the :mod:`calculion.meta` submodule commonly inspected and
thus expected by external automation.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python modules and package-specific
# modules but *NOT* third-party dependencies, which if currently uninstalled
# will only be installed at some later time in the installation. Likewise, to
# avoid circular import dependencies, the top-level of this module should avoid
# importing package-specific modules where feasible.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ IMPORTS                            }....................
from calculion.meta import VERSION, VERSION_PARTS

# ....................{ GLOBALS                            }....................
# Declare PEP 8-compliant version constants expected by external automation.

__version__ = VERSION
'''
Human-readable application version as a ``.``-delimited string.

For :pep:`8` compliance, this specifier has the canonical name ``__version__``
rather than that of a typical global (e.g., ``VERSION_STR``).
'''


__version_info__ = VERSION_PARTS
'''
Machine-readable application version as a tuple of integers.

For :pep:`8` compliance, this specifier has the canonical name
``__version_info__`` rather than that of a typical global (e.g.,
``VERSION_PARTS``).
'''
