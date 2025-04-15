#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Calculion.**

For :pep:`8` compliance, this namespace exposes a subset of the metadata
constants provided by the :mod:`calculion.meta` submodule commonly inspected and
thus expected by external automation.
'''

# ....................{ TODO                               }....................
#FIXME: [SESSION] As time permits, implement most or all of the excellent advice
#at this blog article. Although the author focuses on session auto-save and
#auto-load to improve resiliency in the face of browser timeouts and refreshes
#(which is essential functionality, really), pretty much *ALL* of the advice
#here is awesome:
#    https://towardsdatascience.com/10-features-your-streamlit-ml-app-cant-do-without-implemented-f6b4f0d66d36

# ....................{ IMPORTS                            }....................
# Subject all subsequent imports to @beartype-based hybrid runtime-static
# type-checking *BEFORE* importing anything further.
# from beartype import BeartypeConf
from beartype.claw import beartype_this_package
beartype_this_package()

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
