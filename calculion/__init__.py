#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Calculion.**

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
import beartype
def identity_decorator(func_or_cls, *args, **kwargs):
    return func_or_cls

beartype.beartype = identity_decorator

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

# ....................{ TODO                               }....................
#FIXME: [CONFIG] Define a new ".streamlit/config.toml" file resembling:
#
#    #FIXME: Also shift logging settings here from our "main" script, please.
#
#    [theme]
#    # Set Streamlit's built-in Dark Mode as the default theme for this web app.
#    # By default, Streamlit dynamically detects on app startup whether the
#    # current browser prefers a light or dark theme and conditionally sets that
#    # as the default theme for this web app. Since that yields a
#    # non-deterministic user experience (UX) across browsers, we force
#    # determinism through our preferred default theme: naturally, Dark.
#    base = "dark"

#FIXME: [SESSION] As time permits, implement most or all of the excellent advice
#at this blog article. Although the author focuses on session auto-save and
#auto-load to improve resiliency in the face of browser timeouts and refreshes
#(which is essential functionality, really), pretty much *ALL* of the advice
#here is awesome:
#    https://towardsdatascience.com/10-features-your-streamlit-ml-app-cant-do-without-implemented-f6b4f0d66d36
