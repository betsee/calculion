#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Project-wide **callable caching** (i.e., general-purpose memoization of
function and method calls) utilities.

This private submodule is *not* intended for importation by downstream callers.
'''

# ....................{ IMPORTS                            }....................

# ....................{ IMPORTS ~ publish                  }....................
# Publicize this private decorator for our own use. Although this constitutes a
# fundamental violation of privacy encapsulation:
# * We author the "beartype" package and thus guarantee backward compatibility.
# * We do *NOT* intend to publicize this private decorator, as doing so would
#   expose us to undue maintenance burden.
# * We do *NOT* intend to violate Don't Repeat Yourself (DRY) by
#   copying-and-pasting both:
#   * The entirety of this private "beartype" submodule.
#   * The entirety of unit tests exercising this private "beartype" submodule.
# This is actually the best possible solution, though it may not look like it.
from beartype._util.cache.utilcachecall import (
    callable_cached,
    property_cached,
)
