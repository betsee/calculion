#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Package metadata.**

This submodule exports global constants synopsizing this package -- including
versioning and dependencies.

Python Version
----------
For uniformity between this codebase and the ``setup.py`` setuptools script
importing this module, this module also validates the version of the active
Python 3 interpreter. An exception is raised if this version is insufficient.

As a tradeoff between backward compatibility, security, and maintainability,
this package strongly attempts to preserve compatibility with the first stable
release of the oldest version of CPython still under active development. Hence,
obsolete and insecure versions of CPython that have reached their official End
of Life (EoL) (e.g., Python 3.5) are explicitly unsupported.
'''

# ....................{ TODO                               }....................
#FIXME: *THIS FILE IS NOW OBSOLETE.* By Streamlit Cloud mandate, all project
#metadata now resides in the top-level "pyproject.toml" file. That said, various
#functionality throughout the codebase still expects this submodule to exist.
#Thankfully, we should *NOT* need to entirely eliminate this submodule. Instead,
#this submodule can continue to happily exist -- albeit as a surrogate of
#"pyproject.toml". How so? Simple. Refactor this submodule to either:
#* Defer to the "importlib.metadata" submodule, which internally defers to
#  "pyproject.toml". For example, here's how one would set this project's
#  version below:
#    from importlib.metadata import version
#    VERSION = version('calculion')
#  The issue here is that many globals declared below have *NO* corresponding
#  "importlib.metadata" API. So it goes. That said, there is *BASICALLY* no
#  alternative, because...
#* We can't defer to the "tomllib" submodule. This requires Python >= 3.11,
#  which is non-ideal. But that's not the *REAL* issue. The real issue is that
#  "tomllib" requires a TOML file to read. Makes sense. But "pyproject.toml"
#  *CANNOT* be assumed to exist when this project is installed in the standard
#  way to a "site-packages/" directory.

#FIXME: [QA] Prevalidate that all packages listed under
#"LIBS_RUNTIME_MANDATORY" are all importable at app startup *BEFORE* attempting
#to import from those packages. Doing so will expose whether any underlying
#C[++] dependencies for these packages are also installed. If MKL is *NOT* also
#installed, for example, importing "sparse_dot_mkl" raises this exception:
#    ImportError: Unable to load the MKL libraries through libmkl_rt. Try
#    setting $LD_LIBRARY_PATH. mkl_rt not found.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and package modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: To avoid polluting the public module namespace, external attributes
# should be locally imported at module scope *ONLY* under alternate private
# names (e.g., "from argparse import ArgumentParser as _ArgumentParser" rather
# than merely "from argparse import ArgumentParser").
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
import sys as _sys

# ....................{ METADATA                           }....................
NAME = 'Calculion'
'''
Human-readable package name.
'''

# ....................{ PYTHON ~ version                   }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: Changes to this section *MUST* be synchronized with:
# * Continuous integration test matrices, including:
#   * The top-level "tox.ini" file.
#   * The GitHub CI-specific ".github/workflows/python_test.yml" file.
# * Front-facing documentation (e.g., "README.rst", "doc/md/INSTALL.md").
#
# On bumping the minimum required version of Python, consider also documenting
# the justification for doing so in the "Python Version" section of this
# submodule's docstring above.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# A relatively recent version of Python is required due to:
# * PEP 585-compliant type hints (e.g., "list[str]"), now leveraged through the
#   codebase due both to convenience and PEP 585's timely deprecation of PEP
#   484-compliant type hints (e.g., "typing.list[str]") by 2026 -- which
#   renders PEP 484-compliant type hints genuinely dangerous in 2021.
# * "importlib.metadata", first introduced with Python 3.8.0.
#
# Note that:
# * Our mandatory runtime dependency on "pyvista" transitively requires "vtk".
#   For unknown reasons, "vtk" has yet to publish Python 3.10 wheels. Since
#   "vtk" is *ONLY* installable via binary wheels, this effectively means that
#   "vtk" and thus "pyvista" and thus this package currently *CANNOT* be
#   installed under Python 3.10. See also this unresolved VTK issue:
#       https://gitlab.kitware.com/vtk/vtk/-/issues/18335
PYTHON_VERSION_MIN = '3.10.0'
'''
Human-readable minimum version of Python required by this package as a
``.``-delimited string.

See Also
----------
"Python Version" section of this submodule's docstring for a detailed
justification of this constant's current value.
'''


def _convert_version_str_to_tuple(version_str: str) -> tuple:
    '''
    Convert the passed human-readable ``.``-delimited version string into a
    machine-readable version tuple of corresponding integers.
    '''
    assert isinstance(version_str, str), f'"{version_str}" not version string.'
    return tuple(int(version_part) for version_part in version_str.split('.'))


PYTHON_VERSION_MIN_PARTS = _convert_version_str_to_tuple(PYTHON_VERSION_MIN)
'''
Machine-readable minimum version of Python required by this package as a
tuple of integers.
'''


# Validate the version of the active Python interpreter *BEFORE* subsequent code
# possibly depending on this version. Since this version should be validated
# both at setuptools-based install time and post-install runtime *AND* since
# this module is imported sufficiently early by both, stash this validation here
# to avoid duplication of this logic and hence the hardcoded Python version.
#
# The "sys" module exposes three version-related constants for this purpose:
#
# * "hexversion", an integer intended to be specified in an obscure (albeit
#   both efficient and dependable) hexadecimal format: e.g.,
#    >>> sys.hexversion
#    33883376
#    >>> '%x' % sys.hexversion
#    '20504f0'
# * "version", a human-readable string: e.g.,
#    >>> sys.version
#    2.5.2 (r252:60911, Jul 31 2008, 17:28:52)
#    [GCC 4.2.3 (Ubuntu 4.2.3-2ubuntu7)]
# * "version_info", a tuple of three or more integers *OR* strings: e.g.,
#    >>> sys.version_info
#    (2, 5, 2, 'final', 0)
#
# For sanity, this package will *NEVER* conditionally depend upon the
# string-formatted release type of the current Python version exposed via the
# fourth element of the "version_info" tuple. Since the first three elements of
# that tuple are guaranteed to be integers *AND* since a comparable 3-tuple of
# integers is declared above, comparing the former and latter yield the
# simplest and most reliable Python version test.
#
# Note that the nearly decade-old and officially accepted PEP 345 proposed a
# new field "requires_python" configured via a key-value pair passed to the
# call to setup() in "setup.py" (e.g., "requires_python = ['>=2.2.1'],"), that
# field has yet to be integrated into either disutils or setuputils. Hence,
# that field is validated manually in the typical way.
if _sys.version_info[:3] < PYTHON_VERSION_MIN_PARTS:
    # Human-readable current version of Python. Ideally, "sys.version" would be
    # leveraged here instead; sadly, that string embeds significantly more than
    # merely a version and hence is inapplicable for real-world usage: e.g.,
    #     >>> import sys
    #     >>> sys.version
    #     '3.6.5 (default, Oct 28 2018, 19:51:39) \n[GCC 7.3.0]'
    _PYTHON_VERSION = '.'.join(
        str(version_part) for version_part in _sys.version_info[:3])

    # Die ignominiously.
    raise RuntimeError(
        f'{NAME} requires at least Python {PYTHON_VERSION_MIN}, but '
        f'the active interpreter only targets Python {_PYTHON_VERSION}. '
        f'We feel unbearable sadness for you.'
    )

# ....................{ METADATA ~ package                 }....................
PACKAGE_MAIN_NAME = 'calculion'
'''
Fully-qualified name of the top-level Python package containing this submodule.
'''


PACKAGE_TEST_NAME = f'{PACKAGE_MAIN_NAME}_test'
'''
Fully-qualified name of the top-level Python package exercising this project.
'''

# ....................{ METADATA ~ version                 }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: When modifying the current version of this package below, consider
# adhering to the Semantic Versioning schema. Specifically, the version should
# consist of three "."-delimited integers "{major}.{minor}.{patch}", where:
#
# * "{major}" specifies the major version, incremented only when either:
#   * Breaking backward compatibility in this package's public API.
#   * Implementing headline-worthy functionality (e.g., a GUI). Technically,
#     this condition breaks the Semantic Versioning schema, which stipulates
#     that *ONLY* changes breaking backward compatibility warrant major bumps.
#     But this is the real world. In the real world, significant improvements
#     are rewarded with significant version changes.
#   In either case, the minor and patch versions both reset to 0.
# * "{minor}" specifies the minor version, incremented only when implementing
#   customary functionality in a manner preserving such compatibility. In this
#   case, the patch version resets to 0.
# * "{patch}" specifies the patch version, incremented only when correcting
#   outstanding issues in a manner preserving such compatibility.
#
# When in doubt, increment only the minor version and reset the patch version.
# For further details, see http://semver.org.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: Synchronize changes to this value against the "version" setting in
# the root "buildozer.spec" file.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
VERSION = '0.1.0'
'''
Human-readable package version as a ``.``-delimited string.
'''


VERSION_PARTS = _convert_version_str_to_tuple(VERSION)
'''
Machine-readable package version as a tuple of integers.
'''
