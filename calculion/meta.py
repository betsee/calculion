#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
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
NAME = 'CalculIon'
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
PYTHON_VERSION_MIN = '3.9.0'
'''
Human-readable minimum version of Python required by this package as a
``.``-delimited string.

See Also
----------
"Python Version" section of this submodule's docstring for a detailed
justification of this constant's current value.
'''


PYTHON_VERSION_MINOR_MAX = 12
'''
Maximum minor stable version of this major version of Python currently released
(e.g., ``5`` if Python 3.5 is the most recent stable version of Python 3.x).
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

# ....................{ METADATA ~ license                 }....................
LICENSE = 'MIT'
'''
Human-readable name of the license this package is licensed under.
'''

# ....................{ METADATA ~ package                 }....................
PACKAGE_NAME = 'calculion'
'''
Fully-qualified name of the top-level Python package containing this submodule.
'''


PACKAGE_TEST_NAME = f'{PACKAGE_NAME}_test'
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
VERSION = '0.0.1'
'''
Human-readable package version as a ``.``-delimited string.
'''


VERSION_PARTS = _convert_version_str_to_tuple(VERSION)
'''
Machine-readable package version as a tuple of integers.
'''

# ....................{ METADATA ~ synopsis                }....................
SYNOPSIS = (
    'CalulIon is an open-source cross-platform web-based simulator for '
    'single-cell computational problems in the field of bioelectricity.'
)
'''
Human-readable single-line synopsis of this package.

By PyPI design, this string must *not* span multiple lines or paragraphs.
'''

# ....................{ METADATA ~ authors                 }....................
AUTHOR_EMAIL = 'alexis.pietak@gmail.com'
'''
Email address of the principal corresponding author (i.e., the principal author
responding to public correspondence).
'''


AUTHORS = 'Alexis Pietak, Cecil Curry, et al.'
'''
Human-readable list of all principal authors of this package as a
comma-delimited string.

For brevity, this string *only* lists authors explicitly assigned copyrights.
For the list of all contributors regardless of copyright assignment or
attribution, see the top-level ``AUTHORS.md`` file.
'''


COPYRIGHT = '2022-2023 Alexis Pietak & Cecil Curry.'
'''
Legally binding copyright line excluding the license-specific prefix (e.g.,
``"Copyright (c)"``).

For brevity, this string *only* lists authors explicitly assigned copyrights.
For the list of all contributors regardless of copyright assignment or
attribution, see the top-level ``AUTHORS.md`` file.
'''

# ....................{ METADATA ~ urls                    }....................
URL_HOMEPAGE = 'https://github.com/betsee/calculion'
'''
URL of this package's homepage.
'''


URL_DOWNLOAD = f'{URL_HOMEPAGE}/archive/{VERSION}.tar.gz'
'''
URL of the source tarball for the current version of this package.

This URL assumes a tag whose name is ``v{VERSION}`` where ``{VERSION}`` is the
human-readable current version of this package (e.g., ``v0.4.0``) to exist.
Typically, no such tag exists for live versions of this package -- which
have yet to be stabilized and hence tagged. Hence, this URL is typically valid
*only* for previously released (rather than live) versions of this package.
'''


URL_ISSUES = f'{URL_HOMEPAGE}/issues'
'''
URL of this package's issue tracker.
'''



URL_RELEASES = f'{URL_HOMEPAGE}/releases'
'''
URL of this package's release list.
'''

# ....................{ METADATA ~ libs                    }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: Changes to this section *MUST* be synchronized with:
# * The "requirements" setting in the root "buildozer.spec" file.
# * Conda-specific package dependencies listed in the developer-specific
#   "conda/development/environment.yml" file.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: *AVOID ADDING GPL-LICENSED DEPENDENCIES.* Web backends are
# explicitly omitted from GPL-licensing requirements (notably, "infection" with
# viral open-sourcing), which is exactly why the Affero GPL (AGPL) exists.
# Nonetheless, this specific web backend would prefer to retain the option of
# transitively importing from one or more non-free closed-source Python
# packages (e.g., "mkl") despite not currently doing so. Since GPL-licensed
# dependencies are intrinsically incompatible with closed-source dependencies,
# our usage of the latter effectively precludes our usage of the former.
# There is *NO* reasonable way of circumventing this constraint, sadly.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Relatively recent versions of the standard scientific stack, selected merely
# because these versions are currently marked as stable on the source-based
# Gentoo Linux distribution that *EVERYONE* loves! :p
LIBS_RUNTIME_MANDATORY = (
    # QA stack. Dismantled, this is:
    # * beartype >= 0.11.0, the first release to support class decoration.
    'beartype >=0.11.0',

    # Science stack.
    'numpy >=1.22.0',
    'pandas >=1.5.0',
    'scipy >=1.7.0',
    'sympy >=1.9.0',

    # Web stack.
    'streamlit >=1.12.0',

    #FIXME: Uncomment as needed, please.
    # # 2D stack.
    # 'svgpathtools >=1.4.1',
    #
    # # 3D stack.
    # 'pyvista >=0.30.0',
)
'''
Mandatory runtime package dependencies as a tuple of :mod:`setuptools`-specific
requirements strings of the format ``{project_name}
{comparison1}{version1},...,{comparisonN}{versionN}``, where:

* ``{project_name}`` is a :mod:`setuptools`-specific project name (e.g.,
  ``numpy``, ``scipy``).
* ``{comparison1}`` and ``{comparisonN}`` are :mod:`setuptools`-specific
  version comparison operators. As well as standard mathematical comparison
  operators (e.g., ``==``, ``>=``, ``<``), :mod:`setuptools` also supports the
  PEP 440-compliant "compatible release" operator ``~=`` more commonly denoted
  by ``^`` in modern package managers (e.g., poetry, npm); this operator
  enables forward compatibility with all future versions of this dependency
  known *not* to break backward compatibility, but should only be applied to
  dependencies strictly following the semantic versioning contract.
* ``{version1}`` and ``{version1}`` are arbitrary version strings (e.g.,
  ``2020.2.16``, ``0.75a2``).
'''


#FIXME: Intentionally disabled. Optional dependencies have *NO* relevancy to
#mobile apps and only obfuscate already non-trivial portability complexities.
# LIBS_RUNTIME_OPTIONAL = ()
# '''
# Optional runtime package dependencies as a tuple of :mod:`setuptools`-specific
# requirements strings of the format ``{project_name}
# {comparison1}{version1},...,{comparisonN}{versionN}``.
# '''

# ....................{ METADATA ~ libs                    }....................
#FIXME: Uncomment as needed. We'll probably at least want this when reporting
#GitHub Actions-hosted CI test coverage, for example.
# LIBS_TESTTIME_MANDATORY_COVERAGE = (
#     'coverage >=5.5',
# )
# '''
# **Mandatory test-time coverage package dependencies** (i.e., dependencies
# required to measure test coverage for this package) as a tuple of
# :mod:`setuptools`-specific requirements strings of the format ``{project_name}
# {comparison1}{version1},...,{comparisonN}{versionN}``.
#
# See Also
# ----------
# :data:`LIBS_RUNTIME_OPTIONAL`
#     Further details.
# '''


LIBS_TESTTIME_MANDATORY = (
    # pytest should ideally remain the only hard dependency for testing on
    # local machines. While our testing regime optionally leverages third-party
    # frameworks and pytest plugins (e.g., "tox", "pytest-xdist"), these
    # dependencies are *NOT* required for simple testing.
    #
    # A relatively modern version of pytest is required.
    'pytest >=4.0.0',
)
'''
**Mandatory developer test-time package dependencies** (i.e., dependencies
required to test this package with :mod:`tox` as a developer at the command
line) as a tuple of :mod:`setuptools`-specific requirements strings of the
format ``{project_name} {comparison1}{version1},...,{comparisonN}{versionN}``.

See Also
----------
:data:`LIBS_RUNTIME_MANDATORY`
    Further details.
'''

# ....................{ METADATA ~ libs : doc              }....................
LIBS_DOCTIME_MANDATORY = (
    'sphinx >=4.4.0',
)
'''
**Mandatory developer documentation build-time package dependencies** (i.e.,
dependencies required to manually build documentation for this package as a
developer at the command line) as a tuple of :mod:`setuptools`-specific
requirements strings of the format ``{project_name}
{comparison1}{version1},...,{comparisonN}{versionN}``.

For flexibility, these dependencies are loosely relaxed to enable developers to
build with *any* versions satisfying at least the bare minimum. For the same
reason, optional documentation build-time package dependencies are omitted.
Since our documentation build system emits a non-fatal warning for each missing
optional dependency, omitting these optional dependencies here imposes no undue
hardships while improving usability.

See Also
----------
:data:`LIBS_RUNTIME_MANDATORY`
    Further details.
'''

# ....................{ METADATA ~ libs : dev              }....................
LIBS_DEVELOPER_MANDATORY = LIBS_TESTTIME_MANDATORY + LIBS_DOCTIME_MANDATORY
'''
**Mandatory developer package dependencies** (i.e., dependencies required to
develop and meaningfully contribute pull requests for this package) as a tuple
of :mod:`setuptools`-specific requirements strings of the format
``{project_name} {comparison1}{version1},...,{comparisonN}{versionN}``.

This tuple includes all mandatory test- and documentation build-time package
dependencies and is thus a convenient shorthand for those lower-level tuples.

See Also
----------
:data:`LIBS_RUNTIME_MANDATORY`
    Further details.
'''
