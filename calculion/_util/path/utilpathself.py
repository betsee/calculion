#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Project paths** (i.e., :class:`pathlib.Path` instances encapsulating
project-specific paths relative to the current directory containing this
project).
'''

# ....................{ IMPORTS                            }....................
from calculion.meta import PACKAGE_TEST_NAME
from calculion._util.cache.utilcachecall import callable_cached
from calculion._util.path.utilpathmake import (
    DirRelative,
    FileRelative,
)
from pathlib import Path

# ....................{ GETTERS ~ package                  }....................
@callable_cached
def get_package_dir() -> Path:
    '''
    :mod:`Path` encapsulating the absolute dirname of the **top-level package**
    (i.e., directory providing this package's top-level package containing at
    least an ``__init__.py`` file) if found *or* raise an exception otherwise.
    '''

    # Path encapsulating the current module.
    MODULE_FILE = Path(__file__)

    # Path encapsulating the current module's package.
    MODULE_PACKAGE_DIR = MODULE_FILE.parent

    # Path encapsulating the dirname of this package's directory relative to the
    # dirname of the subpackage defining the current module.
    PACKAGE_DIR = DirRelative(MODULE_PACKAGE_DIR, '../../')

    # If this package's directory either does not contain our package-specific
    # "meta" submodule *OR* does but that path is not a file, raise an
    # exception. This basic sanity check improves the likelihood that this
    # package directory is what we assume it is.
    #
    # Note that we intentionally avoid testing paths *NOT* bundled with release
    # tarballs (e.g., a root ".git/" directory), as doing so would prevent
    # external users and tooling from running tests from release tarballs.
    FileRelative(PACKAGE_DIR, 'meta.py')

    # Return this path.
    return PACKAGE_DIR

# ....................{ GETTERS ~ main                     }....................
@callable_cached
def get_main_dir() -> Path:
    '''
    :mod:`Path` encapsulating the absolute dirname of the **root project
    directory** (i.e., directory containing both a ``.git/`` subdirectory and a
    subdirectory providing this project's package) if found *or* raise an
    exception otherwise.
    '''
    # print(f'current module paths: {__package__} [{__file__}]')

    # Path encapsulating the dirname of this project's directory relative to the
    # dirname of the top-level package defining this project.
    MAIN_DIR = DirRelative(get_package_dir(), '../')

    # If this project's directory either does not contain a test-specific
    # subdirectory *OR* does but that path is not a directory, raise an
    # exception. This basic sanity check improves the likelihood that this
    # project directory is what we assume it is.
    #
    # Note that we intentionally avoid testing paths *NOT* bundled with release
    # tarballs (e.g., a root ".git/" directory), as doing so would prevent
    # external users and tooling from running tests from release tarballs.
    DirRelative(MAIN_DIR, PACKAGE_TEST_NAME)

    # Return this path.
    return MAIN_DIR


@callable_cached
def get_main_readme_file() -> Path:
    '''
    :mod:`Path` encapsulating the absolute filename of the **project readme
    file** (i.e., this project's front-facing ``README.rst`` file) if found
    *or* raise an exception otherwise.

    Note that the :meth:`Path.read_text` method of this object trivially yields
    the decoded plaintext contents of this file as a string.
    '''

    # Perverse pomposity!
    return FileRelative(get_main_dir(), 'README.rst')

# ....................{ GETTERS ~ data : dir               }....................
@callable_cached
def get_data_dir() -> Path:
    '''
    :mod:`Path` encapsulating the absolute dirname of the **project-wide data
    subdirectory** (i.e., directory providing supplementary non-Python paths
    required throughout this package and thus *not* containing an
    ``__init__.py`` file) if found *or* raise an exception otherwise.
    '''

    # Obverse obviation!
    return DirRelative(get_package_dir(), 'data')


@callable_cached
def get_data_png_dir() -> Path:
    '''
    :mod:`Path` encapsulating the absolute dirname of the **project-wide
    portable network graphics (PNG) subdirectory** (i.e., directory containing
    ``.png``-suffixed files describing lossless but *not* scalable images) if
    found *or* raise an exception otherwise.
    '''
    return DirRelative(get_data_dir(), 'png')


@callable_cached
def get_data_svg_dir() -> Path:
    '''
    :mod:`Path` encapsulating the absolute dirname of the **project-wide
    scalable vector graphics (SVG) subdirectory** (i.e., directory containing
    ``.svg``-suffixed files describing losslessly scalable images) if found *or*
    raise an exception otherwise.
    '''

    # Perverse pomposity!
    return DirRelative(get_data_dir(), 'svg')

@callable_cached
def get_data_html_dir() -> Path:
    '''
    :mod:`Path` encapsulating the absolute dirname of the **project-wide
    portable html subdirectory** (i.e., directory containing
    ``.html``-suffixed files describing lossless but *not* scalable images) if
    found *or* raise an exception otherwise.
    '''
    return DirRelative(get_data_dir(), 'html')

# ....................{ GETTERS ~ data : file : png        }....................
#FIXME: Docstring us up, please.
@callable_cached
def get_data_png_cell_network_schematic_0_file() -> Path:

    return FileRelative(
        get_data_png_dir(), 'CellNetworkSchematic_3C_BaseLayer.png')


#FIXME: Docstring us up, please.
@callable_cached
def get_data_png_cell_network_schematic_1_file() -> Path:

    return FileRelative(get_data_png_dir(), 'CellNetworkSchematic_3C_1.png')


#FIXME: Docstring us up, please.
@callable_cached
def get_data_png_cell_network_schematic_2_file() -> Path:

    return FileRelative(get_data_png_dir(), 'CellNetworkSchematic_3C_2.png')


#FIXME: Docstring us up, please.
@callable_cached
def get_data_png_cell_network_schematic_3_file() -> Path:

    return FileRelative(get_data_png_dir(), 'CellNetworkSchematic_3C_3.png')

@callable_cached
def get_data_png_cell_network_schematic_4_file() -> Path:

    return FileRelative(get_data_png_dir(), 'CellNetworkSchematic_3D_BaseLayer.png')

@callable_cached
def get_data_png_cell_network_schematic_5_file() -> Path:

    return FileRelative(get_data_png_dir(), 'CellNetworkSchematic_3D_1.png')

@callable_cached
def get_data_png_cell_network_schematic_6_file() -> Path:

    return FileRelative(get_data_png_dir(), 'CellNetworkSchematic_3D_2.png')

@callable_cached
def get_data_png_cell_network_schematic_7_file() -> Path:

    return FileRelative(get_data_png_dir(), 'CellNetworkSchematic_3D_3.png')

@callable_cached
def get_bioe_network_file() -> Path:

    return FileRelative(get_data_svg_dir(), 'BioENet.svg')




#FIXME: Docstring us up, please.
@callable_cached
def get_data_png_membrane_schematic_file() -> Path:

    return FileRelative(get_data_png_dir(), 'MembraneSchematic_2.png')

@callable_cached
def get_data_png_banner_file() -> Path:

    return FileRelative(get_data_png_dir(), 'calculion_logo_lion_banner_i.png')


# @callable_cached
# def get_data_html_test_file() -> Path:
#
#     return FileRelative(get_data_html_dir(), 'test_1.html')

# ....................{ GETTERS ~ data : file : svg        }....................
# @callable_cached
# def get_data_svg_cell_network_schematic_file() -> Path:
#     '''
#     :mod:`Path` encapsulating the absolute filename of the **project-wide**
#     ``CellNetworkSchematic_3.svg`` file** if found *or* raise an exception
#     otherwise.
#     '''
#
#     # Terrifying terseness!
#     return FileRelative(get_data_svg_dir(), 'CellNetworkSchematic_3.svg')
#
#
# @callable_cached
# def get_data_svg_cell_network_schematic_b_file() -> Path:
#     '''
#     :mod:`Path` encapsulating the absolute filename of the **project-wide**
#     ``CellNetworkSchematic_3B.svg`` file** if found *or* raise an exception
#     otherwise.
#     '''
#
#     # Transverse transgression!
#     return FileRelative(get_data_svg_dir(), 'CellNetworkSchematic_3B.svg')
#
#
# @callable_cached
# def get_data_svg_membrane_schematic_file() -> Path:
#     '''
#     :mod:`Path` encapsulating the absolute filename of the **project-wide**
#     ``MembraneSchematic_2.svg`` file** if found *or* raise an exception
#     otherwise.
#     '''
#
#     # Terrifying terseness!
#     return FileRelative(get_data_svg_dir(), 'MembraneSchematic_2.svg')
