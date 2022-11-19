#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2021-2022 Ionovate.
# See "LICENSE" for further details.

'''
Project-wide **project path** unit tests.

This submodule unit tests the public API of the private
:mod:`calculion._util.path.utilpathself` submodule.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable test errors, avoid importing from
# package-specific submodules at module scope.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ TESTS                              }....................
def test_utilpathself() -> None:
    '''
    Test the entirety of the
    :mod:`calculion._util.path.utilpathself` submodule.
    '''

    # Defer test-specific imports.
    from calculion._util.path.utilpathself import (
        get_package_dir,
        get_main_dir,
        get_main_readme_file,
        get_data_dir,
        get_data_png_dir,
        get_data_svg_dir,
        get_data_png_cell_network_schematic_0_file,
        get_data_png_cell_network_schematic_1_file,
        get_data_png_cell_network_schematic_2_file,
        get_data_png_cell_network_schematic_3_file,
        get_data_png_membrane_schematic_file,
    )

    # Assert each public getter published by the
    # "calculion._util.path.utilpathself" submodule implicitly behaves as
    # expected, thanks to type-checking performed by @beartype.
    get_package_dir()
    get_main_dir()
    get_main_readme_file()
    get_data_dir()
    get_data_png_dir()
    get_data_svg_dir()
    get_data_png_cell_network_schematic_0_file()
    get_data_png_cell_network_schematic_1_file()
    get_data_png_cell_network_schematic_2_file()
    get_data_png_cell_network_schematic_3_file()
    get_data_png_membrane_schematic_file()
