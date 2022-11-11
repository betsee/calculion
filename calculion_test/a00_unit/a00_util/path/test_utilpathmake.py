#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2021-2022 Ionovate.
# See "LICENSE" for further details.

'''
Project-wide **path factory** unit tests.

This submodule unit tests the public API of the private
:mod:`calculion._util.path.utilpathmake` submodule.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable test errors, avoid importing from
# package-specific submodules at module scope.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ TESTS                              }....................
def test_utilpathmake() -> None:
    '''
    Test the entirety of the
    :mod:`calculion._util.path.utilpathmake` submodule, particularly with
    respect to the :func:`calculion._util.path.utilpathmake.DirRelative` and
    :func:`calculion._util.path.utilpathmake.FileRelative` factory functions.
    '''

    # ..................{ IMPORTS                            }..................
    # Defer test-specific imports.
    from calculion.error import (
        CalculionDirException,
        CalculionFileException,
    )
    from calculion._util.path.utilpathmake import (
        DirRelative,
        FileRelative,
    )
    from pathlib import Path
    from pytest import raises

    # ..................{ LOCALS                             }..................
    # Path encapsulating the current module.
    TEST_MODULE_FILE = Path(__file__)

    # Path encapsulating the current module's package.
    TEST_PACKAGE_DIR = TEST_MODULE_FILE.parent

    # ..................{ PASS                               }..................
    # Path encapsulating the relative dirname of the parent directory of the
    # package containing this module. Since this function already validates
    # this path to be an existing directory, no additional validation is
    # warranted.
    DirRelative(TEST_PACKAGE_DIR, '../')

    # Path encapsulating the relative filename of the mandatory "__init__.py"
    # sibling submodule of this module. Since this function already validates
    # this path to be an existing file, no additional validation is warranted.
    FileRelative(TEST_PACKAGE_DIR, '__init__.py')

    # ..................{ FAIL                               }..................
    # Assert that attempting to create a path encapsulating the relative
    # dirname of a non-existent subpackage of this package (guaranteed *NOT* to
    # exist due to its ludicrous name) raises the expected exception.
    with raises(CalculionDirException):
        DirRelative(
            TEST_PACKAGE_DIR,
            'The works and ways of man, their death and birth,',
        )

    # Assert that attempting to create a path encapsulating the relative
    # filename of a non-existent sibling submodule of this module (guaranteed
    # *NOT* to exist due to its ludicrous name) raises the expected exception.
    with raises(CalculionFileException):
        FileRelative(
            TEST_PACKAGE_DIR, 'And that of him and all that his may be;')
