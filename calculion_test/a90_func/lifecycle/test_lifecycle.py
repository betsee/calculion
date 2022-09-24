#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**App lifecycle** (i.e., full process of spinning up, running, and shutting down
this Streamlit-based web app) integration tests.

This submodule integration tests the public API of the
:mod:`calculion` package.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable test errors, avoid importing from
# package-specific submodules at module scope.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ TESTS                              }....................
#FIXME: Consider posting a new feature request on the Streamlit issue tracker
#publicly suggesting that Streamlit publish an official "pytest-streamlit"
#plugin trivially enabling Streamlit apps to be tested under "pytest". That
#plugin would provide support (presumably in the form of one or more pytest
#fixtures) for:
#* Programmatically running any arbitrary Streamlit app.
#* A Streamlit mock API similar to that of the existing third-party
#  "streamlit_mock" package but hopefully more actively maintained against
#  upstream changes in the Streamlit API.
def test_app_lifecycle(monkeypatch) -> None:
    '''
    Integration test exercising the **maximally trivial app lifecycle** (i.e.,
    workflow spinning up, running, and then immediately shutting down this app).

    Parameters
    ----------
    monkeypatch : MonkeyPatch
        Builtin fixture object permitting object attributes to be safely
        modified for the duration of this test.
    '''

    # Defer test-specific imports.
    import sys
    from calculion import core
    from streamlit.web import cli

    # Absolute filename of the "calculion.core" submodule serving as the main
    # entry-point for this web app.
    core_filename = core.__file__

    #FIXME: Non-ideal. This assumes "streamlit" is in the current ${PATH}.
    #Frankly, this entire approach is absurd. Let's inspect the
    #"streamlit.web.cli" submodule a bit more closely to see if we can't tease
    #out a ${PATH}-agnostic approach. After all, Streamlit absolutely *MUST*
    #know where it itself lives, right?
    #
    #Alternately, we could always simply run something resembling:
    #    sys.argv = [python_filename, '-m', 'streamlit', 'run', main_filename]
    #FIXME: *OKAY.* Click is fundamentally insane. That's not necessarily
    #Streamlit's fault, of course. But see this ridiculous StackOverflow answer
    #on this exact subject of how to programmatically call an otherwise trivial
    #Click-decorated function:
    #    https://stackoverflow.com/a/48727139/2809027
    #
    #Thankfully, it turns out that you can just "call" a Click-decorated
    #function by passing that function a POSIX-compliant argument list rather
    #standard Python parameters. See also this StackOverflow answer:
    #    https://stackoverflow.com/a/54259095/2809027
    #
    #The approach below *MIGHT* actually be fine. Might. Let's test whether
    #Click requires that "streamlit" be in the current ${PATH} (e.g., by
    #temporarily moving our "/usr/bin/streamlit" binary aside). If it does, we
    #might still be able to work around this by calling the lower-level
    #streamlit.web.cli.main_run() rather than main() function. *shrug*

    # List of one or more arguments running this web app under the
    # "streamlit run" subcommand in the active Python process.
    STREAMLIT_RUN_ARGS = [
        'streamlit',
        'run',
        core_filename,

        # Sublist of zero or more "--"-prefixed arguments to be passed to this
        # "streamlit run" subcommand.
        #
        # Note this subcommand accepts *ALL* configuration settings accepted by
        # the ".streamlit/config.toml" file, specified as
        # "--{section_name}.{setting_name}={setting_value}".

        #FIXME: For unknown reasons, Streamlit appears to currently ignore this
        #and instead emit a non-fatal warning resembling:
        #    Warning: to view this Streamlit app on a browser, run it with the
        #    following command:
        #
        #    streamlit run /usr/lib/python3.10/site-packages/pytest/__main__.py [ARGUMENTS]

        # Prevent Streamlit from attempting to open a browser window, which also
        # has the harmful side effect of emitting a non-fatal warning
        # resembling:
        #    Warning: to view this Streamlit app on a browser, run it with the
        #    following command:
        '--server.headless', 'true',
    ]

    # Dynamically monkey-patch the list of one or more command-line arguments
    # with which the active Python process was invoked with a list instead
    # invoking the Streamlit-specific "run" subcommand against this entry-point.
    monkeypatch.setattr('sys.argv', STREAMLIT_RUN_ARGS)

    # Run this web app with these command-line arguments.
    cli.main()
