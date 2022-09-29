#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Non-reusable **Python app entry-point wrapper** (i.e., submodule
unconditionally running this app as an intentional side effect of importation
by the active Python interpreter itself and thus *not* intended to be imported
from anywhere within this codebase).

This submodule is implicitly run by the active Python interpreter when the
fully-qualified name of the top-level package directly containing this
submodule is passed as the value of the ``--m`` option to this interpreter on
the command line (e.g., ``python3 -m ionyouapp``).

Caveats
----------
Streamlit does *not* support running high-level submodules contained in Python
packages as Streamlit-based web apps. Currently, Streamlit *only* supports
running low-level Python scripts agnostic of Python packaging. The latter
approach reduces the barrier to entry (which is Streamlit's entire raison
d'être) while significantly complicating real-world logistics like packaging,
distribution, and execution.

Thankfully, intrepid StackOverflowers have already reverse-engineered
Streamlit's invocation mechanism to easily "unlock" support for the former
approach as well. This submodule is architected ala this StackOverflow answer:
    https://stackoverflow.com/a/62775219/2809027
'''

# ....................{ TODO                               }....................
#FIXME: Unsurprisingly, we ended up *NOT* using this file. Streamlit is
#incredibly particular about startup execution. We finally relented and now
#simply launch this app in the standard (albeit non-ideal) Streamlit way.

#FIXME: Consider submitting an answer to this StackOverflow thread:
#    https://stackoverflow.com/a/62775219/2809027
#
#The approach performed below is genuinely novel and useful enough to see a
#wider audience. The currently top-upvoted answer assumes that the script to be
#run resides in the current working directory -- inviting catastrophe in the
#general deployment case. Instead, the approach below dynamically obtains the
#absolute filename of an entry-point submodule residing in a full-blown Python
#package to be run. \o/

#FIXME: An alternative to the approach employed below would be to manually
#inject the top-level "calculion" directory into "sys.path". Of course, someone
#has already thought of that as well. See also:
#    https://www.isticktoit.net/?p=2499

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: Avoid importing anything at module scope *EXCEPT* from official
# Python modules in the standard library guaranteed to exist. Subsequent logic
# in the _main() function called below validates third-party runtime
# dependencies of this package to be safely importable, Before performing that
# validation, *NO* other modules are safely importable from.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
import sys
from calculion import main
from streamlit.web import cli

# ....................{ MAIN                               }....................
# If this submodule is directly invoked by Python as the main module to be run
# (e.g., by being passed as the argument to the "-m" option), run this
# Streamlit-based web app.
if __name__ == '__main__':
    # Absolute filename of the "calculion.main" submodule serving as the
    # entry-point for this web app.
    core_filename = main.__file__

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

    # Dynamically monkey-patch the list of one or more command-line arguments
    # with which the active Python process was invoked with a list instead
    # invoking the Streamlit-specific "run" subcommand against this entry-point.
    sys.argv = ['streamlit', 'run', core_filename]

    # Run this web app, returning the exit code returned by that function as the
    # exit code of the active Python process.
    sys.exit(cli.main())
