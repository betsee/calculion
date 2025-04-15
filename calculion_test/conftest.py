#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


'''
:mod:`pytest` **global test configuration** (i.e., early-time configuration
guaranteed to be run by :mod:`pytest` *after* passed command-line arguments are
parsed).

:mod:`pytest` implicitly imports *all* functionality defined by this module
into *all* submodules of this subpackage.
'''

# ....................{ TODO                               }....................
#FIXME: *HOOBOY.* Testing Streamlit apps is "excitingly" non-trivial, as
#Streamlit has yet to implement a proper solution. Instead, they've just waved
#their hands madly and authored a blog post encouraging Selenium -- which is:
#* Insanely heavyweight.
#* Extremely slow.
#* Occasionally non-deterministic. Okay! It's actually faulty Selenium-driven
#  tests referred to as "flaky" that behave non-deterministically, but the two
#  are effectively one and the same. See this superb 2020 article on the subject
#  for deep insight into what writing non-flaky Selenium tests entails.
#  Basically, you have to embed explicit waits everywhere in Selenium tests to
#  ensure that the DOM has properly updated *BEFORE* actually testing:
#      https://salmonmode.github.io/2020/05/11/is-selenium-flaky.html
#
#Selenium tests are referred to as "end-to-end" (E2E). They're the gold standard
#for web dev testing -- but they're so expensive to author and maintain that
#they typically require a dedicated full-time engineer referred to as a QA
#Automation Engineer. Calculion simply cannot justify that kind of expenditure.
#Every hour is precious here. Let's let it count.
#
#Instead, here's what we're going to do:
#* Define a new "calculion.a90_func.test_lifecycle" submodule. In that:
#  * Define a new test_lifecycle() integration test. This test should either
#    externally fork *OR* possibly just call the streamlit.web.cli.main()
#    function to pretend that we are externally forking a new "streamlit run"
#    subcommand running this web app. The latter approach is quite alluring --
#    and *NOT* simply because it's significantly easier. The latter approach
#    also ensures that uncaught exceptions raised by running our app are
#    implicitly treated as pytest failures, complete with exhaustive tracebacks.
#* Define a medley of unit tests under "calculion.a00_unit" leveraging the
#  third-party "streamlit_mock" package. The idea here is that "streamlit_mock"
#  mocks the *ENTIRETY* of Streamlit, so Streamlit itself is *NEVER* actually
#  involved. That's both good and bad. That's good, because involving Streamlit
#  necessarily involves lots of other obtrusive technology like servers,
#  hostnames, ports, JavaScript, HTML, CSS, and so on that are innately
#  non-trivial to test. (See also: Selenium.) That's bad, because mocks are
#  mocks. Still, that's why we first define our test_lifecycle() integration
#  test. If both pass, we can be reasonably certain that the app at least
#  conforms to superficial expectations. See also:
#      https://github.com/acschofield/streamlit_mock
#FIXME: Also note that Selenium has largely been obsoleted by Cypress, but that
#Cypress is largely unusable in Python and especially under pytest. Enter
#"pylenium", a best-of-breed admixture that aims to expose Selenium to pytest in
#a Cypress-like way. "pylenium" sounds and looks great, but it's *NOT* as
#actively maintained as we'd like. Still, if we ever get into E2E (which seems
#doubtful), then "pylenium" would probably be our first stop along the way.

# ....................{ FIXTURES                           }....................
# Provide app-specific pytest fixtures required by lower-level tests.
