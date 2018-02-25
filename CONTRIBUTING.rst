.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/IngoScholtes/pathpy/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.
However to avoid problems down the line, comment on the issue and let us know how
you intend to solve the issue/implement the feature.

Write Documentation
~~~~~~~~~~~~~~~~~~~

pathpy could always use more documentation, whether as part of the
official pathpy docs, in docstrings, or even on the web in blog posts,
articles, and such.

Pathpy uses the numpy docstring format, you can see the guide
`here <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/IngoScholtes/pathpy/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `pathpy` for local development.

1. Fork the `pathpy` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/pathpy.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv pathpy
    $ cd pathpy/
    $ pip install -e .
    $ pip install -r requirements_dev.txt

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests::

    $ flake8 pathpy tests
    $ python setup.py test or py.test

   To get flake8 and py.test, just pip install them into your virtualenv

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.5 and 3.6. Check
   https://travis-ci.org/IngoScholtes/pathpy/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

$ py.test tests.test_pathpy


Deploying
---------

A reminder for the maintainers on how to deploy.

1. Make sure that all tests pass on https://travis-ci.org/IngoScholtes/pathpy .
2. Test that the README.rst is valid rst, so as to avoid rendering issues on
   https://pypi.python.org , the render engine is very strict::
    $ rstcheck README.rst

3. Add an entry in HISTORY.rst, outlining the Changes from the previous version
4. Make sure all your changes are committed (including an entry in HISTORY.rst).
5. Then run the following to sync the tags and correctly increase the version number::

    $ bumpversion patch # possible: major / minor / patch
    $ git push
    $ git push --tags

6. Push new version of pathpy to https://pypi.python.org
   Install `twine` if you haven't done so yet. Remove the `--repository-url` flag
   to push to pypi proper and not the test staging site ::

    $ make dist
    $ twine upload --repository-url https://test.pypi.org/legacy/ dist/*

7. Check https://test.pypi.org/project/pathpy/ to check if everything is fine.
