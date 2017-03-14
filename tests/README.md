# Unit tests for pypath

This directory contains the unit tests for methods and functions 
in pathpy.
The testing framework [pytest](doc.pytest.org/) 
is required to run the tests.

To run the test suite (without slow tests) run
```bash
$ pytest tests
```

## Slow functions

Slow functions can be decorated with `slow` to mark them 
as skippable if you require only a quick check.
To run all tests add the flag `--runslow`:
```bash
$ pytest --runslow
```

## Coverage report

To compute a coverage report of the tests you need to install 
[coverage.py](https://coverage.readthedocs.io/en/coverage-4.3.4/)
as well as its `pytest` integration 
[pytest-cov][1]
```bash
$ pytest tests/ --runslow --cov=pathpy --cov-report html
```
which will create an html coverage report in the same directory.

[1]: https://pypi.python.org/pypi/pytest-cov
