"""Helpers for context management"""

import contextlib
import io
import os
import pathlib
import sys
import typing


@contextlib.contextmanager
def working_dir(wd: typing.Union[str, pathlib.Path, None] = None):
    """Context manager for performing a task with a different working directory.

    .. rubric:: Usage

    .. code-block:: python

        # do some things in a different working directory, then return to the
        # the current working directory
        with working_dir(wd=some_path):
            ... some code ...

    Parameters
    ----------
    wd: Union[str, pathlib.Path, None] = None
        The working directory to use while in context.
    """
    orig_wd = os.getcwd()
    if wd is None:
        wd = orig_wd
    os.chdir(wd)
    try:
        yield
    finally:
        os.chdir(orig_wd)


@contextlib.contextmanager
def captured_output(wd: typing.Union[str, pathlib.Path, None] = None):
    """Capture standard output and error in StringIO objects.

    .. rubric:: Usage

    .. code-block:: python

        # do some things and capture standard output and error
        # as StringIO objects sout, serr
        # optionally, pass 'wd', a different working directory in
        # which to use while executing the code
        with captured_output(wd=some_path) as (sout, serr):
            ... some code ...
        # print sout and serr nicely:
        print_stringIO(sout)
        print_stringIO(serr)

    Parameters
    ----------
    wd: Union[str, pathlib.Path, None] = None
        The working directory to use while in context.
    """
    with working_dir(wd):
        new_out, new_err = io.StringIO(), io.StringIO()
        old_out, old_err = sys.stdout, sys.stderr
        try:
            sys.stdout, sys.stderr = new_out, new_err
            yield sys.stdout, sys.stderr
        finally:
            sys.stdout, sys.stderr = old_out, old_err


def print_stringIO(strio: io.StringIO):
    """Print a StringIO object value nicely with "----" on the lines before and
    after.
    """
    print("\n----\n", strio.getvalue(), "\n----")
