"""Text file input/output"""

import gzip
import pathlib

from .json_io import (
    printpathstr,
)


def safe_write(
    text: str,
    path: pathlib.Path,
    force: bool = False,
    quiet: bool = False,
    gz: bool = False,
):
    """Write text with overwrite/skipping/write output messaging

    Writes to the temporary file `path + ".tmp"`, then removes `path`
    and renames the temporary file, to avoid losing the original file
    without writing the new file. This method does not avoid race
    conditions.

    If the temporary file already exists an exception is raised.

    Parameters
    ----------
    text : str
        The text to write to the file.
    path : pathlib.Path
        The path to the file where the text will be written.
    force : bool = False
        If True, the file will be overwritten if it already exists. If False, the file
        will not be overwritten if it exists. Default is False.
    quiet : bool = False
        If True, no messages will be printed about writing the file. If False, messages
        will be printed. Default is False.
    gz : bool = False
        If True, the text will be written to a gzip-compressed file. If False, it will
        be written as plain text. Default is False.
    """
    path = pathlib.Path(path)

    def _safe_write(text, path, gz: bool = False):
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = pathlib.Path(str(path) + ".tmp")
        if tmp_path.exists():
            raise Exception("Error: " + str(tmp_path) + " already exists")

        if gz is True:
            with gzip.open(tmp_path, "w") as f:
                f.write(text.encode("utf-8"))
        else:
            with open(tmp_path, "w") as f:
                f.write(text)

        if path.exists():
            path.unlink()
        tmp_path.rename(path)

    if path.exists():
        if force:
            if not quiet:
                print("overwrite:", printpathstr(path))
            _safe_write(text, path, gz=gz)
        elif not quiet:
            print("skipping:", printpathstr(path))
    else:
        if not quiet:
            print("write:", printpathstr(path))
        _safe_write(text, path, gz=gz)
