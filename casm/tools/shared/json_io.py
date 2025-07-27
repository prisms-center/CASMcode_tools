"""JSON input/output"""

import gzip
import json
import os
import pathlib
import tarfile
import typing

import libcasm.xtal as xtal


def pretty_json(
    data: dict,
) -> str:
    """Pretty-print Python dict as a JSON string."""
    return xtal.pretty_json(data)


def printpathstr(path):
    abspath = pathlib.Path(path).resolve()
    try:
        return str(abspath.relative_to(pathlib.Path.cwd()))
    except ValueError:
        return str(abspath)


def read_required(path: pathlib.Path, gz: bool = False):
    path = pathlib.Path(path)
    if path.exists():
        if gz is True:
            with gzip.open(path, "r") as f:
                return json.loads(f.read().decode("utf-8"))
        else:
            with open(path, "r") as f:
                return json.load(f)
    else:
        raise FileNotFoundError(
            "Required file: '" + printpathstr(path) + "' does not exist"
        )


def read_contents(
    parent_dir: pathlib.Path,
    relpath: pathlib.Path,
    default: typing.Any = None,
    quiet: bool = False,
):
    """Read file contents, whether it is an individual or inside an archive

    Parameters
    ----------
    parent_dir: pathlib.Path
        The file path to a Monte Carlo run path directory

    relpath: pathlib.Path
        Relative path, referenced to parent_dir, of a file to open. The `parent_dir`
        directory may be tar gzipped, in which case this file is read from the archive.

        If `relpath.suffix == ".gz"`, the file is read as a gzipped file.

    Returns
    -------
    data: typing.Any
        The JSON contents of the file, or a default value.
    """
    parent_dir = pathlib.Path(parent_dir)
    relpath = pathlib.Path(relpath)
    datafile_path = parent_dir / relpath
    tgz_path = parent_dir.parent / (parent_dir.name + ".tgz")

    try:
        if datafile_path.exists():
            if relpath.suffix == ".gz":
                # gzipped JSON file
                with gzip.open(datafile_path, "rb") as f:
                    return json.loads(f.read().decode("utf-8"))
            else:
                # JSON file
                with open(datafile_path, "r") as f:
                    return json.load(f)
        elif tgz_path.exists():
            with tarfile.open(tgz_path, "r:gz") as f:
                membername = str(pathlib.Path(parent_dir.name) / relpath)
                with f.extractfile(membername) as g:
                    if relpath.suffix == ".gz":
                        # gzipped JSON file
                        return json.loads(gzip.open(g).read().decode("utf-8"))
                    else:
                        # JSON file
                        return json.load(g)
        else:
            if not quiet:
                print(printpathstr(datafile_path) + ": does not exist, skipping")
            return default
    except Exception:
        if not quiet:
            print(printpathstr(datafile_path) + ": error reading, skipping")
        return default


def read_optional(path: pathlib.Path, default: typing.Any = None, gz: bool = False):
    path = pathlib.Path(path)
    if path.exists():
        if gz is True:
            with gzip.open(path, "r") as f:
                return json.loads(f.read().decode("utf-8"))
        else:
            with open(path, "r") as f:
                return json.load(f)
    return default


def read_cascading(paths: list[pathlib.Path], quiet: bool = False, gz: bool = False):
    """Find a required file that may be at multiple locations"""
    for path in paths:
        data = read_optional(path, gz=gz)
        if data is not None:
            if not quiet:
                print("reading:", os.path.relpath(path))
            return data
    print("Required file not found. Checked:")
    for path in paths:
        print("- " + printpathstr(path) + ": does not exist")
    raise Exception("Required file does not exist")


def dump(
    data, path: pathlib.Path, force: bool = False, quiet: bool = False, gz: bool = False
):
    """Json dump with overwrite/skipping/write output messaging"""

    def _write(data, path: pathlib.Path, gz: bool = False):
        path.parent.mkdir(parents=True, exist_ok=True)
        if gz is True:
            with gzip.open(path, "w") as f:
                f.write(json.dumps(data).encode("utf-8"))
        else:
            with open(path, "w") as f:
                json.dump(data, f)

    path = pathlib.Path(path)
    if path.exists():
        if force:
            if not quiet:
                print("overwrite:", printpathstr(path))
            _write(data, path, gz=gz)
        elif not quiet:
            print("skipping:", printpathstr(path))
    else:
        if not quiet:
            print("write:", printpathstr(path))
        _write(data, path, gz=gz)


def safe_dump(
    data,
    path: pathlib.Path,
    force: bool = False,
    quiet: bool = False,
    gz: bool = False,
):
    """Json dump with overwrite/skipping/write output messaging

    Writes to the temporary file `path + ".tmp"`, then removes `path`
    and renames the temporary file, to avoid losing the original file
    without writing the new file. This method does not avoid race
    conditions.

    If the temporary file already exists an exception is raised.
    """
    path = pathlib.Path(path)

    def _safe_write(data, path, gz: bool = False):
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = pathlib.Path(str(path) + ".tmp")
        if tmp_path.exists():
            raise Exception("Error: " + str(tmp_path) + " already exists")

        if gz is True:
            with gzip.open(tmp_path, "w") as f:
                f.write(xtal.pretty_json(data).encode("utf-8"))
        else:
            with open(tmp_path, "w") as f:
                f.write(xtal.pretty_json(data))

        if path.exists():
            path.unlink()
        tmp_path.rename(path)

    if path.exists():
        if force:
            if not quiet:
                print("overwrite:", printpathstr(path))
            _safe_write(data, path, gz=gz)
        elif not quiet:
            print("skipping:", printpathstr(path))
    else:
        if not quiet:
            print("write:", printpathstr(path))
        _safe_write(data, path, gz=gz)


def get(data: typing.Union[dict, list], path: list, default: typing.Any = None):
    """Get item from a heirarchical data structure of dict and lists

    Parameters
    ----------
    data: typing.Union[dict, list]
        A heirarchical data structure of dict and lists
    path: list
        A path of keys and indices into the data
    default: typing.Any = None
        A default value to return if a key does not exist. Invalid indices into
        lists raise Exceptions instead of returning the default value.
    """
    if isinstance(data, dict):
        key = str(path[0])
        if key in data:
            if len(path) == 1:
                return data[key]
            else:
                return get(data[key], path=path[1:], default=default)
        else:
            return default
    elif isinstance(data, list):
        index = int(path[0])
        if index < 0:
            raise Exception("Error in json_io.get: index < 0")
        elif index >= len(data):
            raise Exception("Error in json_io.get: index > len(data)")
        else:
            if len(path) == 1:
                return data[index]
            else:
                return get(data[index], path=path[1:], default=default)
    else:
        raise Exception("Error in json_io.get: not a dict or list")
