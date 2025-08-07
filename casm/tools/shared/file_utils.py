import pathlib
import shutil
import tarfile

from .contexts import working_dir
from .json_io import printpathstr


def _compressed_file(
    dir: pathlib.Path,
    extension: str = ".tgz",
) -> pathlib.Path:
    return pathlib.Path(str(dir) + ".tgz")


def compress(
    dir: pathlib.Path,
    quiet: bool = False,
    remove_dir: bool = True,
    extension: str = ".tgz",
):
    """Compress a directory into a tar gzipped archive file

    Parameters
    ----------
    dir: pathlib.Path
        The directory to compress. The directory must exist.
    quiet: bool = False
        If True, suppresses output messages during compression and removal.
    remove_dir: bool = True
        If True, removes the original directory after compression.
        If False, keeps the original directory.

    """
    dir = dir.resolve()

    if not dir.exists():
        raise FileNotFoundError(f"Directory {dir} does not exist.")

    if not quiet:
        print("compressing:", printpathstr(dir))
    tgz_file = _compressed_file(dir=dir, extension=extension)
    with working_dir(dir.parent):
        tar = tarfile.open(tgz_file, "w:gz")
        tar.add(dir.name)
        tar.close()
    if not quiet:
        print("created:", printpathstr(tgz_file))
    if remove_dir:
        if not quiet:
            print("removing:", printpathstr(dir))
        shutil.rmtree(dir)


def uncompress(
    tgz_file: pathlib.Path,
    quiet: bool = False,
    remove_tgz_file: bool = True,
):
    r"""Uncompress an archive (\*.tar.gz or \*.tgz) file into its parent directory

    Parameters
    ----------
    output_dir: pathlib.Path
        The directory where the tar.gz file is located and where the contents will be
        extracted. The parent directory of `output_dir` should contain the tar.gz file.
    quiet: bool = False
        If True, suppresses output messages.
    remove_tgz_file: bool = True
        If True, removes the archive after extraction. If False, keeps the archive.
    """
    tgz_file = tgz_file.resolve()

    if not tgz_file.exists():
        raise FileNotFoundError(f"Archive {tgz_file} does not exist.")

    if not quiet:
        print("uncompressing:", printpathstr(tgz_file))
    with working_dir(tgz_file.parent):
        tar = tarfile.open(tgz_file, "r:gz")
        tar.extractall()
        tar.close()
    if remove_tgz_file:
        if not quiet:
            print("removing:", printpathstr(tgz_file))
        tgz_file.unlink()
