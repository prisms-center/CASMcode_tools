__version__ = "1.0a1"

# Available at setup time due to pyproject.toml
from setuptools import setup

setup(
    name="casm-map",
    version=__version__,
    packages=[
        "casm",
        "casm.map",
        "casm.map.commands",
        "casm.map.misc",
        "casm.tools",
        "casm.tools.calc",
        "casm.tools.calc.commands",
        "casm.tools.map",
        "casm.tools.map.commands",
        "casm.tools.shared",
    ],
)
