from setuptools import setup, find_packages
from pathlib import Path

with open(Path(__file__) / "cij" / "version.py") as fp: exec(fp.read())

setup(
    name='cij',
    version=__version__,
    description='High temperature thermal elasticity',
    author='Chenxing Luo',
    author_email='chenxing.luo@columbia.edu',
    url='https://github.com/MineralsCloud/qha-cij/',
    packages=find_packages(),
    install_requires=[
        "numpy >= 1.10.0",
        "pandas",
        "scipy",
        "qha",
        "lazy_property",
        "pint >= 0.10",
        "networkx",
        "click",
        "jsonschema",
        "sympy"
    ],
    package_data={
        "cij.data": [
            "schema/*.schema.json", 
            "output/writer_rules.yml", 
            "constraints/*"
        ]
    },
    entry_points = {
        'console_scripts': [
            'cij=cij.cli.cij:main',
            'cij-run=cij.cli.main:main',
            'cij-run-static=cij.cli.static:main',
            'cij-extract=cij.cli.extract:main',
            'cij-plot=cij.cli.plot:main',
            'cij-modes=cij.cli.modes:main',
            'cij-fill=cij.cli.fill:main',
        ],
    }
)
