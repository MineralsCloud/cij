from setuptools import setup, find_packages

setup(
    name='cij',
    version='1.0',
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
        "click"
    ],
    package_data={
        "cij/data/schema": "*.schema.json",
        "cij/data/output": "writer_rules.yml",
        "cij/data/constraints": "*"
    },
    entry_points = {
        'console_scripts': [
            'qha-cij=cij.cli.main:main',
            'qha-cij-extract=cij.cli.extract:main',
            'qha-cij-plot=cij.cli.plot:main',
            'qha-cij-modes=cij.cli.modes:main',
            'qha-cij-static=cij.cli.static:main'
        ],
    }
)
