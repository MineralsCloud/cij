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
        "numpy >= 1.7.0",
        "scipy",
        "qha",
        "lazy_property",
        "pint >= 0.10",
        "networkx"
    ],
    package_data={
        "cij/data/schema": "*.schema.json",
        "cij/data/output": "writer_rules.yml"
    },
    entry_points = {
        'console_scripts': [
            'qha-cij=cij.cli.main:main',
            'qha-cij-plot=cij.cli.plot:main'
        ],
    }
)
