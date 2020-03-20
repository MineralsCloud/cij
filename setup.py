from setuptools import setup, find_packages

setup(
    name='cij',
    version='1.0',
    description='High temperature thermal elasticity',
    author='Chenxing Luo',
    author_email='chenxing.luo@columbia.edu',
    url='https://github.com/chazeon/qha-cij/',
    packages=find_packages(),
    install_requires=[
        "numpy >= 1.7.0",
        "scipy",
        "qha",
        "lazy_property",
        "jsonschema",
        "pint @ git+https://github.com/hgrecco/pint.git#egg=pint"
    ],
    entry_points = {
        'console_scripts': ['qha-cij=cij.cli.main:main'],
    }
)
