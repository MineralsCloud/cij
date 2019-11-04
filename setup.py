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
        "numpy",
        "scipy",
        "qha"
    ]
)
