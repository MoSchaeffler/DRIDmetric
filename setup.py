from setuptools import setup, find_packages

setup(
    name="DRIDmetric",
    version="1.1.0",
    author="Moritz Schaeffler",
    description="Package to calculate the DRID metric from a molecular dynamics trajectory",
    packages=find_packages(),
    install_requires=["numpy","tqdm","MDAnalysis"],
)
