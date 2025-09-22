from setuptools import setup, find_packages

setup(
    name="DRIDmetric",
    version="1.1.0",
    author="Moritz Schaeffler",
    description="Package to calculate the DRID metric from a molecular dynamics trajectory",
    packages=find_packages(
        exclude=[
            "example", "examples", 
            "example.*", "examples.*",
        ]
    ),
    include_package_data=False,   # do not auto-include non-code files
    package_data={}, 
    install_requires=["numpy","tqdm","MDAnalysis"],
)
