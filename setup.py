from setuptools import setup, find_packages

setup(
    name="primelt",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "matplotlib",
        "numpy"
        # add other dependencies here
    ],
    entry_points={
        "console_scripts": [
            "primelt=primelt.cli:main",
        ],
    },
)
