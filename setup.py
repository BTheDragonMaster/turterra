from distutils.core import setup

setup(
    name="turterra",
    version="1.0",
    authors=["Barbara Terlouw", "Janani Durairaj"],
    packages=[
        "turterra",
        "turterra/utils",
        "turterra/dependencies",
        "turterra/dependencies/pikachu",
        "turterra/dependencies/modeller",
    ],
    install_requires=[
        "numpy",
        "numba",
        "scipy",
        "prody",
        "biopython",
        "typer",
        "pyparsing",
        "dash",
        "dash-bio==0.5.*",
        "dash-daq",
        "dash-bio-utils==0.0.6",
        "cryptography",
        "dash-cytoscape==0.2.0",
        "dash-extensions",
        "matplotlib",
        "pandas",
    ],
    extras_require={"model": ["modeller",]},
    scripts=["bin/turterra", "bin/turterra-build"],
)
