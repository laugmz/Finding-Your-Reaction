from setuptools import setup, find_packages

setup(
    name="Finding Your Reaction",
    version='0.1',
    packages=find_packages(),
    scripts=["finding_your_reaction", "functions", "test_cleantext", "test_compare", "test_findrow", "test_isomer","test_issmiles", "test_mapping", "test_molecularweight", "test_namefromespider", "test_nametosmiles", "test_plot3D", "test_randmoness", "test_removepercentage"],
    install_requires=["pandas", "matplotlib", "jupyter lab", "chemspipy", "rdkit", "pathlib", "pubchempy"
        # List dependencies here
    ],
    # Other metadata
    author='Lilou Buffet, Eva Guleac, Laura Gomez',
    author_email='lilou.buffet@epfl,ch, eva-maria.guleac@epfl.ch, laura.gomez@epfl.ch',
    description='This package can be used to find the reaction of formation of a molecule, given it is contained in our database.',
    url='https://github.com/laugmz/Finding-Your-Reaction',
)
