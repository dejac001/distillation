import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="distillation",
    version="0.1.0",
    author="Robert F. De Jaco",
    author_email="dejac001@umn.edu",
    description="Distillation for Chemical Engineers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dejac001/distillation",
    packages=['distillation'],
    package_data={'distillation': [
        'equilibrium_data/depriester.csv', 'equilibrium_data/heat_capacity_liquid.csv',
        'equilibrium_data/heats_of_vaporization.csv',
    ]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["numpy==1.18.1","scipy == 1.4.1"]
)
