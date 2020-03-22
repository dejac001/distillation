import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="distillation",
    version="0.0.1",
    author="Robert F. De Jaco",
    author_email="dejac001@umn.edu",
    description="Distillation for Chemical Engineers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dejac001/distillation",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_reqs=["numpy==1.18.1","scipy == 1.4.1"]
)
