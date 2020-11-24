import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blender-atomic-loader",
    version="0.1.0",
    author="Nanotech@surfaces, Empa",
    description="Package to load atomic data to Blender",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nanotech-empa/blender-atomic-loader",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "ase",
    ]
)