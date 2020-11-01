# Blender Atomic Loader

This is a simple library that allows to load atomic data into blender using ASE and draw spheres for atoms and cylinders for bonds. This README is meant to be a rimender for myself on how to setup blender with ASE, and how to use Blender to render a simple PDB or XYZ (or anything supported by ASE). 

The functions in here are pretty simple and most of them are meant for 2D systems made by a metal substrate with a molecule on top. However, this is just an example and can easily be extend to work with other systems.

**All this have been tested with Blender 2.83 and Blender 2.90**

## Importing ASE from Blender

To load atomic structures from a PDB (or XYZ) file [ASE](https://wiki.fysik.dtu.dk/ase/) is need. ASE is an extremely useful library to manipulate atomic data. Unfortunately, Blender has an internal python interpreter which has nothing to do with that of the system, thus we have to add ASE manually.

The first step is too check what Pyhton version Blender comes with. To check this, just start Blender and open the Python console:

![Python version for Blender 2.90](.imgs_readme/python_version.png)

The idea is to install ase using pip in a local virtual environment having the same exact version of Python as Blender:

```bash
$ pyenv virtualenv 3.8.5 ase
$ pyenv local ase
$ pip install ase
```

If this works out smoothly, we can proceed copying the necessary modules to your local Blender path:


```bash
$ mkdir -p ~/.config/blender/2.90/scripts/modules
$ cp -r $PATH_VIRTUALENV/lib/python3.8/site-packages/ase ~/.config/blender/2.90/scripts/modules
$ cp -r $PATH_VIRTUALENV/lib/python3.8/site-packages/scipy ~/.config/blender/2.90/scripts/modules
$ cp -r $PATH_VIRTUALENV/lib/python3.8/site-packages/scipy.libs ~/.config/blender/2.90/scripts/modules
```

Everything should work now. You can remove the local virtual environment, open blender and load ase:

![Test ASE import](.imgs_readme/test_ase_import.png)

## Importing this library from Blender

Clone this repository and open Blender. The library can be loaded using `importlib`:

```python
import importlib.util
 
spec = importlib.util.spec_from_file_location("blender_atomic_loader", "$PATH_TO_Blender_atomic_loader/blender_atomic_loader.py")
baloader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(baloader)
```

Simply copy and paste the lines above in your Blender's python console (changing the correct to the git folder) and you will be able to use the functions need to parse atomic structures and draw the corresponding objects from the command line. N.B: always use `bloader.` before the function's name:


```python
# Import ASE and Numpy
from ase.io import read
import numpy as np

# Read an example system
frame=read('example.pdb')

# Extract the molecule and discard the substrate
molecule=bloader.get_molecule(frame)
```

## Example usage

Follow this simple example to render an image of a small molecule.

## Possible issues

When rendering from a laptop it can happen that the memory is not enough to render the image. The first suggestion is to remove all the atoms not visble from the camera view and the second is to render the image from the command line, without using the GUI (even better if you ssh to a larger machine with blender installed!):

```bash
blender -b test.blend -o output_name -f 1
```
