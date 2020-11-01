# Blender Atomic Loader

This is a simple library that allows to load atomic data into blender using ASE and draw spheres for atoms and cylinders for bonds. This README is meant to be a rimender for myself on how to setup blender with ASE, and how to use Blender to render a simple PDB or XYZ (or anything supported by ASE).

**All this have been tested with Blender 2.83 and 2.90**

## Requirements

To load atomic structures from a PDB (or XYZ) file ASE is need. Blender has an internal python interpreter which has nothing to do with that of the system, thus we need to it manually.
