# spartan-tpms
Code to generate TPMS structures for the SPARTAN collaboration between Oceanit, UCB, and Marquette.

## Installation
This code is written and tested in Python 3.11. Packages required to run this code are listed in `requirements.txt`. To install these packages, run `pip install -r requirements.txt`. Consider using a virtual environment to avoid conflicts with other packages.

## Usage
1. Command line tool

STL files of Diamond and Gyroid "boxes" can be generated with the command line tool `tpmsgenerator.py` like so:

```
> python tpmsgenerator.py test_gyroid.stl gyroid 1 1 1 0 0 0 0.5
```

These argument must be provided in order to generate an stl:
1. filename of the stl to be generated
2. type of TPMS ("gyroid", "sheet_gyroid", or "diamond")
3. lambda_x -- the "wavelength" in the x dimension
4. lambda_y
5. lambda_z
6. theta_x -- rotation about the x axis (degrees)
7. theta_y
8. theta_z
9. porosity -- fraction of free volume in the structure (0-1)

There are optional arguments as well:
10. `-s` or `--step_size` can control the resolution of the STL file. Defaults to 0.2. This is in the same units as lambda_x,y,z.
11. `-n` or `--num_periods` can control the number of unit cells in the STL file. 
Defaults to 4. This is a multiple of lambda_x,y,z. It can be a float value.
12. `-nx`, `-ny`, `-nz` can control the number of unit cells independently in the x, y, and z directions.  
13. `-sx`, `-sy`, `-sz` can control the size of the rectangular prism to be rendered in millimeters. If this is used, `-nx`, `-ny`, `-nz`, and `-n` cannot be used.

```
> python tpmsgenerator.py test_diamond.stl diamond 1 2 0.4 0 10 34 0.5 -s 0.01 -n 10
```

The line above will generate a diamond TPMS with wavelengths of 1, 2, and 0.4 mm in the x, y, and z dimensions. This lattice is rotated by 10 degrees about y and 34 degrees about z. Finally, a rectangular prism of size 10 x 20 x 4 mm is rendered to an STL with a step size of 0.01 mm (this will take a long time).


```
> python tpmsgenerator.py test_diamond.stl diamond 1 1 1 0 0 0 -nx 1 -ny 1 -nz 10
```
This line will generate a diamond structure with a single unit cell in the x and y directions, and 10 unit cells in the z direction. You must provide all three of `-nx`, `-ny`, and `-nz` to use this option. 

You can always run `> python tpmsgenerator.py --help` for more guidance.

2. Python interface

You can look at the Jupyter notebook `test.ipynb` for some usage examples of this code.