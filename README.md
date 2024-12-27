# cofimaker <img src="cofimaker.svg" width="150" height="150"/>
command-line tool for manipulating atomic coordinate files and more

license: AGPLv3

language: fortran

usage: cofima --task [options]

help: cofima --help

installation on  a linux computer:
- download source code
- go to directory with source code
- (if not done yet, load or install a fortran compiler (for example gfortran) and linear algebra library "Lapack") 
- (if needed, edit the file "Makefile")
- type "make"

# Examples
cofima --cut: cut surface/slab/sphere from crystal, saturate dangling bonds with H             
cofima --vasp_bs, cofima --vasp_bs_pr: (atom-projected) electronic band structure (vasp postprocessing)

<img src="https://github.com/skoerbel/cofimaker/blob/main/pictures/HYDR.ZnO_sphere.svg" width="150" height="150"/> <img src="https://github.com/skoerbel/cofimaker/blob/main/pictures/HYDR.ZnO.svg" width="200" height="200"/> <img src="https://github.com/skoerbel/cofimaker/blob/main/pictures/BS_PROJ_UP_Hf.svg" width="200" height="200"/>

# Acknowledgement
Inspiration: Pierre Hirel
Testing: Eva M. Kalivoda, Anais Colibaba
Name: Pavel Márton
