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
cofima --cut <img src="https://github.com/skoerbel/misc/CZTSe_15fu.svg" width="50" height="50"/> <img src="https://github.com/skoerbel/misc/HYDR.ZnO.svg" width="50" height="50"/>

