      module help

      contains
    
      subroutine print_help(morehelpon)
      implicit none
      character(len=*) morehelpon
        print '(a)', 'usage: cofima [OPTIONS]'
        print '(a)', ''
        print '(a,x,a,x,a,x,a)', 
     &   'Without further options, cofima reads input from',
     &   'INPUT.COFIMA and writes output to OUTPUT.COFIMA.',
     &   'The cofima documentation tells you how to set up an input',
     &   'file.'
        print '(a)', ''
        print '(a)', 'cofima options:'
        print '(a)', ''
        print '(a)', 
     &   '  -v, --version     print version information and exit'
        print '(a)', 
     &   '  -h, --help        print usage information and exit'
        print '(a)', '  -t, --time        print time and exit'
        print '(a)', 
     &   '  --input file1     read input from file file1'
        print '(a)', 
     &   '  --output file2    write output to file file2'

        print '(a)',
     &   '  --add number1 number 2 ',
     &   '                     adds number 1 to number 2.'
        print '(a)',
     &   '  --addatoms file format species nadd outfile coords(1,1:3) ..  &
     &. coords (nadd,1:3) (optional: abs/frac) ',
     &   '                     adds nadd atoms of species species to fil  &
     &e file with format format and writes to outfile. Default: frac'
        print '(a)',
     &   '  --angs2bohr numberinangs',
     &   '                     converts numberinangs to from Angstrom to
     & Bohr.'
        print '(a)',
     &   '  --addeqv infile informat symmformat outfile (top)',
     &   '                     adds atoms to a set of irreducible atoms   &
     &by applying symmetries read from a file "SYMM.DAT".',
     &   '                     "infile" with format "informat" contains   &
     &atom coordinates, "symmformat" indicates if the symmetries in "SYM  &
     &M.DAT" have the mbpp or gulp format, ',
     &   '                     and outfile will contain the added atoms'
        print '(a)',
     &   '  --angs2hart numberinangs',
     &   '                     converts numberinangs from Angstrom to Ha
     &rtree.'
        print '(a)',
     &   '  --av file format firststep step laststep',
     &   '                     calculate averaged lattice parameters',
     &   '                     and write them to the file AV.DAT.',
     &   '                     Example: --av myfile cfg 40 20 100',
     &   '                     if you have myfile40.cfg, myfile60.cfg,',
     &   '                     ...myfile100.cfg.'
        print '(a)',
     &   '  --bandgap infile informat',
     &   '                     reads the bandgap from infile with format
     &informat (supported: nwchem,outcar,outcar_gw).',
     &   '                     usage: e.g. cofima --bandgap nwchem.out n
     &wchem'
        print '(a)',
     &   '  --basis4fiesta ele zetamin zetamax nzeta lmax (rmax "rmax"',
     &   '                     npoints "npoints")',
     &   '                     creates an auxiliary basis (el2.ion) to b
     &e read by fiesta, where el is an element name, zetamin and ',
     &   '                     zetamax are the minimum and mmaximum gaus
     &sian exponents that are to be taken into account, ',
     &   '                     nzeta is the number of different exponent
     &s, lmax the largest angular momentum. rmax (optional, default20)',
     &   '                     is the largest radial gridpoint, npoints 
     &(optional, default 10000) is the number of radial gridpoints.',
     &   '                     The ratio of coeffients is constant.',   
     &   '                     E.g.: ',
     &   '                     cofima --basis4fiesta He 0.8 4.0 8 2 ',
     &   '                     cofima --basis4fiesta He 0.8 4.0 8 2 rmax
     &30.0 npoints 5000'
        print '(a)',
     &   '  --bohr2angs numberinbohr',
     &   '                     converts numberinbohr from Bohr to Angstr
     &om.'
        print '(a)',
     &   '  --broaden file x-column value-column periodic period brd (lc  &
     &ontinuous)',
     &   '                     reads the x values as column "x-column" a
     &nd the function values from column "value-column" of file "file"',
     &   '                     and applies a Gaussian broadening by "brd
     &". If "periodic"==T, periodicity with length "period" is accounted
     & for. If lontinuous==T, a continuous distr. is assumed'
!        print '(a)',
!     &   '  --broaden_continuous file x-column value-column periodic per  &
!     &iod brd',
!     &   '                     same as --broaden but for a continuous di  &
!     &stribution'
        print '(a)',
     &   '  --broaden_lin file x-column value-column periodic period brd  &
     &1 brd2 (lcontinuous)',
     &   '                     same as --broaden but with a linearly cha  &
     &nging width (from brd1 to brd2)'
!        print '(a)',
!     &   '  --broaden_lin_continuous file x-column value-column periodic  &
!     &period brd1 brd2',
!     &   '                     same as --broaden_lin but for a continuou  &
!     &s distribution' 
        print '(a)',
     &   '  --broaden_quad file x-column value-column periodic period br  &
     &d1 brd2 (lcontinuous)',
     &   '                     same as --broaden_lin but broadening incr  &
     &eases quadratically'
        print '(a)',
     &   '  --com file format',
     &   '                    calculates the Center Of Molecule from coo
     &rdinates in "file" which has format "format".'
        print '(a)',
     &   '  --comass file format',
     &   '                    calculates the center of mass from coordin
     &ates in "file" which has format "format".'
        print '(a)',
     &   '  --cageZnO nunits',
     &   '                     creates a cage-like cluster of nunit
     &s ZnO units.'
        print '(a)',
     &   '  --compare file1 file2 ',
     &   '                     calculates overlap <f1|f2> / (|f1|*|f2|)
     &of two fields on a 1d grid, read from file1 and file2. ',
     &   '                     The grids should be similar, though not n
     &ecessarily equal. Each line i in file_j should be rgrid_j(i), f_j(
     &(i)'
        print '(a)',
     &   '  --overlaps file1 file2 ',
     &   '                     calculates overlap <f1|f2> / (|f1|*|f2|)
     &of two fields on a 1d grid, read from file1 and file2. ',
     &   '                     The grids should be similar, though not n
     &ecessarily equal. Each line i in file_j should be rgrid_j(i), f_j(
     &(i)'
        print '(a)',
     &   '  --convert file1 file2 format1 format2 (v11 v12 v13 v21...)',
     &   '                     convert file1 from format1 to format2',
     &   '                     and write to file2.',
     &   '                     Supported formats so far:',
     &   '                     coorat (MBPP), OUT (MBPP), xsf, gin (GULP
     &), data (LAMMPS), cfg, lxyz (LAMMPS), abinit, nwchem, xyz, (cif)',
     &   '                    poscar, outcar, siesta, siesta_out, ascii   &
     &(only format2), poscar_forces, xsf_forces',
     &   '                    findsym, vesta. If format1=coorat/nwchem/x  &
     &yz, give lattice vectors v11...v33.'
        print '(a)',
     &   '  --coord file format (vecs) (-cutoff cutoff) ',
     &   '                     determines min, max, av coord. numbers. I  &
     &f the format is one without lattice vectors, supply them.',         &
     &   '                     Optionally supply the cutoff bond length   &
     &(default 3 Angs)'
        print '(a)',
     &   '  --corr file species property corrf rcut (binwidth rmax brd) 
     &',
     &   '                     calculate a correlation function using ',
     &   '                     an auxiliary property number "property"',
     &   '                     of the species "species" in a cfg file',
     &   '                     "file" using a predefined correlation  ',
     &   '                     function "corrf" with cutoff "rcut" and',
     &   '                     binwidth "binwidth" up to a maximum ',
     &   '                     distance "rmax" and gaussian broadening',  
     &   '                     by "brd".',   
     &   '                     possible corrf: "theta"'   
        print '(a)',
     &   '  --create_neb_chain file1 file2 format N (vecs1 vecs2)',       &
     &   '                     read file1 and file2, which have format "  &
     &format", if needed read vectors, read number N of images inbetween  &
     &', '                     create chain of linearly interpolated ima  &
     &ges'  
        print '(a)',
     &   '  --cross_product vec1(1:3) vec2(1:3)',
     &   '                     does the expected.' 
        print '(a)',
     &   '  --crysanpero ',
     &   '                     execute the routine "crysanpero".',
     &   '                     An input file with instructions for crysa
     &npero is required.',
     &   '                     By default, the program reads "INPUT.CRYS
     &ANPERO". A different input file can be specified via the "input" o
     &ption.',
     &   '                     To see an example of an input file for cr
     &ysanpero, type "help crysanpero".'  
        print '(a)',
     &   '  --cut file1 file2 format1 format2 above/below 1/2/3 cut ',
     &   '                     (nndist (hdist) -hydr)               ',
     &   '                     remove atoms above/below x/y/z=cuthere',
     &   '                     and write to file2.',
     &   '  --cut file1 file2 format1 format2 sphere 0/1 radius x y z',
     &   '                     remove atoms outside(0)/inside(1) of a',
     &   '                     sphere at origin x y z'   ,
     &   '                     A file "HYDR.file2" with hydrogenated ',
     &   '                     dangling bonds is written if requested.',
     &   '                     nndist is the maximum broken bondlength',
     &   '                     which is hydrogenated, hdist is the',
     &   '                     distance along the former bond at which',
     &   '                     the H is placed. Defaults: nndist=3.0 ',
     &   '                     Angstrom, hdist=1.25 Angs. If you wish ',
     &   '                     to specify hdist, you need to specify ',
     &   '                     nndist first. hdist>0: absolute dist. ',
     &   '                     in Angs, hdist<0: fraction of former ',
     &   '                     bondlength.',
     &   '  --cut file1 file2 format1 format2 plane 0/1 distance n1 n2 n  &
     &3',
     &   '                     remove atoms above(0)/below(1) of a plane  &
     &with distance "distance" from origin and direction (n1,n2,n2)',  
     &   '  --cut file1 file2 format1 format2 dirplane 0/1 distance n1 n  &
     &2 n3',
     &   '                     remove atoms above(0)/below(1) of a plane  &
     &with distance "distance" from origin and direction n1*a + n2*b + n  &
     &3 * c (a,b,c = direct lattice vectors)',   
     &   '  --cut file1 file2 format1 format2 hklplane 0/1 distance h k   &
     &l',
     &   '                     remove atoms above(0)/below(1) of a plane  &
     &with distance "distance" from origin specified by Miller indices' 
        print '(a)',
     &   '  --deriv file1 deriv_mode',
     &   '                     calculates first derivative of second col
     &umn in file1 wrt the first one using finite differences.',
     &   '                     "deriv_mode" can be "cd" (central differe
     &nce), "cd2" (second derivative, central diff.), or fps (5-point-st
     &encil)'      
        print '(a)',
     &   '  --detmat3 matrix',
     &   '                     calculate determinant of 3x3 '
        print '(a)',
     &   '  --dfh dim point direction ',
     &   '                     calculates the distance of "point" from t  &
     &he convex hull. The hull facets are read from "HULL_FACETS", they   &
     &must be line segments in 2d and triangles in 3d.',
     &   '                     dim is the dimension (2 or 3), and "direc  &
     &tion" is the direction along which the distance is required. ',     &
     &   '                     If the energy is the z component, use 0 0  &
     & 1 for "direction".   '
        print '(a)',
     &   '  --dipole file format ',
     &   '                     calculates the dipole moment of a ',
     &   '                     cluster of point charges.',
     &   '                     Caution works only for finite systems.',
     &   '                     Result in eV*Angs.'
!        print '(a)',
!     &   '  --dipoleen0 lattice charge dipole epsil alat',
!     &   '                     same as dipoleen, but epsil is a scalar'
!        print '(a)',
!     &   '  --dipoleen lattice charge dipole epsil alat',
!     &   '                     calculates dipole-dipole energy for some
!     &lattices of point dipoles.',
!     &   '                     lattice is one of (cub,tet,fcc,bcc). In t
!     &he case of tet, give a and c, else a is enough. epsil isa matrix',
!     &   '                     dipole is the distance vector between the
!     & two point charges (d/q).',
!     &   '                     Use Bohrs for lengths!'
        print '(a)',
     &   '  --dipen cellvecs dipolevector epsilscalar ',
     &   '                     calculates dipole-dipole energy of a &lat
     &tice of point dipoles.',
     &   '                     Use Bohrs for lengths!'
        print '(a)',
     &   '  --diffatoms file1 file2 format (vecs) (-tol mytol)',
     &   '                     determines common and different ',
     &   '                     atoms for two structure files.',
     &   '                     tolerance mytol for distinction (default
     &0.01)' 
        print '(a)',
     &   '  --dist file format 1/2/3 (-bin 0.01) (-broad 0.1)',
     &   '                     calculate distrib. along x1/x2/x3',
     &   '                     and apply Gaussian broadening.',
     &   '                     default: bin=0.01, broad=0.1 (frac).',
     &   '                     x1 are lattice directions (not necess.',
     &   '                     cartesian).',
     &   '                     A file LATTSPAC.DAT containing the ',
     &   '                     spacings of layers is written as well.',
     &   '                     Caution: no correct results for',
     &   '                     non-orthogonal vectors.'
        print '(a)',
     &   '  --divide number1 number 2 ',
     &   '                     divides number 1 by number 2.'
        print '(a)',
     &   '  --dot_product vec1(1:3) vec2(1:3)',
     &   '                     does the expected.' 
        print '(a)',
     &   '  --ecoulomb file format ',
     &   '                     calculates the electrostatic energy ',
     &   '                     of a cluster of point charges.',
     &   '                     Caution works only for finite systems.',
     &   '                     Result in eV.'
        print '(a)',
     &   '  --ekin file1 format (-step step)',
     &   '                    calculates the kinetic energy of a density
     & read from cube file "file1". Currently "format" must be cube.',
     &   '                    Only every stepth grid point in file is us
     &ed (faster). CAUTION: Works only if WF is real. ',
     &   '                    Default: step=1.'
        print '(a)',
     &   '  --ev2hart number',
     &   '                     converts number from eV to Hartree.'
        print '(a)',
     &   '  --ev2ryd number',
     &   '                     converts number from eV to Rydberg.'
        print '(a)',
     &   '  --eV2nm number',
     &   '                     converts photon energies in eV to light    &
     &wavelengths in nm'
        print '(a)',
     &   '  --expand',
     &   '                     expands (shrinks) cell vectors by a facto
     &r read from command line. Usage:' ,
     &   '                     cofima --expand FILE format OUTFILE forma
     &t factor '
        print '(a)',
     &   '  --expo_range',
     &   '                     determines the range of gaussian exponent
     &s zeta (range=zetamax-zetamin)', 
     &   '                     from a file which should have the fiesta 
     &basis format, e.g. NWCHEM_Ga.ion , He.ion, Cu2.ion etc.',
     &   '                     usage: cofima --expo_range filename'  
        print '(a)',
     &   '  --fourier_1d file1 ftmode (-nq integer -qmax float)',
     &   '                     calculates FT of 2nd column in file1 usin  &
     &g ftmode (options: simple). Optional: number of qpoints nq, maximu  &
     &m q qmax'      
        print '(a)',
     &   '  --fourier_1d_reverse file1 ftmode',
     &   '                     calculates reverse FT of 2nd column in fi  &
     &le1 using ftmode (options: simple)'      
        print '(a)',
     &   '  --gulpnebana infile ',
     &   '                     read GULP NEB data from GULP NEB outputfi
     &le "infile", such as coordinates and energies of the NEB images',
     &   '                     and write the structures to "NEBGEO.XSF"
     & and the energies to "NEBENERGIES.DAT".',
     &   '                     Deprecated. Use nebana instead.'
        print '(a)',
     &   '  --hart2ev number',
     &   '                     converts number from Hartree to eV.'
        print '(a)',
     &   '  --hart2ryd number',
     &   '                     converts number from Hartree to Rydberg.'
        print '(a)',
     &   '  --hcover file1 file2 format',
     &   '                     reads cluster atom coordinates from "file
     &1" and coordinates of embedding atoms (and cluster atoms)',
     &   '                     from "file2". Replaces first the embeddin
     &g atoms (the ones that are in "file2" but not in "file1") by H, th
     &en strips off the redundant H.' 
        print '(a)',
     &   '  --integrate file1',
     &   '                     calculates integral of second column in f  &
     &ile1 wrt the first one'      
        print '(a)',
     &   '  --invmat3 matrix ',
     &   '                     inverts 3x3 "matrix".'
        print '(a)',
     &   '  --ipol file1 ipolstep ipolmeth ',
     &   '                     interpolates second column in file1 wrt t
     &he first, using ipolmeth=lin,polint,cspline,gauss,',
     &   '                     rdfi (or fi, or fourier), or qspline, or 
     &akima',
     &   '                     and finegrid=grid/step.'
        print '(a)', 
     &   '  --isperovskite file format (-cutoff cutoff)',
     &   '                     check whether file "file" with format "fo
     &rmat "format" has perovskite structure.'                    
        print '(a)', 
     &   '  --isperovskite_loose file format ',
     &   '                     check whether file "file" with format "fo
     &rmat "format" has perovskite structure. Same as isperovskite, but   &
     &with looser criteria (more distortion allowed)'                   
        print '(a)', 
     &   '  --isperovskite_loose_analyze file format ',
     &   '                     same as isperovskite_loose, but with addi  &
     &tional analysis (polarization, octahedral tilts) and additional op  &
     &tion -axes ax11 ax12 ax13 ... ax33 (provide axes for octahedral ',  &
     &   '                     rotation, or determine axes automatically  &
     & from A sublattice (-autoaxes A) or B sublattice (-autoaxes B).'   
        print '(a)',
     &   '  --kcalpermole floatnumber',  
     &   '                     converts kcal/mole to eV/molecule.'  
        print '(a)',
     &   '  --latsys file format (tol1 tol2)',
     &   '                     determine lattice system. Use tolerance 
     &tol1 for cell edges, tol2 for angles (cubic if angles=90+-tol2)'  
        print '(a)',
     &   '  --lmpthermo filename',
     &   '                     read MD data from LAMMPS logfile and',
     &   '                     write them to "MD.DAT".'
        print '(a)',
     &   '  --locpotav filename (direction avlength)',
     &   '                     average vasp potential file "filename"',
     &   '                     (usually LOCPOT) along direction "directi
     &on". avlength is the length for the macroscopic average. Ex.:',
     &   '                     cofima --locpot LOCPOT 2 0.5'  
        print '(a)',
     &   '  --matchatoms file1 format1 file2 format2 outfile format3',    &
     &   '                     reorders atoms in file2 such that their p  &
     &ositions match best those in file1'
        print '(a)',
     &   '  --madeldipen lattice charge dipole epsil alat',
     &   '                     calculates Madelung energy for some latti
     &ces of point charges +- q in a homogeneous background charge.',
     &   '                     lattice is one of (cub,tet,fcc,bcc). In t
     &he case of tet, give a and c, else a is enough. Like madelen1,but'
     &,   '                     a dipole basis (like NaCl).',
     &   '                      Use Bohrs for lengths!'
!        print '(a)',
!     &   '  --madelen lattice charge epsil alat',
!     &   '                     calculates Madelung energy for some latti
!     &ces of point charges in a homogeneous background charge.',
!     &   '                     lattice is one of (cub,tet,fcc,bcc). In t
!     &he case of tet, give a and c, else a is enough. Epsilon=scalar.',
!     &   '                     Use Bohrs for lengths!'
        print '(a)',
     &   '  --madelen1 cellvecs charge epsilontensor ',
     &   '                     calculates the screened Madelung energy o  &
     &f a lattice of point charges with Z=1 in a dielectric medium with   &
     &a compensating homogeneous background charge.',
     &   '                     Formula from Rurali & Cartoixa, Nano Lett  &
     &. 9, 975 (2009); Murphy & Hine, PRB 87, 094111 (2013).',
     &   '                     Use Bohrs for lengths!'
        print '(a)',
     &   '  --madelen file format epsilontensor ',
     &   '                     calculates the screened Madelung energy o  &
     &f a lattice of point charges Z_i in a dielectric medium with a hom  &
     &ogeneous compensating background charge.',                          &
     &   '                     Formula from Richargd Martin, Eq. F.5'
        print '(a)',
     &   '  --mateigs n matrix  ',
     &   '                     calculates eigenvalues of nxn matrix.'
        print '(a)',
     &   '  --mat_times_vec matrix(1:3,1:3) vec(1:3) ',                   &
     &   '                     calculates matriy times vector.'
        print '(a)',
     &   '  --matrix_mult n matrix1 matrix2  ',
     &   '                     calculates matrix3=matrix1*matrix2 (nxn m  &
     &atrices only).'
        print '(a)',
     &   '  --mbppnebana mbppout nimg ',
     &   "                     writes NEB data from an MBPP NEB run
     &performed using Bernd Meyer's NEB scripts." ,  
     &   '                     "mbppout" is the name of the MBPP OUT fil
     &e to read from.',
     &   '                     For example, "cofima --mbppnebana OUT.5 8
     &reads structures and energies from the files "img0/OUT.5" to "img8
     &/OUT.5".',
     &   '                     "nimg" is the number of the final NEB ima
     &ge. The routine writes the image structures t
     &o "NEBGEO.XSF" and the energies to "NEBENERGIES.DAT".',
     &   '                     However, the above is still to be impleme
     &nted. At the moment the routine reads energies from "ETOT" and
     &structures from "GEOMETRY.xsf".',
     &   '                     Deprecated. Use nebana instead.'   
        print '(a)',
     &   '  --matrix infile outfile format m11 m12 ... m33',
     &   '                     multiply all coordinates with a matrix'
        print '(a)',
     &   '  --mdana trjfile format (GULPresfile)',
     &   '                     (-bndry 0.09 0.08 0.07)',
     &   '                     analyze trajectory file "file" with ',
     &   '                     format "format".',
     &   '                     Supported input formats so far:',
     &   '                     trg (GULP), HISTORY (DL_POLY).',   
     &   '                     "GULPresfile" is only needed for GULP ',
     &   '                     trajectory files because they do not',
     &   '                     contain the atomic names.',   
     &   '                     "GULPresfile" stands for the name of',
     &   '                     the GULP restart file.' ,     
     &   '                     Output will be written to',
     &   '                     "MD.DAT".',
     &   '                     "-bdry" changes the default for the',    
     &   '                     "cell region that is treated as the ',  
     &   '                     "cell boundary when calculating the',  
     &   '                     "electric polarization.',  
     &   '                     "Only matters in charged cells.',
     &   '                     If a GULP output file from MD named', 
     &   '                     "GOUT" exists, cofima extracts data ',
     &   '                     and writes them to GULPMDOUT.DAT',
     &   '                     and GULPMDAV.DAT.'   
        print '(a)',
     &   '  --merge file1 format1 file2 format2 outfile outfileformat',
     &   '                     merge structure file1 with file 2 and wri
     &te result to file3.'
        print '(a)',
     &   '  --molecule file format iatom (-cutoff -maxiter) ',
     &   '                     detect all atoms that are in  the same mo  &
     &lecule as atom iatom. Connectivity is detected based on distance (  &
     &cutoff, optional) using an iterative scheme (maxiter, optional).'
        print '(a)',
     &   '  --mult n1 n2 n3 file1 file2 format1 format2 ',
     &   '                     create a n1xn2xn3 supercell.'
        print '(a)',
     &   '  --multiply number1 number 2 ',
     &   '                     multiplies number 1 by number 2.'
        print '(a)',
     &   '  --nebana infile informat ',
     &   '                     read NEB data from NEB outputfile "infile
     &". Formats supported so far are gulp (gin) only.',
     &   '                     These are written to "NEBGEO.XSF" and to 
     &"NEBENERGIES.DAT", respectively.'
        print '(a)',
     &   '  --nebchain infile informat outformat nimg',
     &   '                     read chain of structure files and generat
     &e a NEB input file.',
     &   '                     Example: cofima --nebchain GEO xsf gulp 5
     &',
     &   '                     reads GEO0.xsf to GEO5.xsf and writes the
     &m to the gulp input file "CHAIN.GIN". Supported outformat: (gin).'
        print '(a)',
     &   '  --nicecluster file format radius coordmin coordmean dipolema
     &x species1 n1 species2 n2...',
     &   '                     reads bulk file "file" with format "forma
     &t"',                       
     &   '                     and builds cluster with maximum radius',
     &   '                     "radius" in Angstrom, which contains n1',
     &   '                     atoms of element "species1" and so on.',
     &   '                     If radius<0, it is assumed that the',
     &   '                     provided structure file is already ',
     &   '                     that of a sphere or another finite ',
     &   '                     piece of bulk.',
     &   '                     Cluster geometries are printed if the ',
     &   '                     minimum coordination number is at least',
     &   '                     "coordmin" (integer), the averaged ',
     &   '                     coordination number is at least "coordmea
     &n"',                 
     &   '                     and the dipole moment (e*Angs) is not',
     &   '                     larger than "dipolemax"',
     &   '                     Ex.: cofima --nicecluster bulk.xsf xsf 5.
     &0 2 3.0 6.0 Cu 2 Zn 1 Sn 1 Se 4' 
        print '(a)',
     &   '  --nicecluster2 (-test) -file file -format format -radius0 ra
     &dius0 -radius radius -origin x y z-coordmin coordmin -coordmean co
     &ordmean -dipmax dipolemax -atoms species1 n1 species2 n2...',
     &   '                     The same as nicecluster except for an inn
     &er sphere of radius radius0 that is fixed. "test" is optional and
     &causes the routine to stop after the number of configurations has
     & been calculated.'
        print '(a)',
     &   '  --nm2eV number',
     &   '                     converts light wavelengths in nm to photo  &
     &n energies in eV'
        print '(a)',
     &   '  --overlap file_wf1 file_wf2 file_smat',
     &   '                     calculates the overlap between two sets o
     &f wavefunctions wf1 and wf2, whose coefficient wrt a Gaussian',
     &   '                     basis are in the files file_wf1/2, and wh
     &ere the overlap matrix is in file_smat.',
     &   '                     Overlaps are written to overlaps.dat'
        print '(a)',
     &   '  --overlapwf file_wf1 file_wf2 file_smat',
     &   '                     calculates the overlap between two sets o
     &f wavefunctions wf1 and wf2, whose coefficient wrt a Gaussian',
     &   '                     basis are in the files file_wf1/2, and wh
     &ere the overlap matrix is in file_smat.',
     &   '                     Overlaps are written to overlaps.dat'
        print '(a)',
     &   '  --phonopy_print_bands infile outfile ',
     &   '                     reads phonon band structure from specifie
     &d file (e.g. band.yaml) and prints phonon eigenvalues',
     &   ' and if present the eigenvectors to specified outputfile'
        print '(a)',
     &   '  --phonopy_print_eigenvecs infile outfile (-qpoint 7 -band 5   &
     &)',
     &   '                     reads phonon band structure from specifie
     &d file (e.g. band.yaml) and prints phonon eigenvectors',
     &   ' of qpoint 7, band 5 to specified outputfile'
        print '(a)',
     &   '  --pol file format ',
     &   '                     calculate polarization using from positio
     &ns and charges as given in the geometry file'
        print '(a)',
     &   '  --pol_vol ndim',
     &   '                     calculates polyhedral volume. Polyhedron   &
     &facets are read from file "HUlL_FACETS" (see quickhull).'
        print '(a)',
     &   '  --power number1 number 2 ',
     &   '                     result: number1**number2.'
        print '(a)',
     &   '  --quickhull filename dimension',                              &
     &   '                     determine the convex hull of a set of poi  &
     &nts in "dimension" dimensions using the algorithm of Barber, Barbe  &
     &r, Dobkin, Huhdanpaa, ',                                            &
     &'ACM Transactions on Mathematical software 22 no 4, 469-483, 1996'
        print '(a)',
     &   '  --randisp infile informat outfile outformat maxdisp(1:3) ',   &
     &   '                     apply random displacements to atom positi  &
     &ons. Maxdisp >= 0 : fractional, <0: absolute maximal displacements  &
     &.'   
        print '(a)',
     &   '  --rdf filename format (vecs) rmax binwidth broad',
     &   '                     calculate radial distribution functions',
     &   '                     of atoms in file "filename" with format',
     &   '                     "format". Output is written to "RDF.DAT".  &
     & If the format contains no vectors, supply them after the format.'
        print '(a)',
     &   '  --read_symops_vasp',                                          &
     &   '                     read symmetry operations from vasp OUTCAR  &
     & and transform into matrix-translation format'
        print '(a)',
     &   '  --readxcme filename format ',
     &   '                     read xc matrix elements from "infile" ',
     &   '                     which has format "informat"',
     &   '                     (currently only nwchem is supported) '
        print '(a)',
     &   '  --relax_atoms inputfile format outputfile format step',
     &   '                     read structure and forces from inputfile,
     & relax atoms by forces multiplied by step. Formats: poscar_forces,  &
     & xsf_forces'
        print '(a)',
     &   '  --relax_rigid_molecule inputfile format outputfile format st
     &ep',
     &   '                     read structure and forces from inputfile,
     & translate atoms rigidly by sum of forces multiplied by step and',
     &   '                     rotate rigidly by sum of torques * step' 
        print '(a)',
     &   '  --relax_rigid_molecule_rot inputfile format outputfile forma
     &t step',
     &   '                     read structure and forces from inputfile,
     & rotate molecule rigidly by sum of torques * step' 
        print '(a)',
     &   '  --relax_rigid_molecule_rot2 inputfile format outputfile form
     &at step',
     &   '                     read structure and forces from inputfile,
     & rotate molecule rigidly by sum of torques * step. Same as ',
     &   '                     relax_rigid_molecule_rot, but different a
     &verageing'   
        print '(a)',
     &   '  --relax_rigid_molecule_trala inputfile format outputfile for
     &mat step',
     &   '                     read structure and forces from inputfile,
     & translate atoms rigidly by sum of forces multiplied by step'
        print '(a)',
     &   '  --relax_rigid_molecule_breath inputfile format outputfile fo
     &mat step',
     &   '                     read structure and forces from inputfile,
     & translat end or shrink molecule rigidly by average radial force'
        print '(a)',
     &   '  --replace filename format species nreplace newspecies outfil
     &e',
     &   '                     read coordinates from file "filename,"',
     &   '                     replace nreplace "species" atoms by "',
     &   '                     "newspecies" atoms and write "',
     &   '                     resulting structure to "outfile".',
     &   '                     "The atoms that were replaced are',
     &   '                     written to "REPLACED.DAT".',
     &   '                     If nreplace<0, atom no. -nreplace will ',
     &   '                     be replaced.'  
        print '(a)',
     &   '  --rmeqv infile format symformat outfile', 
     &   '                     read coordinates from file "infile"',
     &   '                     which has format "informat", read ',
     &   '                     number of symmetry operations and', 
     &   '                     symmetry operations from "SYMM.DAT",',
     &   '                     and write irred. atoms to "outfile". ',
     &   '                     Reduc. atoms are written to "REDUCIBLE"
     &,',
     &   '                     "symmformat" can be "mbpp" or "gulp".'
        print '(a)',
     &   '  --rotate infile outfile format iaxis angle',
     &   '                     rotate around cartesian axis "iaxis" by',
     &   '                     angle (in degree). iaxis=1 or 2 or 3.',
     &   '                     If iaxis=0, read arbitrary axis.',
     &   '                     example: cofima --rotate POSCAR.in POSCAR  &
     &.out poscar 1 60', 
     &   '                     or cofima --rotate POSCAR.in POSCAR.out p  &
     &oscar 0 1.0 1.0 0.0 60.0'    
        print '(a)',
     &   '  --rotate_matrix dimension matrix rotationmatrix (row-wise)'
        print '(a)',
     &   '  --ryd2ev number',
     &   '                     converts number from Rydberg to eV.'
        print '(a)',
     &   '  --ryd2hart number',
     &   '                     converts number from Rydberg to Hartree.'
        print '(a)',
     &   '  --silent',
     &   '                     the program still talks, but less.'
        print '(a)',
     &   '  --slidwinav infile windowlength ncolumns icolumn',
     &   '                     averages data in column icolumn according
     &over a window of windowlength*(length of x interval infile). File 
     &"infile" has ncolumns data columns.'
        print '(a)',
     &   '  --sortdata infile ncolumns ncolumn',
     &   '                     sorts data in column ncolumn according to
     & size. File "infile" has ncolumns data columns.'
        print '(a)',
     &   '  --spectrum infile informat (-brd broadening -bin binwidth)',
     &   '                     extracts excitation energies from NWCHEM,
     & abinit or fiesta_bse file and writes them to SPECTRUM.DAT.',
     &   '                     broadening and binwidth are optional sinc
     &e they have default values (0.05 eV and 0.01).'
        print '(a)',
     &   '  --traj file format firststep step laststep species number nn
     &dist (-nnana nnana nnatom)',   
     &   '                     calculate trajectory of the atom number',
     &   '                     "number" of species "species" from a',
     &   '                     series of coordinate files.',
     &   '                     Example: --traj myfile cfg 40 20 100 Li 1
     &0 2.9',
     &   '                     if you have myfile40.cfg, myfile60.cfg,',
     &   '                     ...myfile100.cfg. and want the ',
     &   '                     trajectory of the 10th Li atom.',
     &   '                     nndist is the maximum distance of NN',
     &   '                     for calculating the coordination number.'
     &,
     &   '                     nnana (=octet) is optional. If "-nnana" i
     &s given, nnatom is the NN species that is to be analyzed.'
        print '(a,/a)', '  --transform trans transoptions infile outfile
     & informat',
     &   '                     apply transformation "trans" with options
     & "transoptions" to 
     &coordinates in "filein" with format "informat". Fur further help
     & type "cofima --help transform".' 
        print '(a)',
     &   '  --vacs filename format species nvac outfile',
     &   '                     read coordinates from file "filename,"',
     &   '                     create nvac "species" vacancies and"',
     &   '                     write remaining coordinates to"',
     &   '                     "outfile. The atoms that were removed "',
     &   '                     are written to "VACS.DAT".' ,              &
     &   '                     If nvac<0, remove the nvacth atom of that  &
     & species.'
        print '(a)',
     &   '  --vasp_bs (-lattice lattice -nskip nskip)',
     &   '                     get band structure from EIGENVAL. Optiona  &
     &l: specifiy lattice (to get special kpoints right) and number of k  &
     &points to skip.'
        print '(a)',
     &   '  --vasp_bs_pr ',
     &   '                     get projected band structure from PROCAR'
        print '(a)',
     &   '  --vasp_BSE_trans_dens (lambda)',
     &   '                     reads binary vasp output file WAVEDER and  &
     & formatted file BSEFATBANDS and calculates transition densities', 
     &   '                     from Eq. (3) in Körbel, "Optical signatur  &
     &es of defects in BiFeO3, PRM 7 (10), 104402 (2023),',             
     &   '                     url:  https://link.aps.org/doi/10.1103/Ph  &
     &ysRevMaterials.7.104402 , Default: lambda=1'
             ! 
        print '(a)',
     &   '  --vasp_BSE_dens (lambda)',
     &   '                     similar to BSE_trans_dens, but w/o the di  &
     &pole matrix elements' ,
     &                         'Default: lambda=1'
             ! 
        print '(a)',
     &   '  --vasp_CHG_cut_sphere -file file -origin origin -radius radi  &
     &us',   
     &   '                     read charge density file "file" output by  &
     & VASP, write density inside and outside to files. "origin" must be  &
     & three fractional numbers, "radius" the radius in Angs'
        print '(a)',
     &   '  --vasp_CHG_cut_plane -file file -origin origin -nvec nvec', 
     &   '                     read charge density file "file" output by  &
     & VASP, write density below and above to files. "origin" must be th  &
     &ree fractional numbers, "nvec" the normal vector of the plane'
        print '(a)',
     &   '  --vasp_CHG_overlap file1 file2',
     &   '                     read two charge densities from file1 and   &
     &file2 (CHGCAR format) and calculate their 3D overlap integral. Wri  &
     &tes also file CD_abs_diff.vasp with |charge difference|.'   
        print '(a)',
     &   '  --vasp_CHG_lin_comb file1 file2 fac1 fac2',
     &   '                     read two charge densities from file1 and   &
     &file2 (CHGCAR format) and calculate their linear combination fac1   &
     &* CHG1 + fac2 * CHG2, which is written to CHG_LIN_COMB.vasp.'
        print '(a)',
     &   '  --vasp_dos (-dos_tol)',
     &   '                     get DOS from DOSCAR and OUTCAR with optio  &
     &nal parameters for layer properties. dos_tol',
     &   '                     for band edge detection, default 1E-4' 
        print '(a)',
     &   '  --vasp_dos_pr (-layers nlayers -direction idir -origin origi  &
     &n -dos_tol)',
     &   '                     get projected DOS from DOSCAR, PROCAR, an  &
     &d OUTCAR with optional parameters for layer properties. dos_tol',
     &   '                     for band edge detection, default 1E-4' 
        print '(a)',
     &   '  --vasp_get_BSE_EV_xml',
     &   '                     read BSE eigenvalues and oscillator stren  &
     &gths from vasprun.xml'
        print '(a)',
     &   '  --vasp_get_eps_vs_omega (-rotate rotation matrix)',
     &   '                     get frequency-dependent dielectric matrix  &
     & from OUTCAR and write to EPS1/2.DAT. Optional: rotate epsilon',    &
     &   '                     to epsrot = S⁺ eps S'   
        print '(a)',
     &   '  --vasp_get_eps_vs_omega_from_chi (-rotate rotation matrix)',
     &   '                     get frequency-dependent dielectric matrix  &
     & from OUTCAR and write to EPS1/2.DAT. Optional: rotate epsilon',    &
     &   '                     to epsrot = S⁺ eps S. Same as --vasp_get_  &
     &eps_vs_omega but now from ALGO=CHI, w and w/o local field effects'   
        print '(a)',
     &   '  --vasp_get_eps_vs_omega_from_chi_flip_ndiag (-rotate rotatio  &
     &n matrix)',
     &   '                     same as vasp_get_eps_vs_omega_from_chi, b  &
     &ut with the nondiagonal elements of epsilon-with-LFE multiplied by  &
     & -1 '
        print '(a)',
     &   '  --vasp_get_eps_vs_omega_xml (rotation matrix)',
     &   '                     get frequency-dependent dielectric matrix  &
     & from vasp.xml and write to EPS1/2.DAT. Optional: rotate epsilon',  &
     &   '                     to epsrot = S⁺ eps S'   
        print '(a)',
     &   '  --vasp_get_pol ',
     &   '                     get ferroelectric polarization from a Ber  &
     &ry phase calculation with vasp (LCALCPOL). Needs the OUTCAR file'
        print '(a)',
     &   '  --vasp_plot_CHG (-avdir 1/2/3 -cutdir 1/2/3 -cutpos pos -fil  &
     &e file)',
     &   '                     get CHG from VASP charge file "file" and   &
     &write to gnuplottable file CHG4PLOTTING.',
     &   '                     defaults: average along a1 cell vector, c  &
     &ut along a1 at fractional cut position 0.0, file=CHGCAR'   
        print '(a)',
     &   '  --vasp_dE_dk n (-cutdir 1/2/3 -cutpos 0.123)',
     &   '                     prints a file VASP_E_and_dE_dK_n with the  &
     & energy eigenvalues and energy gradients in k-space for the nth ba  &
     &nd.',                                                               &
     &   '                     Needs as input parameter the band index n  &
     &, and the files EIGENVAL, OUTCAR need to be present.',              &
     &   '                     Another file is written with the energy a  &
     &nd energy gradient on a k plane. Defaults: plane at k1=0 (cutdir=1  &
     &, cutpos=0).',                                                      &
     &   '                     This file can be plotted with the splot c  &
     &ommand in gnuplot' 
        print '(a)',
     &   '  --vasp_eps2_from_WAVEDER (-nomega nomega -omegamin omegamin   &
     &-omegamax omegamax -broad broadening -vbands nvbands vbands -cband  &
     &s ncbands cbands -spins nspins spins -vbandrange firstband lastban  &
     &d -cbandrange firstband lastband -kpoints ... -smearing ...)',
     &   '                     reads WAVEDER, EIGENVAL and OUTCAR files   &
     &                         and calculates epsilon2. Optionally speci  &
     &fy number of frequencies, ',                                        &
     &   '                     min and max frequency, number of valence   &
     &bands, conduction bands, and spins and a vector of the requested b  &
     &ands and spins, range of valence or conduction bands. ' ,           &
     &   '                     Example: -nomega 1001 -omegamax 0.0 -omeg  &
     &amin 10.0 -broad 0.025 -vbands 4  1 2 3 4 -cbands 4  5 6 7 8 -spin  &
     &s 1  2',                                                            &
     &   ' ',                                                             &
     &   '                     WARNING: only use without symmetries (ISY  &
     &M=0 or spacegroup #1 / P1 / C_1).',                                 &
     &   '                     smearing=fermi or gaussian'   
        print '(a)',
     &   '  --vasp_eps2_IPA_local_av (-nomega nomega -omegamin omegamin   &
     &-omegamax omegamax -broad broadening -vbands nvbands vbands -cband  &
     &s ncbands cbands -spins nspins spins -vbandrange firstband lastban  &
     &d -cbandrange firstband lastband -kpoints ... -smearing ... -dir d  &
     &ir -nperp nperp)',
     &   '                     reads WAVEDER, EIGENVAL, OUTCAR, and PARC  &
     &HARG.NNNN.KKKK.Ridir files from DensityTool                         &
     &                         and calculates spatially projected epsilo  &
     &n2. Optionally specify number of frequencies, ',                    &
     &   '                     min and max frequency, number of valence   &
     &bands, conduction bands, and spins and a vector of the requested b  &
     &ands and spins, range of valence or conduction bands. ' ,           &
     &   '                     dir: direction ov er which is not average  &
     &d, nperp: number of points in that direction',                      & 
     &   '                     Example: -nomega 1001 -omegamax 0.0 -omeg  &
     &amin 10.0 -broad 0.025 -vbands 4  1 2 3 4 -cbands 4  5 6 7 8 -spin  &
     &s 1  2 -dir 3 -nperp 500',                                          &
     &   ' ',                                                             &
     &   '                     WARNING: only use without symmetries (ISY  &
     &M=0 or spacegroup #1 / P1 / C_1).',                                 &
     &   '                     smearing=fermi or gaussian'   
        print '(a)',
     &   '  --vasp_eps2_from_WAVEDER_LEH (-nomega nomega -omegamin omega  &
     &min -omegamax omegamax -broad broadening -vbands nvbands vbands -c  &
     &bands ncbands cbands -spins nspins spins -vbandrange firstband las  &
     &tband -cbandrange firstband lastband -kpoints ... -smearing ...)',
     &   '                     -efermi_e efermi_e -efermi_h efermi_h',  
     &   '                     reads WAVEDER, EIGENVAL and OUTCAR files   &
     &                         and calculates epsilon2. Optionally speci  &
     &fy number of frequencies, ',                                        &
     &   '                     min and max frequency, number of valence   &
     &bands, conduction bands, and spins and a vector of the requested b  &
     &ands and spins, range of valence or conduction bands. ' ,           &
     &   '                     Example: -nomega 1001 -omegamax 0.0 -omeg  &
     &amin 10.0 -broad 0.025 -vbands 4  1 2 3 4 -cbands 4  5 6 7 8 -spin  &
     &s 1  2',                                                            &
     &   ' ',                                                             &
     &   '                     WARNING: only use without symmetries (ISY  &
     &M=0 or spacegroup #1 / P1 / C_1).',                                 &
     &   '                     smearing=fermi or gaussian'   
        print '(a)',
     &   '  --vasp_kaux KPOINTfile disp1 disp2 disp3',  
     &   '                     reads kpoints from KPOINTSfile and prints  &
     & a file KPOINTS.AUX with the original kpoints plus 6 auxiliary ',   &
     &   '                     kpoints shifted by +- disp1 in the cartes  &
     &ian x direction, +- disp2 in y direction, etc. If OUTCAR is not',   &
     &   '                     present, reciprocal lattice vectors canno  &
     &t be read, and fractional shifts are used.'
        print '(a)',
     &   '  --vasp_kgrid n1 n2 n3 (-shift s1 s2 s3 )',  
     &   '                     prints a file KPOINTS.GRID with a Gamma-c  &
     &entered (default) or shifted (if specified) Monkhorst-Pack k grid   &
     &with n1 x n2 x n3 kpoints' 
        print '(a)',
     &   '  --vasp_kpath npoints startk(1:3) endk(1:3)',  
     &   '                     prints a file KPOINTS.PATH containing a p  &
     &ath from startk to endk. The path has in total npoints+1 points ',  &
     &   '                     (including startk and endk).'
        print '(a)',
     &   '  --vasp_read_WAVEDER ',
     &   '                     print content of WAVEDER to formatted fil  &
     &e WAVEDER.dat'   
        print '(a)',
     &   '  --vasp_write_WAVEDER_ext ',
     &   '                     reads binary WAVEDER file and prints dipo  &
     &le ME and eigenvalues to formatted file WAVEDER_EXTENDED.dat'
        print '(a)',
     &   '  --WAVECAR_ana ',
     &   '                     reads binary WAVECAR file header'
        print '(a)',
     &   '  --vasp_write_BORN (-natoms 2 -atoms 1 5)',
     &   '                     reads dielectric tensor and Born effectiv  &
     &e charges from OUTCAR file and writes them to file BORN',           &
     &   '                     (used by phonopy). Optional: print only B  &
     &EC for atoms 1 and 5'
        print '(a)',
     &   '  --xcoulomb file1 file2 format (-step step -finestep finestep
     & -rmin rmin)',
     &   '                    calculates the classical Coulomb energy of
     & an exciton from wave functions read from cube file1(e) and file2(
     &h). Currently format must be cube.',
     &   '                    Only every stepth grid point in file is us
     &ed (faster). finestep is number of fine mesh points per coarse poi
     &nt near r=0 where Coulomb diverges.',
     &   '                    rmin is radius below which fine grid is us
     &ed. Default: step=4, finestep=2, rmin~=grid spacing.'
        
!        print '(x/a)','If you need further help on an option, type cofim
!     &a --help option, e.g. cofima --help transform.'
        !read(*,*) morehelpon
        select case(morehelpon)
        case ("transform") 
              print '(" ")'
              print '(a)','Examples for using the option "transform":'
              print '(" ")'
              print '(a)','--transform matrix 0.0 1.0 0.0  1.0 0.0 0.0  
     &0.0 0.0 1.0 COORAT.IN COORAT.OUT coorat -mult 1 1 2 (mult has bug,
     & better run cofima --mult, see help)'
              print '(a)','reads coordinates from file "COORAT.IN", whic
     &h has the MBPP (COORAT) format, exchanges x and y coordinates and
     & writes the output to "COORAT.OUT". The "-mult" option is optional
     & and returns an extended supercell, e.g. a 1x1x2 cell here.',
     &              '(mult has bug, better run cofima --mult, see help)'
              print '(" ")'
              print '(a)','--transform vectors 0.0 1.0 1.0  1.0 0.0 1.0 
     &1.0 1.0 0.0   1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 1.0 COORAT.IN COOR
     &AT.OUT coorat'
              print '(a)','reads coordinates from file "COORAT.IN", whic
     &h has the MBPP (COORAT) format, transforms from fcc to sc lattice
     &vectors and writes the output to "COORAT.OUT".'
              print '(" ")'
              print '(a)','--transform shift 0.0 0.0 0.25 COORAT.IN COOR
     &AT.OUT coorat'
              print '(a)','--transform shiftabs 5.0 0.0 1.0 in.xyz out.x
     &yz xyz'
              print '(a)','reads coordinates from file "COORAT.IN", whic
     &h has the MBPP (COORAT) format, shifts the z coordinate by 0.25 an
     &d writes the output to "COORAT.OUT".'
              print '(" ")'
              print '(a)','--transform sc2hex COORAT.IN COORAT.OUT coora
     &t'
              print '(a)','reads coordinates from file "COORAT.IN", whic
     &h has the MBPP (COORAT) format, transforms from sc to hex and writ
     &es the output to "COORAT.OUT".'
              print '(" ")'
              print '(a)','--transform mirror COORAT.IN COORAT.OUT coora  &
     &t 1/2/3'
              print '(a)','reads coordinates from file "COORAT.IN", whic
     &h has the MBPP (COORAT) format, mirrors the y coordinates at x/y/z
     &&=0 and writes the output to "COORAT.OUT".'
              print '(" ")'
        print '(a)',
     &   '  --transform mapback infile outfile format (-tol tol1 tol2 to
     &l3)',
     &   '                     mapback to region roughly between 0 and 1
     & (between -tol and 1-tol, actually).'
              print '(a)','--transform mult 1 2 2 COORAT.IN COORAT.OUT 
     &coorat.'
              print '(a)','reads coordinates from file "COORAT.IN", whic
     &h has the MBPP (COORAT) format, multiplies the unit cell by 2 in t
     &he y and z direction and writes the output to "COORAT.OUT".',
     &   '                     Deprecated and buggy, use cofima --mult i
     &nstead.'
              print '(" ")'
              print '(a)','--transform sort COORAT.IN COORAT.OUT c
     &oorat'
              print '(a)','reads coordinates from file "COORAT.IN", whic
     &h has the MBPP (COORAT) format, sorts the coordinates and writes t
     &he output to "COORAT.OUT".'
              print '(a)','--transform sort COORAT.IN COORAT.OUT c
     &oorat -tol 0.1 0.1 0.1'
              print '(a)','The option "-tol t1 t2 t3" specifies
     & the position difference above which coordinates are assumed to
     & differ. If the option is omitted, the default tolerance is used.
     & It is defined in defs.f.'
        case ("crysanpero") 
        print '(a)',' '
        print '(a)','The following is an example of an input file for cr
     &ysanpero. Required lines are marked with "(r)", optional ones with
     &"(o)":',
     &   ' ',
     &'crysanpero                # command to execute crysanpero (r)',
     &'# CrySAnPero Version 0.17 input file    # comments (o)',
     &' ',
     &'corner         # specify the atoms on the corners of the
     &perovskite unit cell (A site), (r)',
     &'1          # number of species that occupy A (corner) sites (r)',
     &'K       1.12 0.0 0.0    0.0 1.12 0.0    0.0 0.0 1.12     
     &# element, charge tensor, weight (1/#unit cells it belongs to),max
     &. dist. from cell center (as fraction of unit cell edge) (r)',
     &' ',
     &'center    # specify the B sites atoms (r)',
     &'2                                       # number of species on B 
     &sites  (r)',
     &'Nb  5.00 0.0 0.0    0.0 5.00 0.0    0.0 0.0 5.00                 
     &# element name, charge tensor, weight, radius (every element in a 
     &new line)  (r)',
     &'Cu  2.0  0.0 0.0    0.0 2.0  0.0    0.0 0.0 2.0                  
     &# Cu_Nb: 3.0, Cu_K: 1.15',
     &'',
     &'face    # (r)',
     &'1               # number of species on cube faces, next lines: na
     &me, charge tensor,weight, radius (new line for each element) (r)',
     &'O 1 -2.00000000 0.0 0.0  0.0 -2.00000000 0.0   0.0 0.0 -2.0000000
     &0      # element name, charge tensor, weight, radius (every elemen
     &t in a new line) (r)',
     &'O 2 -2.00000000 0.0 0.0  0.0 -2.00000000 0.0   0.0 0.0 -2.0000000
     &0      # ',
     &'O 3 -2.00000000 0.0 0.0  0.0 -2.00000000 0.0   0.0 0.0 -2.0000000
     &0       # ',
     &'',
     &'box     # (r)',
     &'7.424 7.481 7.481 90.0 90.0 90.0        # box edge length of sing
     &le cell in Bohr and angles alpha beta gamma in degree (r)',
     &'4 4 4                                   # number of unit cells in
     &each cartesian direction (r)    ',                               
     &'#7.424 7.424 7.424 85.0 80.0 75.0       # box edge length of sing
     &le cell in Bohr and angles alpha beta gamma in degree ',
     &'',
     &'',
     &'#octahedral_rotation  # print a COORAT file with rotated octahedr
     &a (o)',
     &'# 1.155 1.155 1.155      1 1 1        # rotation angles around x,
     &y,z axis (vector), alternating directions of x,y,z rotations (inte
     &ger vector, 0:no, 1:yes), (r) if octahedral rotation is supposed
     &to be performed',
     &'',
     &'# AFE   # print a COORAT file with AFE displacements (o)',
     &'# 001 0 0 1   0.01  -0.01         1 1 1          # name of AFE mo
     &de (string(3)), direction of displacements (vector), e.g., 0 0 1 m
     &eans displacement along c axis (not necess. along z), (r) if AFE
     &displacements are supposed to be calculated',
     &'                                                 # amplitude (in 
     &single cell lattice constants) for center and face atoms, alternat
     &ing direction of disps (0:no, 1:yes), (r) if AFE displacements are
     &supposed to be calculated',
     &'',
     &'print_displacements             # displacements w.r.t. cubic posi
     &tions are printed to "DISP.vscub.DAT" (o)',
     &'#reference                       # if this keyword is present, th
     &e program reads a reference COORAT file (o)',
     &'polarization                    # If this keyword is present, the
     &program calculates dipole moments and polarization (o)',
     &'',
     &'print_xsf    # print an xsf file (o)'
        !case("q")
        !      return
        case("nomorehelp")
              return
        case default
              print '("Sorry, no further help for this option.")'
              return
        end select
        return
      
       end subroutine print_help

       end module
