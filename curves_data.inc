      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

c...Maximum number of atoms in complex
      parameter (n1=20000)

c...Maximum number of atoms including hydrogens
      parameter (n1p=30000)

c...Maximum number of subunits
      parameter (n2=1500)

c...Maximum number of nucleotides
      parameter (n3=500)

c...Number of backbone axis points/step for groove calculations
      parameter (n4=5)

c...Maximum number of data stored per level for trajectory analysis
      parameter (n5=100)

c...Maximum number of ligands to analyse
      parameter (n6=100)

c...Delta change of variables for numerical gradient test 
      parameter (deltv=1.d-8)

c...Minimal value of dot product to accept an ion/ligand position (1°)
      parameter (acrit=0.0174d0)

c...Conversion factor from degrees to radians 
      parameter (cdr=0.017453293d0)	

c...Conversion factor from radians to degrees 
      parameter (crd=57.29577951d0)

c...Constant: Pi
      parameter (pi=3.141592654d0)

c...Constant: Precision for Eigenvalue calculation
      parameter (range=1.d-8)
