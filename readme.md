# Dual-Site Targeting by Peptide Inhibitors of the N-Terminal Domain of Hsp90: Mechanism and Design

This repository contains the implementation of the machine learning, MD and MMPBSA for the paper:
Dual-Site Targeting by Peptide Inhibitors of the N-Terminal Domain of Hsp90: Mechanism and Design.

## Machine Learning steps

Follow the steps below to reproduce the CAMP prediction results:

1. Write peptides and proteins sequences into two files: peps.txt and prots.txt.
2. Run gen_input.py to generate input data.
3. Run infer.py to predict with five-fold models.
4. Run higher.py to analyze the results.


## MD steps

Follow the steps below to reproduce the MD results:

1. Prepare model files according to files in MD_example.
2. Run sub.sh to generate trajectory.

## MMPBSA steps

Follow the steps below to reproduce the MMPBSA results:

Prepare the following files to run MMPBSA:
    
1. The input file.

The input file MMPBSA_example\mmpbsa.in is already included in this repository

2. The MD Structure+mass(db) and the trajectory files

    cpptraj program from Amber was used to process the psf and dcd files:

        cpptraj -p HRDMYDD_model_1.psf
        >trajin gmxHRDMYDD_model_1.dcd
        >strip :POT,CLA,TIP3,SOD
        >trajout gromacs.pdb onlyframes 1
        >trajout traj.xtc
        >run
        >exit

    After executing these commands, there should be two new files in te folder, _i.e._ `gromacs.pdb` that is 
    going to be used as the MD Structure+mass(db) file, and `traj.xtc` that is going to be used as the 
    trajectory file.

3. The topology file

    We are going to use `ParmEd` to convert the *.psf file into a GROMACS topology file. 

    To do so, use the script that is already included in the `MMPBSA_example` folder. 

```bash
python script.py
```

4. The index file

    The last file we need to generate is the index file containing the groups with the receptor and ligand atoms.

    To do so, just use make_ndx from GROMACS and the MD Structure+mass(db) that was generated previously.

        gmx make_ndx -f gromacs.pdb -o index.ndx
        >a 1-3371
        >a 3372-3489
        >q

    After this, there should be two new groups (10, 11) in the index file that contains atoms from 1 to 3371 for 
    the receptor and atoms 3372 to 3489 for the ligand.

Once the files have been generated, the program can be run either in serial or in parallel:

Serial:
```bash
        gmx_MMPBSA -O -i mmpbsa.in -cs gromacs.pdb -ct traj.xtc -ci index.ndx -cg 10 11 -cp gromacs.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv
```

Parallel:
```bash
        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs gromacs.pdb -ct traj.xtc -ci index.ndx -cg 10 11 -cp gromacs.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv
```