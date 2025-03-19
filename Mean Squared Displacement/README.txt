Mean Squared Displacement (MSD) is a key metric used in Molecular Dynamics (MD) simulations to analyse the dynamics of water (or any molecule) over time. It measures how far a molecule (or atom) moves from its initial position as a function of time, providing insights into the diffusion behaviour and mobility of the molecules. It measures the average squared distance that a molecule travels from its starting point over time.

Requirements:
1) MDAnalysis
2) numpy 

Additional:
1) Here we have considered oxygen atoms (OW) of water molecule that are around 5.5 Angstrom distance from protein molecule. Change the cutoff distance as per need.
2) The name of the datafiles considered here are msd.pdb and msd.xtc. Change the names accordingly.
3) Datafile created is named as msd.dat
4) Make sure the input files does not contain hydrogen atoms of both protein and water molecule. Process the pdb and xtc files accordingly. If the hydrogen atoms of protein molecule would have been present, then the OW atoms would have been selected 5.5 Angstrom from the hydrogen atoms of protein molecule. Hence, remove the hydrogen atoms of protein molecule.
5) Here we have analyzed all the frames of the entire trajectory. Change it as per needed in the variable name "frames needed".
