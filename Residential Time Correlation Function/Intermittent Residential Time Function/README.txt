Here, the water molecule is allowed to temporarily leave and return within a certain time window. If the molecule revisits the same hydration site, the correlation is still considered intact.

There are 2 ways of calculating Intermittent Residential Time Correlation Function:
1) Run the Python Script named "intermittent_residential_time_correlation_function.py"
2) First run the TCL script named "protein_hydration_water_index.tcl". Then run the Python script named "intermittent_residential_processed.py"

METHOD 1:

Takes longer time to complete. It will create output file named "intermittent_residential_time.dat". Change the name of input files and cutoff values if needed. Process input files without hydrogen atoms from both protein and water molecules.
Requirements:
1) MDAnalysis
2) Numpy

METHOD 2:

Takes quicker time to complete. TCL script creates output file named "hydration_water_index.dat" and then the python script creates final output file named "intermittent_residential_time.dat". Change the name of input files and cutoff values if needed. Process input files without hydrogen atoms from both protein and water molecules.
Requirements:
1) VMD TCL Scripting
2) Python numpy library
