In this case, the water molecule must remain in contact with the protein continuously without interruption. If the water molecule leaves the hydration shell and returns at a later time, it is considered to have left, and the correlation is broken.

There are 2 ways of calculating Continuous Residential Time Correlation Function:
1) Run the Python Script named "continuous_residential_time_correlation_funcrion.py"
2) First run the TCL script named "protein_hydration_water_index.tcl". Then run the Python script named "continuous_residential_processed.py"

METHOD 1:

Takes longer time to complete. It will create output file named "continuous_residential_time.dat".
Requirements:
1) MDAnalysis
2) Numpy

METHOD 2:

Takes quicker time to complete. TCL script creates output file named "hydration_water_index.dat" and then the python script creates final output file named "continuous_residential_time.dat".
Requirements:
1) VMD TCL Scripting
2) Python numpy library
