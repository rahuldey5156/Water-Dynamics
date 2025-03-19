import numpy as np

# Paths to data files
file1 = 'donor.dat'
file2 = 'acceptor.dat'
file3 = 'donor_hydrogen.dat'

# Read data from each file
with open(file1, 'r') as f1, open(file2, 'r') as f2, open(file3, 'r') as f3:
    data1 = [list(map(int, line.split())) for line in f1]
    data2 = [list(map(int, line.split())) for line in f2]
    data3 = [list(map(int, line.split())) for line in f3]

# Create triplet arrays
triplet_arrays = []
for line1, line2, line3 in zip(data1, data2, data3):
    # Create a list of triplets for the current line
    triplet = [(v1, v2, v3) for v1, v2, v3 in zip(line1, line2, line3)]
    triplet_arrays.append(triplet)

# Initialize an array to hold the final values
results = []

for i in range(0, len(triplet_arrays)-1):
    current_triplet = triplet_arrays[i]
    no_of_current_triplet = len(triplet_arrays[i])

    # Initialize an array to hold protein-water hbond values for this frame
    protein_water_hbond = []

    # Second loop starting from the next frame
    for j in range(i+1, len(triplet_arrays)):

        # Find the intersection of triplets between 'current_triplet' and triplet in this frame
        next_triplet = triplet_arrays[j]
        common_triplet = set(current_triplet) & set(next_triplet)

        if len(common_triplet) == 0:
            break # No common triplets; stop the loop

        # Get the no of common triplets
        no_of_triplets = len(common_triplet)

        div = no_of_triplets / no_of_current_triplet # Divide the no of triplets
        protein_water_hbond.append(div) # Store the divided values

        # Update current_triplet to keep only common ones for the next iteration
        current_triplet = common_triplet

    # Store the results for this frame
    results.append(protein_water_hbond)

# Get the maximum length of the array    
max_length = max(len(arr) for arr in results)

# Pad the arrays with np.nan to make them the same length
padded_results = np.array([np.pad(arr, (0, max_length - len(arr)), constant_values=np.nan) for arr in results])

# Calculate the mean along columns, ignoring nan values
means = np.nanmean(padded_results, axis=0)

# Create an array with frame numbers and corresponding mean values
data_to_save = np.column_stack((np.arange(len(means)), means))

# Save the data to a file without a header
np.savetxt('continuous_protein_water_hbond.dat', data_to_save, fmt="%d %.4f", delimiter='\t')
