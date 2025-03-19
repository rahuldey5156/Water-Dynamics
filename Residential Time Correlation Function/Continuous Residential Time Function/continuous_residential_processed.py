import numpy as np

# Read data from the file
data = []

# Path to data file
with open('hydration_water_index.dat', 'r') as file:
    # Skip the first two lines
    next(file)
    next(file)
    
    # Read the remaining lines and convert each line into an array of integers
    for line in file:
        array = list(map(int, line.split()))
        data.append(array)
    
# Initialize an array to hold the final values
results = []

for i in range(0, len(data)-1):
    current_index = data[i]
    no_of_current_index = len(data[i])

    # Initialize an array to hold residential values for this frame
    residential = []

    # Second loop starting from the next frame
    for j in range(i+1, len(data)):

        # Find the intersection of indices between 'current_index' and index in this frame
        next_indices = data[j]
        common_indices = set(current_index) & set(next_indices)

        if len(common_indices) == 0:
            break # No common indices; stop the loop

        # Get the no of common indices
        no_of_indices = len(common_indices)

        div = no_of_indices / no_of_current_index # Divide the no of indices
        residential.append(div) # Store the divided values

        # Update current_index to keep only common ones for the next iteration
        current_index = common_indices

    print(f"Frame {i} completed")

    # Store the results for this frame
    results.append(residential)

# Get the maximum length of the array    
max_length = max(len(arr) for arr in results)

# Pad the arrays with np.nan to make them the same length
padded_results = np.array([np.pad(arr, (0, max_length - len(arr)), constant_values=np.nan) for arr in results])

# Calculate the mean along columns, ignoring nan values
means = np.nanmean(padded_results, axis=0)

# Create an array with frame numbers and corresponding mean values
data_to_save = np.column_stack((np.arange(1, (len(means)+1)), means))

# Save the data to a file without a header
np.savetxt('continuous_residential_time.dat', data_to_save, fmt="%d %.4f", delimiter='\t')
