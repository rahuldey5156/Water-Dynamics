import numpy as np
import MDAnalysis as mda

def calculate_residential_time(u, frame_needed):
    """Calculate residential time of oxygen atoms over specified frames."""
    results = []

    # Loop through each frame in the trajectory
    for ts in u.trajectory[:frame_needed - 1]:
        current_frame = ts.frame
        oxygens = u.select_atoms("name OW and around 5.5 protein")

        # Store indices and number of selected oxygen atoms for this frame
        selected_indices = oxygens.indices
        no_of_initial_oxygens = len(selected_indices)

        # Initialize an array to hold residential values for this frame
        residential = []

        # Second loop starting from the next frame
        for next_ts in u.trajectory[current_frame + 1:frame_needed]:
            current_oxygen = u.select_atoms("name OW and around 5.5 protein")

            # Find the intersection of indices between `selected_indices` and `next_indices`
            common_indices = np.intersect1d(selected_indices, current_oxygen.indices)
            
            if len(common_indices) == 0:
                break # No common indices; stop the loop

            #Get the no of common oxygen atoms
            no_of_current_oxygens = len(common_indices)

            div = no_of_current_oxygens / no_of_initial_oxygens  # Divide the no of oxygens
            residential.append(div) # Store the divided values

            # Update selected_indices to keep only common ones for the next iteration
            selected_indices = common_indices

        print(f"Frame {current_frame} completed")

        # Store the results for this frame
        results.append(residential)

    return results

def pad_and_calculate_means(results):
    """Pad results with NaN and calculate means."""
    max_length = max(len(arr) for arr in results)

    # Pad the arrays with np.nan to make them the same length
    padded_results = np.array([np.pad(arr, (0, max_length - len(arr)), constant_values=np.nan) for arr in results])

    # Calculate the mean along columns, ignoring nan values
    means = np.nanmean(padded_results, axis=0)
    
    return means

# Load pdb and xtc files
u = mda.Universe('test.pdb', 'test.xtc')
frame_needed = len(u.trajectory)

# Calculate squared differences
results = calculate_residential_time(u, frame_needed)

# Calculate means from the results
means = pad_and_calculate_means(results)

# Create an array with frame numbers and corresponding mean values
data_to_save = np.column_stack((np.arange(1, (len(means)+1)), means))

# Save the data to a file without a header
np.savetxt('continuous_residential_time.dat', data_to_save, fmt="%d %.4f", delimiter='\t')
