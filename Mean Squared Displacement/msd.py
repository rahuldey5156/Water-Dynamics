import numpy as np
import MDAnalysis as mda

def calculate_squared_diffs(u, frame_needed):
    """Calculate mean squared differences of oxygen atom positions over specified frames."""
    results = []

    # Loop through each frame in the trajectory
    for ts in u.trajectory[:frame_needed - 1]:
        current_frame = ts.frame
        oxygens = u.select_atoms("name OW and around 5.5 protein")

        # Store indices and coordinates of selected oxygen atoms for this frame
        selected_indices = oxygens.indices
        coords_oxygens = oxygens.positions.copy()

        # Initialize an array to hold squared differences for this frame
        squared_diffs = []

        # Second loop starting from the next frame
        for next_ts in u.trajectory[current_frame + 1:frame_needed]:
            current_oxygen = u.select_atoms("name OW and around 5.5 protein")

            # Find the intersection of indices between `selected_indices` and `next_indices`
            common_indices = np.intersect1d(selected_indices, current_oxygen.indices)
            
            if len(common_indices) == 0:
                break # No common indices; stop the loop

            #Select the common oxygen atoms and get their positions
            filtered_oxygen = u.select_atoms(f"index {' '.join(map(str, common_indices))}")
            current_coords_oxygens = filtered_oxygen.positions

            # Get original coordinates based on common indices
            original_coords_oxygens = coords_oxygens[np.isin(oxygens.indices, common_indices)]

            diff = current_coords_oxygens - original_coords_oxygens  # Subtract coordinates
            squared_diff = np.sum(diff ** 2, axis=1)  # Square and sum differences

            squared_diffs.append(np.mean(squared_diff))  # Average over all oxygens

            # Update selected_indices to keep only common ones for the next iteration
            selected_indices = common_indices

        print(f"Frame {current_frame} completed")

        # Store the results for this frame
        results.append(squared_diffs)

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
u = mda.Universe('msd.pdb', 'msd.xtc')
frame_needed = len(u.trajectory)

# Calculate squared differences
results = calculate_squared_diffs(u, frame_needed)

# Calculate means from the results
means = pad_and_calculate_means(results)

# Create an array with frame numbers and corresponding mean values
data_to_save = np.column_stack((np.arange(len(means)), means))

# Save the data to a file without a header
np.savetxt('msd.dat', data_to_save, fmt="%d %.4f", delimiter='\t')
