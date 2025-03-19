import numpy as np
import MDAnalysis as mda

def calculate_van_hove(u, frame_needed):
    """Calculate van hove of oxygen atom positions over specified frames."""
    results = []

    # Loop through each frame in the trajectory
    for ts in u.trajectory[:frame_needed - 250]:
        current_frame = ts.frame
        oxygens = u.select_atoms("name OW and around 5.5 protein")

        # Store indices and coordinates of selected oxygen atoms for this frame
        selected_indices = oxygens.indices
        coords_oxygens = oxygens.positions.copy()

        # Next frame taken after 200 frame
        next_frame = current_frame + 250
        # Access the target frame
        u.trajectory[next_frame]

        if next_frame > frame_needed:
            break # Next frame exceeds the total no of frames; stop the loop
        
        current_oxygen = u.select_atoms("name OW and around 5.5 protein")

        # Find the intersection of indices between `selected_indices` and `next_indices`
        common_indices = np.intersect1d(selected_indices, current_oxygen.indices)
            
        if len(common_indices) == 0:
            continue  # Move to the next frame

        #Select the common oxygen atoms and get their positions
        filtered_oxygen = u.select_atoms(f"index {' '.join(map(str, common_indices))}")
        current_coords_oxygens = filtered_oxygen.positions

        # Get original coordinates based on common indices
        original_coords_oxygens = coords_oxygens[np.isin(selected_indices, common_indices)]

        diff = current_coords_oxygens - original_coords_oxygens  # Subtract coordinates
        squared_diff = (np.sum(diff ** 2, axis=1))  # Square root and sum differences

        # Take square root of summed squared differences
        displacements = np.sqrt(squared_diff)

        #diffs = np.mean(displacements)  # Average over all oxygens

        for idx, values in zip(common_indices, displacements):
            results.append((idx, values))

        # Append the current_frame, next_frame, and average difference to results
        #results.append([current_frame, next_frame, diffs])

        print(f"Frame {current_frame} , {next_frame} completed")

    return np.array(results)


# Load pdb and xtc files
u = mda.Universe('van_hove.pdb', 'van_hove.xtc')
frame_needed = 5000

# Calculate van hove
results = calculate_van_hove(u, frame_needed)

# Save the data to a file without a header
np.savetxt('van_hove_5000.dat', results, fmt="%d %.4f", delimiter='\t')
