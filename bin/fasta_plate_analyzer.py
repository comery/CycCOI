import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from collections import defaultdict
import os
import argparse # Added for command-line arguments

def parse_fasta_and_aggregate(fasta_filepath):
    """
    Parses a FASTA file and aggregates sequence counts per plate and well.
    Assumes ID format: ..._PlateID_WellNumber (e.g., ..._P1_080)
    WellNumber is an integer from 1 to 96.
    """
    # Using defaultdict to simplify initialization
    # plate_data structure:
    # {
    #   "P1": { "well_counts": {1: count, 2: count, ...}, "total_sequences": total_count },
    #   "P2": { ... }
    # }
    plate_data = defaultdict(lambda: {"well_counts": defaultdict(int), "total_sequences": 0})

    try:
        with open(fasta_filepath, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]  # Remove '>' and newline
                    parts = header.split('_')
                    if len(parts) >= 2:
                        plate_id = parts[-2]
                        try:
                            well_num_str = parts[-1]
                            well_num = int(well_num_str)  # Well number should be 1-96
                            if not (1 <= well_num <= 96):
                                print(f"Warning: Invalid well number {well_num} in header '{header}' (should be 1-96). Skipping.")
                                continue
                        except ValueError:
                            print(f"Warning: Could not parse well number from '{parts[-1]}' in header '{header}'. Skipping.")
                            continue

                        plate_data[plate_id]["well_counts"][well_num] += 1
                        plate_data[plate_id]["total_sequences"] += 1
                    else:
                        print(f"Warning: Incorrect header format, not enough underscore-separated fields: '{header}'. Skipping.")
    except FileNotFoundError:
        print(f"Error: File {fasta_filepath} not found.")
        return None
    except Exception as e:
        print(f"An error occurred while reading or parsing the file: {e}")
        return None
        
    return plate_data

def print_statistics(plate_data):
    """
    Prints statistics for each plate.
    """
    if not plate_data:
        print("No data to report.")
        return
    
    print("\n--- Plate Statistics ---")
    for plate_id, data in plate_data.items():
        num_used_wells = len(data["well_counts"])  # Count how many wells are actually used
        total_seqs = data["total_sequences"]
        print(f"Plate {plate_id}: Used Wells = {num_used_wells}, Total Sequences = {total_seqs}")

def plot_plate_layout(plate_id, well_counts, output_dir=".", numbering_order='column'):
    """
    Plots a 96-well plate layout for a single plate.
    numbering_order can be 'column' or 'row'.
    """
    fig, ax = plt.subplots(figsize=(12, 8.5)) # Adjust figure size to accommodate title and colorbar

    current_max_count = max(well_counts.values()) if well_counts else 0
    norm = mcolors.Normalize(vmin=0, vmax=max(1, current_max_count)) 
    
    cmap = plt.get_cmap('Reds')

    for well_num_1_indexed in range(1, 97): # Iterate through logical wells 1 to 96
        # Determine the (row, column) on the plot based on numbering_order
        if numbering_order == 'column':
            # Column-major: A1, B1, ..., H1 (wells 1-8 for col 1), then A2, B2, ... (wells 9-16 for col 2)
            col_plot_idx = (well_num_1_indexed - 1) // 8  # 0-11 for plot columns
            row_plot_idx = (well_num_1_indexed - 1) % 8   # 0-7 for plot rows (A-H)
        elif numbering_order == 'row':
            # Row-major: A1, A2, ..., A12 (wells 1-12 for row A), then B1, B2, ... (wells 13-24 for row B)
            row_plot_idx = (well_num_1_indexed - 1) // 12 # 0-7 for plot rows (A-H)
            col_plot_idx = (well_num_1_indexed - 1) % 12  # 0-11 for plot columns
        else:
            # This case should ideally be prevented by argparse choices
            print(f"Error: Invalid numbering_order '{numbering_order}'. Defaulting to 'column'.")
            col_plot_idx = (well_num_1_indexed - 1) // 8
            row_plot_idx = (well_num_1_indexed - 1) % 8
            
        count = well_counts.get(well_num_1_indexed, 0) # Get sequence count for this logical well
        color = cmap(norm(count)) 

        circle = plt.Circle((col_plot_idx + 0.5, row_plot_idx + 0.5),
                              0.45,
                              facecolor=color,
                              edgecolor='grey',
                              linewidth=0.5)
        ax.add_patch(circle)

        if count > 0:
            luminance = 0.299*color[0] + 0.587*color[1] + 0.114*color[2]
            text_color = 'black' if luminance > 0.5 else 'white'
            
            ax.text(col_plot_idx + 0.5, row_plot_idx + 0.5, str(count),
                    ha='center', va='center', fontsize=8, color=text_color, weight='bold')

    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.set_xticks(np.arange(12) + 0.5)
    ax.set_xticklabels(np.arange(1, 13))
    ax.set_xlabel("Column")
    ax.set_yticks(np.arange(8) + 0.5)
    ax.set_yticklabels(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
    ax.set_ylabel("Row")
    ax.invert_yaxis()
    ax.set_title(f"Plate: {plate_id} - Sequence Counts per Well (Total: {sum(well_counts.values())}, Order: {numbering_order})", fontsize=14)
    ax.set_aspect('equal', adjustable='box')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label('Sequence Count', fontsize=10)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    output_filename = os.path.join(output_dir, f"plate_{plate_id}_visualization_{numbering_order}_order.pdf")
    try:
        plt.savefig(output_filename)
        print(f"Plot saved to: {output_filename}")
    except Exception as e:
        print(f"Error saving plot {output_filename}: {e}")
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Analyze FASTA sequences from 96-well plates and generate visualizations.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output_dir", default="plate_plots", help="Directory to save output plots (default: plate_plots).")
    parser.add_argument("-no", "--numbering_order", choices=['column', 'row'], default='column', 
                        help="Order to number wells for plotting: 'column' (A1,B1..H1,A2..) or 'row' (A1,A2..A12,B1..). Default: column.")
    
    args = parser.parse_args()

    fasta_file = args.input
    output_plot_dir = args.output_dir
    numbering_order = args.numbering_order
    
    if not os.path.exists(fasta_file):
        print(f"Error: Input FASTA file '{fasta_file}' not found. Please check the path.")
        return

    plate_data = parse_fasta_and_aggregate(fasta_file)

    if plate_data:
        print_statistics(plate_data)
        
        plots_generated = False
        for plate_id, data in plate_data.items():
            if data["well_counts"]: 
                 plot_plate_layout(plate_id, data["well_counts"], output_dir=output_plot_dir, numbering_order=numbering_order)
                 plots_generated = True
            else:
                print(f"Plate {plate_id} has no sequence data, skipping plot generation.")
        
        if plots_generated:
             print(f"\nAll plots have been saved in the '{os.path.abspath(output_plot_dir)}' directory.")
        else:
            print("\nNo valid plate data found to generate plots.")
    else:
        print("\nNo data was processed from the FASTA file.")

if __name__ == "__main__":
    main()