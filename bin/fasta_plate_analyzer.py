import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from collections import defaultdict
import os
import argparse

def parse_fasta_and_aggregate(fasta_filepath):
    """
    Parses a FASTA file and aggregates sequence counts per plate and well.
    Assumes ID format: ..._PlateID_WellNumber (e.g., ..._P1_080)
    WellNumber is an integer from 1 to 96.
    """
    plate_data = defaultdict(lambda: {"well_counts": defaultdict(int), "total_sequences": 0})

    try:
        with open(fasta_filepath, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]
                    parts = header.split('_')
                    if len(parts) >= 2:
                        plate_id = parts[-2]
                        try:
                            well_num_str = parts[-1]
                            well_num = int(well_num_str)
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
    for plate_id, data in sorted(plate_data.items()): # Sort by plate_id for consistent output
        num_used_wells = len(data["well_counts"]) 
        total_seqs = data["total_sequences"]
        print(f"Plate {plate_id}: Used Wells (with sequences) = {num_used_wells}, Total Sequences = {total_seqs}")

def plot_plate_layout(plate_id, well_counts, output_dir=".", numbering_order='column'):
    """
    Plots a 96-well plate layout for a single plate.
    numbering_order can be 'column' or 'row'.
    """
    fig, ax = plt.subplots(figsize=(12, 8.5)) 

    current_max_count = max(well_counts.values()) if well_counts else 0
    norm = mcolors.Normalize(vmin=0, vmax=max(1, current_max_count)) 
    
    cmap = plt.get_cmap('Reds')

    for well_num_1_indexed in range(1, 97): 
        if numbering_order == 'column':
            col_plot_idx = (well_num_1_indexed - 1) // 8
            row_plot_idx = (well_num_1_indexed - 1) % 8
        elif numbering_order == 'row':
            row_plot_idx = (well_num_1_indexed - 1) // 12
            col_plot_idx = (well_num_1_indexed - 1) % 12
        else:
            print(f"Error: Invalid numbering_order '{numbering_order}'. Defaulting to 'column'.")
            col_plot_idx = (well_num_1_indexed - 1) // 8
            row_plot_idx = (well_num_1_indexed - 1) % 8
            
        count = well_counts.get(well_num_1_indexed, 0)
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

def plot_overall_statistics_by_plate(plate_data, output_dir="."):
    """
    Plots a single bar chart summarizing statistics for all plates.
    X-axis: Plate IDs.
    Left Y-axis: Number of wells with sequences per plate.
    Right Y-axis: Total number of sequences per plate.
    """
    if not plate_data:
        print("No plate data available to generate overall statistics plot.")
        return

    # Sort data by plate ID for consistent plotting order
    sorted_plate_ids = sorted(plate_data.keys())
    
    num_wells_with_seqs_list = [len(plate_data[pid]["well_counts"]) for pid in sorted_plate_ids]
    total_sequences_list = [plate_data[pid]["total_sequences"] for pid in sorted_plate_ids]

    if not sorted_plate_ids: # Should not happen if plate_data is not empty, but good check
        print("No plates found in data for overall statistics plot.")
        return

    x_indices = np.arange(len(sorted_plate_ids))
    
    fig, ax1 = plt.subplots(figsize=(max(10, len(sorted_plate_ids) * 0.8), 7)) # Adjust width based on number of plates

    # Plotting Number of Wells with Sequences (Left Y-axis)
    bar_width = 0.35
    bars_wells = ax1.bar(x_indices - bar_width/2, num_wells_with_seqs_list, bar_width, color='dodgerblue', alpha=0.8, label='Wells with Sequences')
    ax1.set_xlabel('Plate ID', fontsize=12)
    ax1.set_ylabel('Number of Wells with Sequences', color='dodgerblue', fontsize=12)
    ax1.tick_params(axis='y', labelcolor='dodgerblue')
    ax1.set_xticks(x_indices)
    ax1.set_xticklabels(sorted_plate_ids, rotation=45, ha="right")
    ax1.grid(axis='y', linestyle='--', alpha=0.6)

    # Creating a second Y-axis for Total Sequences
    ax2 = ax1.twinx()
    # Plotting Total Sequences (Right Y-axis - can be bars or line)
    # Using bars for total sequences as well, slightly offset
    bars_seqs = ax2.bar(x_indices + bar_width/2, total_sequences_list, bar_width, color='tomato', alpha=0.8, label='Total Sequences')
    # Alternatively, for a line plot:
    # line_seqs = ax2.plot(x_indices, total_sequences_list, color='tomato', marker='o', linestyle='-', linewidth=2, markersize=6, label='Total Sequences')
    ax2.set_ylabel('Total Sequences', color='tomato', fontsize=12)
    ax2.tick_params(axis='y', labelcolor='tomato')
    
    # Add text labels on top of bars
    for bar in bars_wells:
        yval = bar.get_height()
        if yval > 0:
            ax1.text(bar.get_x() + bar.get_width()/2.0, yval + 0.01 * max(num_wells_with_seqs_list, default=1), 
                     int(yval), ha='center', va='bottom', fontsize=8, color='dodgerblue')

    for bar in bars_seqs:
        yval = bar.get_height()
        if yval > 0:
            ax2.text(bar.get_x() + bar.get_width()/2.0, yval + 0.01 * max(total_sequences_list, default=1), 
                     int(yval), ha='center', va='bottom', fontsize=8, color='tomato')

    # Add a combined legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    fig.legend(handles1 + handles2, labels1 + labels2, loc='upper center', bbox_to_anchor=(0.5, 0.98), ncol=2)

    fig.suptitle("Overall Plate Statistics", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.93]) 

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    output_filename = os.path.join(output_dir, "overall_plate_statistics_summary.pdf")
    try:
        plt.savefig(output_filename)
        print(f"Overall statistics plot saved to: {output_filename}")
    except Exception as e:
        print(f"Error saving overall statistics plot {output_filename}: {e}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Analyze FASTA sequences from 96-well plates and generate visualizations.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output_dir", default="plate_plots", help="Directory to save output plots (default: plate_plots).")
    parser.add_argument("-no", "--numbering_order", choices=['column', 'row'], default='column', 
                        help="Order to number wells for plotting individual plate layouts: 'column' or 'row'. Default: column.")
    
    args = parser.parse_args()

    fasta_file = args.input
    output_plot_dir = args.output_dir
    numbering_order = args.numbering_order # Still used for individual plate layouts
    
    if not os.path.exists(fasta_file):
        print(f"Error: Input FASTA file '{fasta_file}' not found. Please check the path.")
        return

    plate_data = parse_fasta_and_aggregate(fasta_file)

    if plate_data:
        print_statistics(plate_data)
        
        # Generate individual plate layout plots
        individual_plots_generated = False
        for plate_id, data in plate_data.items():
            if data["well_counts"]: 
                 plot_plate_layout(plate_id, data["well_counts"], output_dir=output_plot_dir, numbering_order=numbering_order)
                 individual_plots_generated = True
            else:
                print(f"Plate {plate_id} has no sequence data for layout plot, skipping.")
        
        if individual_plots_generated:
             print(f"\nIndividual plate layout plots saved in '{os.path.abspath(output_plot_dir)}'.")
        else:
            print("\nNo data to generate individual plate layout plots.")

        # Generate the overall statistics plot
        plot_overall_statistics_by_plate(plate_data, output_dir=output_plot_dir)

    else:
        print("\nNo data was processed from the FASTA file.")

if __name__ == "__main__":
    main()