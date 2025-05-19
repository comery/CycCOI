import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import sys

def read_sequences(filename):
    """Read sequence file"""
    sequences = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                sequences[parts[0]] = parts[1]
    return sequences

def edit_distance(seq1, seq2):
    """Calculate edit distance between two sequences"""
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize boundary conditions
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
    
    # Dynamic programming to calculate edit distance
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1]) + 1
    
    return dp[m][n]

def calculate_all_distances(sequences):
    """Calculate edit distances for all sequence pairs"""
    distances = []
    for seq1_id, seq2_id in combinations(sequences.keys(), 2):
        dist = edit_distance(sequences[seq1_id], sequences[seq2_id])
        distances.append(dist)
    return distances

def plot_histogram(distances):
    """Plot histogram of edit distances"""
    plt.figure(figsize=(10, 6))
    plt.hist(distances, bins='auto', edgecolor='black')
    plt.title('Distribution of Base Sequence Edit Distances')
    plt.xlabel('Edit Distance')
    plt.ylabel('Frequency')
    plt.grid(True, alpha=0.3)
    plt.savefig('edit_distance_histogram.pdf')  # 将文件扩展名改为.pdf
    plt.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 EditDistance.py <sequence_file_path>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    
    try:
        # Read input file
        sequences = read_sequences(input_file)
        
        # Calculate edit distances for all sequence pairs
        distances = calculate_all_distances(sequences)
        
        # Print statistics
        print(f"Minimum edit distance: {min(distances)}")
        print(f"Maximum edit distance: {max(distances)}")
        print(f"Average edit distance: {np.mean(distances):.2f}")
        
        # Plot histogram
        plot_histogram(distances)
        print("Histogram has been saved as 'edit_distance_histogram.pdf'")  # 更新输出消息
        
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()