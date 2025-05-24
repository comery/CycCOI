import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from collections import defaultdict, Counter
import subprocess
import os
import shutil
import tempfile

def parse_fasta_and_group(fasta_filepath):
    """
    Parses a FASTA file and groups sequences by plate and well.
    ID format: ..._PlateID_WellID
    Returns a dictionary: {(plate_id, well_id): [SeqRecord, ...]}
    """
    grouped_sequences = defaultdict(list)
    try:
        for record in SeqIO.parse(fasta_filepath, "fasta"):
            try:
                parts = record.id.split('_')
                if len(parts) >= 2:
                    plate_id = parts[-2]
                    well_id_str = parts[-1]
                    # Attempt to convert well_id to int for consistent sorting/handling if needed, but store as str
                    # well_id = int(well_id_str) # We'll keep it as string as per FASTA ID
                    grouped_sequences[(plate_id, well_id_str)].append(record)
                else:
                    print(f"Warning: Could not parse plate/well from ID: {record.id}. Skipping sequence.")
            except Exception as e:
                print(f"Warning: Error processing ID {record.id}: {e}. Skipping sequence.")
        if not grouped_sequences:
            print("Warning: No sequences were successfully grouped. Check FASTA IDs and format.")
        return grouped_sequences
    except FileNotFoundError:
        print(f"Error: Input FASTA file {fasta_filepath} not found.")
        return None
    except Exception as e:
        print(f"Error reading FASTA file {fasta_filepath}: {e}")
        return None

def run_vsearch_cluster(sequences, vsearch_path, id_threshold, temp_dir_sample):
    """
    Runs vsearch to cluster sequences for a single sample.
    sequences: list of SeqRecord objects for the current sample.
    vsearch_path: path to vsearch executable.
    id_threshold: clustering identity threshold (e.g., 0.97).
    temp_dir_sample: temporary directory for this sample's files.
    Returns a list of clusters, where each cluster is a list of SeqRecord objects.
    Or path to centroid file and uclust file for further parsing.
    """
    if not sequences:
        return []

    sample_fasta_path = os.path.join(temp_dir_sample, "sample_input.fasta")
    SeqIO.write(sequences, sample_fasta_path, "fasta")

    centroids_path = os.path.join(temp_dir_sample, "centroids.fasta")
    uclust_path = os.path.join(temp_dir_sample, "clusters.uc")

    # VSEARCH command for clustering
    # Using --cluster_fast as it's generally good for OTU picking
    # --iddef 1 might be needed if sequences are very similar and you want exact matches for centroids
    cmd = [
        vsearch_path,
        "--cluster_fast", sample_fasta_path,
        "--id", str(id_threshold),
        "--centroids", centroids_path,
        "--uc", uclust_path,
        "--quiet", # Add quiet to reduce verbosity, remove for debugging
        # "--sizeout" # Optionally, to include abundance in centroid headers
    ]

    print(f"Running vsearch for sample: {'_'.join(sequences[0].id.split('_')[-2:])} ...")
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # print(f"VSEARCH stdout:\n{process.stdout}") # For debugging
        # print(f"VSEARCH stderr:\n{process.stderr}") # For debugging
    except subprocess.CalledProcessError as e: # 捕获vsearch命令执行失败的异常
        print(f"Error running vsearch for {sample_fasta_path}:")
        print(f"Command: {' '.join(e.cmd)}")
        print(f"Return code: {e.returncode}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return None, None # Indicate error
    except FileNotFoundError:
        print(f"Error: vsearch executable not found at {vsearch_path}. Please check the path.")
        return None, None

    return centroids_path, uclust_path

def parse_uclust_and_group_sequences(uclust_path, all_sample_sequences_dict):
    """
    Parses a .uc file to group sequences into clusters based on centroid.
    all_sample_sequences_dict: A dictionary mapping original sequence IDs to SeqRecord objects for the current sample.
    Returns a dictionary: {centroid_id: [SeqRecord_member1, SeqRecord_member2, ...]}.
    """
    clusters = defaultdict(list)
    try:
        with open(uclust_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                record_type = parts[0]
                # query_id = parts[8]
                # target_id = parts[9]

                if record_type == 'H': # Hit line, query is a member of target (centroid)
                    member_id = parts[8]
                    centroid_id = parts[9]
                    if member_id in all_sample_sequences_dict:
                        clusters[centroid_id].append(all_sample_sequences_dict[member_id])
                    else:
                        print(f"Warning: Sequence ID {member_id} from .uc file not found in original sample sequences.")
                elif record_type == 'S': # Centroid declaration by vsearch (not always a hit to itself)
                    # This line declares a sequence as a centroid. We ensure it's added to its own cluster.
                    # However, --cluster_fast centroids are actual sequences from input.
                    # We can also get centroids directly from centroids.fasta and then fill with members from 'H' lines.
                    pass # Handled by 'C' or by ensuring centroids from centroids.fasta are processed
                elif record_type == 'C': # Cluster line (alternative way to get cluster size, not used here for members)
                    pass 

        # Ensure centroids themselves are part of their clusters if not already added via 'H' to self
        # This is more robustly handled by reading centroids.fasta and then populating with .uc hits
        # For now, this structure assumes .uc provides enough info or centroids.fasta is primary source

    except FileNotFoundError:
        print(f"Error: uclust file {uclust_path} not found.")
        return {}
    except Exception as e:
        print(f"Error parsing uclust file {uclust_path}: {e}")
        return {}
    return clusters

def run_mafft_msa(sequences_for_msa, mafft_path, temp_dir_cluster):
    """
    Runs mafft to perform MSA on a list of sequences (a cluster).
    sequences_for_msa: list of SeqRecord objects.
    mafft_path: path to mafft executable.
    temp_dir_cluster: temporary directory for this cluster's files.
    Returns path to the MSA file (FASTA format) or None on error.
    """
    if not sequences_for_msa or len(sequences_for_msa) < 2: # MSA needs at least 2 sequences
        # If only 1 sequence, it is its own consensus, no MSA needed.
        # We can return a path to a file containing just this one sequence, or handle it upstream.
        if sequences_for_msa:
            single_seq_path = os.path.join(temp_dir_cluster, "single_sequence.fasta")
            SeqIO.write(sequences_for_msa, single_seq_path, "fasta")
            return single_seq_path # This will be read as an "alignment" of one sequence
        return None

    cluster_fasta_path = os.path.join(temp_dir_cluster, "cluster_input.fasta")
    SeqIO.write(sequences_for_msa, cluster_fasta_path, "fasta")

    msa_output_path = os.path.join(temp_dir_cluster, "aligned.fasta")

    # MAFFT command (basic auto strategy)
    # --auto: Automatically selects an appropriate strategy
    # --quiet: Suppresses verbose output
    cmd = [
        mafft_path,
        "--auto",
        "--quiet",
        cluster_fasta_path
    ]

    print(f"Running mafft for cluster ({len(sequences_for_msa)} sequences) ...") 
    try:
        with open(msa_output_path, "w") as msa_file:
            process = subprocess.run(cmd, stdout=msa_file, stderr=subprocess.PIPE, text=True, check=True)
        # print(f"MAFFT stderr:\n{process.stderr}") # For debugging
    except subprocess.CalledProcessError as e:
        print(f"Error running mafft for {cluster_fasta_path}:")
        print(f"Command: {' '.join(e.cmd)}")
        print(f"Return code: {e.returncode}")
        # Stderr might be in e.stderr if captured, or check mafft's own logging
        # print(f"Stdout: {e.stdout}") # stdout is redirected to file
        print(f"Stderr: {e.stderr}")
        return None
    except FileNotFoundError:
        print(f"Error: mafft executable not found at {mafft_path}. Please check the path.")
        return None

    return msa_output_path

def get_consensus_from_msa(msa_filepath, original_centroid_id, abundance, plate_id, well_id, cluster_num):
    """Generates a consensus sequence from an MSA file with custom consensus logic."""
    try:
        alignment = AlignIO.read(msa_filepath, "fasta")
        if not alignment:
            print(f"Warning: MSA file {msa_filepath} is empty or unreadable by AlignIO.")
            return None
        
        if len(alignment) == 1: # If MSA was of a single sequence
            consensus_seq = alignment[0].seq
            consensus_id = f"{plate_id}_{well_id}_Cluster{cluster_num}_Abundance{abundance}_CentroidIsSelf_{alignment[0].id}"
            return SeqRecord(consensus_seq, id=consensus_id, description="Consensus from single-sequence cluster")

        # 自定义consensus逻辑，替代原来的dumb_consensus
        con_len = alignment.get_alignment_length()
        consensus = ''
        
        # 遍历每个位置
        for n in range(con_len):
            # 统计碱基频率
            base_counts = {}
            valid_bases = 0
            
            for record in alignment:
                if n < len(record.seq):
                    base = record.seq[n]
                    if base != '-' and base != '.':
                        valid_bases += 1
                        if base not in base_counts:
                            base_counts[base] = 1
                        else:
                            base_counts[base] += 1
            
            # 如果没有有效碱基，跳过（这会在最终序列中删除该位点）
            if valid_bases == 0:
                continue
                
            # 找出频率最高的两个碱基
            sorted_bases = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
            
            # 如果只有一种碱基
            if len(sorted_bases) == 1:
                consensus += sorted_bases[0][0]
            # 如果有多种碱基
            else:
                top1_base, top1_count = sorted_bases[0]
                top2_base, top2_count = sorted_bases[1] if len(sorted_bases) > 1 else (None, 0)
                
                # 如果top1频率超过50%
                if top1_count > valid_bases * 0.5:
                    consensus += top1_base
                # 如果top1和top2频率相同
                elif top1_count == top2_count:
                    # 这里可以实现简并碱基的逻辑
                    # 例如：如果是A和G，使用R；如果是C和T，使用Y等
                    # 简单示例：使用IUPAC简并碱基表示
                    bases_set = set([b[0] for b in sorted_bases if b[1] == top1_count])
                    if bases_set == {'A', 'G'}:
                        consensus += 'R'  # A或G
                    elif bases_set == {'C', 'T'}:
                        consensus += 'Y'  # C或T
                    elif bases_set == {'G', 'C'}:
                        consensus += 'S'  # G或C
                    elif bases_set == {'A', 'T'}:
                        consensus += 'W'  # A或T
                    elif bases_set == {'G', 'T'}:
                        consensus += 'K'  # G或T
                    elif bases_set == {'A', 'C'}:
                        consensus += 'M'  # A或C
                    else:
                        # 其他组合使用N
                        consensus += 'N'
                else:
                    consensus += top1_base
        
        # 创建序列对象
        consensus_seq = Seq(consensus)
        consensus_id = f"{plate_id}_{well_id}_Cluster{cluster_num}_Abundance{abundance}_Centroid_{original_centroid_id}"
        return SeqRecord(consensus_seq, id=consensus_id, description="Custom consensus sequence")
        
    except ValueError as e:
        print(f"Error generating consensus for {msa_filepath}: {e}. This might be due to an empty or invalid MSA file.")
        return None
    except Exception as e:
        print(f"Unexpected error generating consensus for {msa_filepath}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Process FASTA sequences: group, cluster, MSA, and generate consensus.")
    parser.add_argument("-i", "--input_fasta", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save output consensus FASTA files and temp files.")
    parser.add_argument("--vsearch_path", default="vsearch", help="Path to the vsearch executable (default: 'vsearch').")
    parser.add_argument("--mafft_path", default="mafft", help="Path to the mafft executable (default: 'mafft').")
    parser.add_argument("--id_threshold", type=float, default=0.97, help="VSEARCH clustering identity threshold (default: 0.97).")
    parser.add_argument("--min_abundance", type=int, default=2, help="Minimum number of sequences in a cluster to be processed (default: 2).")
    parser.add_argument("--keep_temp_files", action="store_true", help="Keep temporary files after execution (for debugging).")

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Create a main temporary directory within the output directory
    main_temp_dir = os.path.join(args.output_dir, "temp_processing_files")
    if os.path.exists(main_temp_dir):
        shutil.rmtree(main_temp_dir) # Clean up from previous run if necessary
    os.makedirs(main_temp_dir)

    print(f"Parsing FASTA file: {args.input_fasta}")
    grouped_sequences_by_sample = parse_fasta_and_group(args.input_fasta)

    if not grouped_sequences_by_sample:
        print("Exiting due to errors in FASTA parsing.")
        if not args.keep_temp_files:
            shutil.rmtree(main_temp_dir)
        return

    all_final_consensus_seqs = []

    for (plate_id, well_id), sequences in grouped_sequences_by_sample.items():
        sample_id_str = f"{plate_id}_{well_id}"
        print(f"\nProcessing sample: {sample_id_str} ({len(sequences)} sequences)")

        # Create a temporary directory for this specific sample's processing
        temp_dir_sample = tempfile.mkdtemp(dir=main_temp_dir, prefix=f"sample_{sample_id_str}_")
        
        # Create a dictionary of sequence ID to SeqRecord for easy lookup within the sample
        sample_sequences_dict = {seq.id: seq for seq in sequences}

        centroids_fasta_path, uclust_filepath = run_vsearch_cluster(sequences, args.vsearch_path, args.id_threshold, temp_dir_sample)

        if not centroids_fasta_path or not uclust_filepath:
            print(f"Skipping sample {sample_id_str} due to vsearch clustering error.")
            continue
        
        # Read centroids directly from centroids.fasta
        # This gives us the centroid sequences and their vsearch-generated IDs
        try:
            centroids_records = {rec.id: rec for rec in SeqIO.parse(centroids_fasta_path, "fasta")}
        except FileNotFoundError:
            print(f"Error: Centroids file {centroids_fasta_path} not found for sample {sample_id_str}. Skipping.")
            continue
        
        # Parse .uc file to get cluster members for each centroid
        # clusters_by_centroid_id will be: {centroid_id_from_uc: [SeqRecord_member1, SeqRecord_member2, ...]}
        clusters_by_centroid_id = parse_uclust_and_group_sequences(uclust_filepath, sample_sequences_dict)
        
        # Refine clusters: ensure centroid sequence is first and map to actual centroid SeqRecords
        # final_clusters will be: {centroid_seq_record_id: [SeqRecord_centroid, SeqRecord_member1, ...]}
        final_clusters_for_sample = defaultdict(list)
        for vsearch_centroid_id, member_records in clusters_by_centroid_id.items():
            if vsearch_centroid_id in centroids_records:
                centroid_seq_record = centroids_records[vsearch_centroid_id]
                # Add centroid first, then other members (if not already the centroid itself)
                current_cluster_seqs = [centroid_seq_record] 
                for member_rec in member_records:
                    if member_rec.id != centroid_seq_record.id: # Avoid duplicating the centroid
                        current_cluster_seqs.append(member_rec)
                final_clusters_for_sample[centroid_seq_record.id] = current_cluster_seqs
            else:
                # This case might happen if a centroid in .uc is not in centroids.fasta, which is unusual
                # Or if a sequence hits a centroid that wasn't itself a sequence in the input (e.g. --cluster_smallmem)
                # For --cluster_fast, centroids are from input, so this should be rare.
                print(f"Warning: Centroid ID {vsearch_centroid_id} from .uc not found in {centroids_fasta_path} for sample {sample_id_str}.")

        # Also, any centroid from centroids.fasta that didn't get members from .uc (a cluster of 1)
        for cid, crec in centroids_records.items():
            if cid not in final_clusters_for_sample:
                 final_clusters_for_sample[cid] = [crec]

        if not final_clusters_for_sample:
            print(f"No clusters found or parsed for sample {sample_id_str}.")
            continue

        print(f"Found {len(final_clusters_for_sample)} initial clusters for sample {sample_id_str}.")

        cluster_counter = 0
        for centroid_original_id, cluster_sequences in sorted(final_clusters_for_sample.items(), key=lambda item: len(item[1]), reverse=True):
            abundance = len(cluster_sequences)
            if abundance < args.min_abundance:
                print(f"  Cluster with centroid {centroid_original_id} (Abundance: {abundance}) below threshold {args.min_abundance}. Skipping.")
                continue
            
            cluster_counter += 1
            print(f"  Processing Cluster {cluster_counter} (Centroid: {centroid_original_id}, Abundance: {abundance})")

            # Create a temporary directory for this specific cluster's MSA
            temp_dir_cluster = tempfile.mkdtemp(dir=temp_dir_sample, prefix=f"cluster_{cluster_counter}_")

            msa_filepath = run_mafft_msa(cluster_sequences, args.mafft_path, temp_dir_cluster)

            if msa_filepath:
                consensus_record = get_consensus_from_msa(msa_filepath, centroid_original_id, abundance, plate_id, well_id, cluster_counter)
                if consensus_record:
                    all_final_consensus_seqs.append(consensus_record)
            else:
                print(f"    MSA failed or no sequences for cluster {cluster_counter}. No consensus generated.")

    # Write all consensus sequences to a single output file
    output_consensus_fasta = os.path.join(args.output_dir, "all_consensus_sequences.fasta")
    if all_final_consensus_seqs:
        SeqIO.write(all_final_consensus_seqs, output_consensus_fasta, "fasta")
        print(f"\nSuccessfully generated {len(all_final_consensus_seqs)} consensus sequences.")
        print(f"Output saved to: {output_consensus_fasta}")
    else:
        print("\nNo consensus sequences were generated.")

    # Clean up main temporary directory if not keeping files
    if not args.keep_temp_files:
        print(f"Cleaning up temporary directory: {main_temp_dir}")
        shutil.rmtree(main_temp_dir)
    else:
        print(f"Temporary files kept at: {main_temp_dir}")

    print("Processing complete.")

if __name__ == "__main__":
    main()