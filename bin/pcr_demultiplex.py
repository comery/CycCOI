from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import defaultdict
import argparse
import sys
import os
import time
from multiprocessing import Pool
import multiprocessing
from icecream import ic

class PCRDemultiplexer:
    def __init__(self, forward_primer, reverse_primer, plate_indices, well_indices,
                 primer_max_mismatch=3, primer_max_indel=1,
                 index_max_mismatch=4, index_max_indel=1):
        """
        Initialize demultiplexer
        """
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.forward_primer_rc = str(Seq(forward_primer).reverse_complement())
        self.reverse_primer_rc = str(Seq(reverse_primer).reverse_complement())
        self.plate_indices = plate_indices
        self.well_indices = well_indices
        self.index_length = 13
        self.primer_max_mismatch = primer_max_mismatch
        self.primer_max_indel = primer_max_indel
        self.index_max_mismatch = index_max_mismatch
        self.index_max_indel = index_max_indel

        # 预计算所有可能需要的反向互补序列
        self.index_rc_cache = {}
        for _, index_seq in plate_indices + well_indices:
            self.index_rc_cache[index_seq] = self.reverse_complement(index_seq)

        # 创建index序列的精确匹配字典
        self.plate_dict = {}
        self.well_dict = {}

        for plate_id, index_seq in plate_indices:
            self.plate_dict[index_seq] = plate_id
            self.plate_dict[self.reverse_complement(index_seq)] = plate_id

        for well_id, index_seq in well_indices:
            self.well_dict[index_seq] = well_id
            self.well_dict[self.reverse_complement(index_seq)] = well_id

    def reverse_complement(self, sequence):
        """
        Calculate the reverse complement sequence
        """
        return str(Seq(sequence).reverse_complement())


    def smith_waterman(self, seq1, seq2, match_score=2, mismatch_score=-1, indel_score=-2):
        """
        Perform local sequence alignment using Smith-Waterman algorithm
        返回 (max_score, align_start, align_end)
        """
        m, n = len(seq1), len(seq2)
        score_matrix = np.zeros((m+1, n+1))
        pointer_matrix = np.zeros((m+1, n+1), dtype=int)  # 0: none, 1: diag, 2: up, 3: left

        max_score = 0
        max_pos = (0, 0)

        for i in range(1, m+1):
            for j in range(1, n+1):
                match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
                delete = score_matrix[i-1][j] + indel_score
                insert = score_matrix[i][j-1] + indel_score
                scores = [0, match, delete, insert]
                best = np.argmax(scores)
                score_matrix[i][j] = scores[best]
                pointer_matrix[i][j] = best
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)

        # 回溯获得起点
        i, j = max_pos
        align_end = i
        while pointer_matrix[i][j] != 0 and i > 0 and j > 0:
            if pointer_matrix[i][j] == 1:
                i -= 1
                j -= 1
            elif pointer_matrix[i][j] == 2:
                i -= 1
            elif pointer_matrix[i][j] == 3:
                j -= 1
        align_start = i

        return max_score, align_start, align_end


    def find_index_in_region(self, search_region, index_sequences, is_plate=True):
        """
        在搜索区域中寻找最佳匹配的索引序列
        使用Smith-Waterman局部比对算法
        """
        best_id = None
        best_score = float('-inf')
        best_seq = None
        best_mismatches = 0
        best_indels = 0

        max_mismatch = self.index_max_mismatch
        max_indel = self.index_max_indel

        for id_value, index_seq in index_sequences:
            score, start_pos, end_pos = self.smith_waterman(search_region, index_seq)
            ic(search_region, id_value, index_seq, score, end_pos)
            score_rc, start_pos_rc, end_pos_rc = self.smith_waterman(search_region, self.index_rc_cache[index_seq])

            threshold = len(index_seq) * 2 - (max_mismatch * 1 + max_indel * 2)

            if score > score_rc and score >= threshold:
                if score > best_score:
                    best_score = score
                    best_id = id_value
                    start_pos = max(0, start_pos)
                    best_seq = search_region[start_pos:end_pos]
                    mismatches = len(index_seq) * 2 - score
                    indels = mismatches // 2
                    mismatches = mismatches - (indels * 2)
                    best_mismatches = mismatches
                    best_indels = indels
            elif score_rc >= threshold:
                if score_rc > best_score:
                    best_score = score_rc
                    best_id = id_value
                    start_pos_rc = max(0, start_pos_rc)
                    best_seq = search_region[start_pos_rc:end_pos_rc]
                    mismatches = len(index_seq) * 2 - score_rc
                    indels = mismatches // 2
                    mismatches = mismatches - (indels * 2)
                    best_mismatches = mismatches
                    best_indels = indels
        return best_id, best_seq, best_score, best_mismatches, best_indels

    def process_sequence(self, seq_record):
        sequence = str(seq_record.seq)
        
        # Use a list to collect log information instead of printing directly
        logs = []
        logs.append("==========================================")
        logs.append(f"Processing sequence: {seq_record.id}")

        # Initialize plate_id and well_id as None
        plate_id = None
        well_id = None

        # Calculate the minimum required index length (considering allowed indels)
        min_index_length = self.index_length - self.index_max_indel

        # Forward-reverse structure
        fp_pos = -1
        rp_rc_pos = -1

        # Try exact match for forward primer
        try:
            fp_pos = sequence.index(self.forward_primer)
            logs.append(f"Forward primer: EXACT match at position {fp_pos}")
        except ValueError:
            max_score, fp_pos, end_pos = self.smith_waterman(sequence, self.forward_primer)
            #ic(max_score, end_pos)
            if max_score >= len(self.forward_primer) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_indel * 2):
                mismatches = len(self.forward_primer) * 2 - max_score
                indels = mismatches // 2
                mismatches = mismatches - (indels * 2)
                ic(max_score, mismatches, indels)
                logs.append(f"Forward primer: FUZZY match at position {fp_pos}")
                logs.append(f"Score: {max_score}, Mismatches: {mismatches}, Indels: {indels}")
            else:
                logs.append("Forward primer: NO match found")

        # If forward primer found, try to find reverse primer
        if fp_pos != -1:
            try:
                rp_rc_pos = sequence.index(self.reverse_primer_rc)
                logs.append(f"Reverse primer RC: EXACT match at position {rp_rc_pos}")
            except ValueError:
                max_score, rp_rc_pos, end_pos = self.smith_waterman(sequence, self.reverse_primer_rc)
                if max_score >= len(self.reverse_primer_rc) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_indel * 2):
                    mismatches = len(self.reverse_primer_rc) * 2 - max_score
                    indels = mismatches // 2
                    mismatches = mismatches - (indels * 2)
                    logs.append(f"Reverse primer RC: FUZZY match at position {rp_rc_pos}")
                    logs.append(f"Score: {max_score}, Mismatches: {mismatches}, Indels: {indels}")
                else:
                    logs.append("Reverse primer RC: NO match found")

        # If both primers found, try to find index
        if fp_pos != -1 and rp_rc_pos != -1:
            # 正向结构：plate index 紧挨着 forward primer 前面，well index 紧挨着 reverse primer RC 后面
            forward_index_search_end = fp_pos
            forward_index_search_start = max(0, forward_index_search_end - (self.index_length + self.index_max_indel))
            forward_index_search_region = sequence[forward_index_search_start:forward_index_search_end]

            reverse_index_search_start = rp_rc_pos + len(self.reverse_primer_rc)
            reverse_index_search_end = min(len(sequence), reverse_index_search_start + self.index_length + self.index_max_indel)
            reverse_index_search_region = sequence[reverse_index_search_start:reverse_index_search_end]

            # Check if search region is long enough
            if len(forward_index_search_region) < min_index_length:
                logs.append("Sequence too short for forward index extraction")
                return None

            if len(reverse_index_search_region) < min_index_length:
                logs.append("Sequence too short for reverse index extraction")
                return None
            
            logs.append(f"Forward index search region: {forward_index_search_region}")
            logs.append(f"Reverse index search region: {reverse_index_search_region}")
            
            # Find best matching index sequence in the search region
            plate_id, plate_seq, plate_score, plate_mismatches, plate_indels = self.find_index_in_region(
                forward_index_search_region, self.plate_indices, is_plate=True)
            
            if plate_id is not None:
                logs.append(f"Plate index: FUZZY match -> {plate_id}")
                logs.append(f"Matched sequence: {plate_seq}")
                logs.append(f"Score: {plate_score}, Mismatches: {plate_mismatches}, Indels: {plate_indels}")
            else:
                for seq in forward_index_search_region.split():
                    if seq in self.plate_dict:
                        plate_id = self.plate_dict[seq]
                        logs.append(f"Plate index: EXACT match -> {plate_id}")
                        break
                
                if plate_id is None:
                    logs.append("Plate index: NO match found")
                    return None
            
            if plate_id is not None:
                well_id, well_seq, well_score, well_mismatches, well_indels = self.find_index_in_region(
                    reverse_index_search_region, self.well_indices, is_plate=False)
                
                if well_id is not None:
                    logs.append(f"Well index: FUZZY match -> {well_id}")
                    logs.append(f"Matched sequence: {well_seq}")
                    logs.append(f"Score: {well_score}, Mismatches: {well_mismatches}, Indels: {well_indels}")
                else:
                    for seq in reverse_index_search_region.split():
                        if seq in self.well_dict:
                            well_id = self.well_dict[seq]
                            logs.append(f"Well index: EXACT match -> {well_id}")
                            break
                    
                    if well_id is None:
                        logs.append("Well index: NO match found")
            
            if plate_id and well_id:
                return {
                    'plate_id': plate_id,
                    'well_id': well_id,
                    'sequence': sequence[forward_index_search_start:reverse_index_search_end],
                    'read_id': seq_record.id,
                    'logs': logs
                }
            
            return None

        # Reverse-forward structure
        rp_pos = -1
        fp_rc_pos = -1

        try:
            rp_pos = sequence.index(self.reverse_primer)
            logs.append(f"Reverse primer: EXACT match at position {rp_pos}")
        except ValueError:
            max_score, rp_pos, end_pos = self.smith_waterman(sequence, self.reverse_primer)
            if max_score >= len(self.reverse_primer) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_indel * 2):
                mismatches = len(self.reverse_primer) * 2 - max_score
                indels = mismatches // 2
                mismatches = mismatches - (indels * 2)
                logs.append(f"Reverse primer: FUZZY match at position {rp_pos}")
                logs.append(f"Score: {max_score}, Mismatches: {mismatches}, Indels: {indels}")
            else:
                logs.append("Reverse primer: NO match found")

        if rp_pos != -1:
            try:
                fp_rc_pos = sequence.index(self.forward_primer_rc)
                logs.append(f"Forward primer RC: EXACT match at position {fp_rc_pos}")
            except ValueError:
                max_score, fp_rc_pos, end_pos = self.smith_waterman(sequence, self.forward_primer_rc)
                if max_score >= len(self.forward_primer_rc) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_indel * 2):
                    mismatches = len(self.forward_primer_rc) * 2 - max_score
                    indels = mismatches // 2
                    mismatches = mismatches - (indels * 2)
                    logs.append(f"Forward primer RC: FUZZY match at position {fp_rc_pos}")
                    logs.append(f"Score: {max_score}, Mismatches: {mismatches}, Indels: {indels}")
                else:
                    logs.append("Forward primer RC: NO match found")

        if rp_pos != -1 and fp_rc_pos != -1:
            # 反向结构：well index 紧挨着 reverse primer 前面，plate index 紧挨着 forward primer RC 后面
            reverse_index_search_end = rp_pos
            reverse_index_search_start = max(0, reverse_index_search_end - (self.index_length + self.index_max_indel))
            reverse_index_search_region = sequence[reverse_index_search_start:reverse_index_search_end]
            
            forward_index_search_start = fp_rc_pos + len(self.forward_primer_rc)
            forward_index_search_end = min(len(sequence), forward_index_search_start + self.index_length + self.index_max_indel)
            forward_index_search_region = sequence[forward_index_search_start:forward_index_search_end]
            
            if len(reverse_index_search_region) < min_index_length:
                logs.append("Sequence too short for reverse index extraction")
                return None
                
            if len(forward_index_search_region) < min_index_length:
                logs.append("Sequence too short for forward index extraction")
                return None
            
            logs.append(f"Reverse index search region: {reverse_index_search_region}")
            logs.append(f"Forward index search region: {forward_index_search_region}")
            
            plate_id, plate_seq, plate_score, plate_mismatches, plate_indels = self.find_index_in_region(
                forward_index_search_region, self.plate_indices, is_plate=True)
            
            if plate_id is not None:
                logs.append(f"Plate index: FUZZY match -> {plate_id}")
                logs.append(f"Matched sequence: {plate_seq}")
                logs.append(f"Score: {plate_score}, Mismatches: {plate_mismatches}, Indels: {plate_indels}")
            else:
                for seq in forward_index_search_region.split():
                    if seq in self.plate_dict:
                        plate_id = self.plate_dict[seq]
                        logs.append(f"Plate index: EXACT match -> {plate_id}")
                        break
                
                if plate_id is None:
                    logs.append("Plate index: NO match found")
                    return None
            
            if plate_id is not None:
                well_id, well_seq, well_score, well_mismatches, well_indels = self.find_index_in_region(
                    reverse_index_search_region, self.well_indices, is_plate=False)
                
                if well_id is not None:
                    logs.append(f"Well index: FUZZY match -> {well_id}")
                    logs.append(f"Matched sequence: {well_seq}")
                    logs.append(f"Score: {well_score}, Mismatches: {well_mismatches}, Indels: {well_indels}")
                else:
                    for seq in reverse_index_search_region.split():
                        if seq in self.well_dict:
                            well_id = self.well_dict[seq]
                            logs.append(f"Well index: EXACT match -> {well_id}")
                            break
                    
                    if well_id is None:
                        logs.append("Well index: NO match found")
            
            if plate_id and well_id:
                return {
                    'plate_id': plate_id,
                    'well_id': well_id,
                    'sequence': sequence[reverse_index_search_start:forward_index_search_end],
                    'read_id': seq_record.id,
                    'logs': logs
                }
            
            return None

        return None

    def demultiplex(self, input_file, output_dir, num_processes=None, chunk_size=10000):
        """
        按块读取并处理FASTA文件，支持多进程
        参数:
            input_file: 输入文件路径
            output_dir: 输出目录
            num_processes: 进程数，默认为CPU核心数
            chunk_size: 每次读取的序列数量，默认为10000
        """
        if num_processes is None:
            num_processes = multiprocessing.cpu_count()

        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)

        # 创建进程池
        pool = Pool(processes=num_processes)

        # 用于存储所有结果的字典
        results = defaultdict(list)

        # 分块读取并处理序列
        sequences = []
        for record in SeqIO.parse(input_file, "fasta"):
            sequences.append(record)

            # 当累积的序列达到chunk_size时，进行处理
            if len(sequences) >= chunk_size:
                # 处理当前批次的序列
                chunk_results = pool.map(self.process_sequence, sequences)

                # 将结果添加到总结果中
                for result in chunk_results:
                    if result:
                        # 在主进程中打印日志
                        if 'logs' in result:
                            for log in result['logs']:
                                print(log)
                        
                        # 将结果添加到对应的plate_id列表中
                        results[result['plate_id']].append((result['well_id'], result['sequence'], result['read_id']))

                # 清空序列列表，准备下一批
                sequences = []

        # 处理剩余的序列
        if sequences:
            chunk_results = pool.map(self.process_sequence, sequences)
            for result in chunk_results:
                if result:
                    # 在主进程中打印日志
                    if 'logs' in result:
                        for log in result['logs']:
                            print(log)
                    
                    # 将结果添加到对应的plate_id列表中
                    results[result['plate_id']].append((result['well_id'], result['sequence'], result['read_id']))

        # close pool
        pool.close()
        pool.join()

        # output
        for plate_id, sequences_list in results.items():
            output_file = os.path.join(output_dir, f"{plate_id}.fa")
            with open(output_file, "w") as f:
                for well_id, seq, original_read_id in sequences_list:
                    # 构建所需的 ID 格式 >{read_id}_{plate_id}_{well_id}
                    f.write(f">{original_read_id}_{plate_id}_{well_id}\n{seq}\n")

        return results

def read_plate_index_file(plate_index_file):
    """
    Read plate index file
    Format:
    P1    TCGGTCTTAGACG
    P2    TGTGAAGTTGCCA
    ...
    """
    plate_indices = []
    with open(plate_index_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) == 2:
                    plate_id, index_seq = parts
                    plate_indices.append((plate_id, index_seq))
    return plate_indices

def read_well_index_file(well_index_file):
    """
    Read well index file
    Format:
    001    AGCAATCGCGCAC
    002    AGTGGACCAACAA
    ...
    """
    well_indices = []
    with open(well_index_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) == 2:
                    well_id, index_seq = parts
                    well_indices.append((well_id, index_seq))
    return well_indices

def read_primer_file(primer_file):
    """
    Read primer file
    Format:
    for     TAAACTTCTGGATGTCCAAAAAATCA
    rev     TTTCAACAAATCATAAAGATATTGG
    """
    forward_primer = None
    reverse_primer = None

    with open(primer_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 2:
                    if parts[0].lower() == 'for':
                        forward_primer = parts[1]
                    elif parts[0].lower() == 'rev':
                        reverse_primer = parts[1]

    return forward_primer, reverse_primer

def main():
    parser = argparse.ArgumentParser(description='PCR Product Sequence Demultiplexing Tool')

    parser.add_argument('-p', '--primer', required=True,
                      help='File containing forward and reverse primer sequences')
    parser.add_argument('--plate-index', required=True,
                      help='File containing plate index sequences')
    parser.add_argument('--well-index', required=True,
                      help='File containing well index sequences')
    parser.add_argument('-f', '--fa', required=True,
                      help='Input FASTA format sequence file')
    parser.add_argument('--primer_max_mismatch', type=int, default=3,
                      help='Maximum number of mismatches allowed for primers')
    parser.add_argument('--primer_max_indel', type=int, default=1,
                      help='Maximum number of indels allowed for primers')
    parser.add_argument('--index_max_mismatch', type=int, default=4,
                      help='Maximum number of mismatches allowed for index')
    parser.add_argument('--index_max_indel', type=int, default=1,
                      help='Maximum number of indels allowed for index')
    parser.add_argument('-o', '--output', required=True,
                      help='Output directory path')
    parser.add_argument('-n', '--num_processes', type=int, default=None,
                      help='Number of processes for parallel processing, default uses all available CPU cores')
    parser.add_argument('-c', '--chunk_size', type=int, default=10000,
                      help='Number of sequences to process in each chunk')

    args = parser.parse_args()

    # check file exists
    input_files = {
        'Primer file': args.primer,
        'Plate index file': args.plate_index,
        'Well index file': args.well_index,
        'Input FASTA file': args.fa
    }


    for file_desc, file_path in input_files.items():
        if not os.path.exists(file_path):
            print(f"Error: {file_desc} not found: {file_path}")
            sys.exit(1)

    if os.path.exists(args.output) == True:
        print(f"Error: {args.output} has existed, please remove it first")
        sys.exit(1)

    # Read primer sequences
    forward_primer, reverse_primer = read_primer_file(args.primer)
    if not forward_primer or not reverse_primer:
        print("Error: Failed to read sequences from primer file")
        sys.exit(1)

    # Read plate and well index sequences
    plate_indices = read_plate_index_file(args.plate_index)
    well_indices = read_well_index_file(args.well_index)
    if not plate_indices or not well_indices:
        print("Error: Failed to read sequences from index files")
        sys.exit(1)

    print(f"Processing file: {args.fa}")
    print(f"Forward Primer: {forward_primer}")
    print(f"Reverse Primer: {reverse_primer}")
    print(f"Loaded {len(plate_indices)} plate indices and {len(well_indices)} well indices")

    # Create demultiplexer instance
    demultiplexer = PCRDemultiplexer(forward_primer, reverse_primer,
                                    plate_indices, well_indices,
                                    args.primer_max_mismatch, args.primer_max_indel,
                                    args.index_max_mismatch, args.index_max_indel
                                    )

    # Process sequences
    try:
        results = demultiplexer.demultiplex(args.fa, args.output, args.num_processes)
        print("==========================================")
        print(f"Sequence demultiplexing completed, processed {sum(len(records) for records in results.values())} sequences")
        for combined_id, records in results.items():
            print(f"ID {combined_id}: {len(records)} sequences")
    except Exception as e:
        print(f"Error occurred during processing: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
