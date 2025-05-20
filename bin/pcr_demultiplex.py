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
                 primer_max_mismatch=3, primer_max_gap=1,
                 index_max_mismatch=4, index_max_gap=1):
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
        self.primer_max_gap = primer_max_gap
        self.index_max_mismatch = index_max_mismatch
        self.index_max_gap = index_max_gap

        # 预计算所有可能需要的反向互补序列
        self.index_rc_cache = {}
        for _, index_seq in plate_indices + well_indices:
            self.index_rc_cache[index_seq] = self.reverse_complement(index_seq)

        # 创建index序列的精确匹配字典
        self.plate_dict = {}
        self.well_dict = {}

        # plate index dict
        for plate_id, index_seq in plate_indices:
            self.plate_dict[index_seq] = plate_id
            self.plate_dict[self.reverse_complement(index_seq)] = plate_id

        # well index dict
        for well_id, index_seq in well_indices:
            self.well_dict[index_seq] = well_id
            self.well_dict[self.reverse_complement(index_seq)] = well_id

    def reverse_complement(self, sequence):
        """
        Calculate the reverse complement sequence
        """
        return str(Seq(sequence).reverse_complement())  # 修改这里

    def needleman_wunsch(self, seq1, seq2, match_score=2, mismatch_score=-1, gap_score=-2):
        """
        Optimized Needleman-Wunsch algorithm that terminates early for impossible matches
        """
        m, n = len(seq1), len(seq2)
        max_allowed_penalty = min(self.primer_max_mismatch, self.index_max_mismatch)
        min_required_score = len(seq2) - max_allowed_penalty

        # 如果序列长度差异太大，直接返回一个很低的分数
        if abs(m - n) > max_allowed_penalty:
            return float('-inf')

        score_matrix = np.zeros((m+1, n+1))

        for i in range(m+1):
            score_matrix[i][0] = i * gap_score
        for j in range(n+1):
            score_matrix[0][j] = j * gap_score

        # 提前终止检查
        for i in range(1, m+1):
            if score_matrix[i][0] < min_required_score - (n * match_score):
                return float('-inf')

            for j in range(1, n+1):
                match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
                delete = score_matrix[i-1][j] + gap_score
                insert = score_matrix[i][j-1] + gap_score
                score_matrix[i][j] = max(match, delete, insert)

        return score_matrix[m][n]

    def smith_waterman(self, seq1, seq2, match_score=2, mismatch_score=-1, gap_score=-2):
        """
        Perform local sequence alignment using Smith-Waterman algorithm
        """
        m, n = len(seq1), len(seq2)
        score_matrix = np.zeros((m+1, n+1))

        # 初始化第一行和第一列为0（与Needleman-Wunsch不同）
        for i in range(m+1):
            score_matrix[i][0] = 0
        for j in range(n+1):
            score_matrix[0][j] = 0

        # 记录最高分及其位置
        max_score = 0
        max_pos = (0, 0)

        # 填充得分矩阵
        for i in range(1, m+1):
            for j in range(1, n+1):
                match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
                delete = score_matrix[i-1][j] + gap_score
                insert = score_matrix[i][j-1] + gap_score
                score_matrix[i][j] = max(0, match, delete, insert)  # 与全局比对的关键区别

                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)

        return max_score, max_pos[0]

    def find_primer_position(self, sequence, primer):
        """
        Find the best matching position for primer in sequence
        First try exact match, if failed then use Smith-Waterman algorithm for local alignment
        """
        # 首先尝试精确匹配
        try:
            pos = sequence.index(primer)
            return pos
        except ValueError:
            pass

        # 如果精确匹配失败，使用局部比对算法
        max_score, end_pos = self.smith_waterman(sequence, primer)

        # 计算起始位置（end_pos是匹配结束的位置）
        start_pos = end_pos - len(primer)

        # 检查得分是否满足阈值要求
        min_required_score = len(primer) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2)
        if max_score < min_required_score:
            return -1

        return start_pos

    def match_plate_index(self, sequence):
        """
        Match plate index sequence
        First try exact match, if failed then perform mismatch-tolerant matching
        """
        # 首先尝试精确匹配
        if sequence in self.plate_dict:
            return self.plate_dict[sequence]

        # 如果精确匹配失败，进行容错匹配
        best_match = None
        best_score = float('-inf')

        for plate_id, index_seq in self.plate_indices:
            score = self.needleman_wunsch(sequence, index_seq)
            score_rc = self.needleman_wunsch(sequence, self.index_rc_cache[index_seq])

            max_score = max(score, score_rc)
            if max_score >= len(index_seq) - self.index_max_mismatch - self.index_max_gap:
                if max_score > best_score:
                    best_score = max_score
                    best_match = plate_id

        return best_match

    def match_well_index(self, sequence):
        """
        Match well index sequence
        First try exact match, if failed then perform mismatch-tolerant matching
        """
        # 首先尝试精确匹配
        if sequence in self.well_dict:
            return self.well_dict[sequence]

        # 如果精确匹配失败，进行容错匹配
        best_match = None
        best_score = float('-inf')

        for well_id, index_seq in self.well_indices:
            score = self.needleman_wunsch(sequence, index_seq)
            score_rc = self.needleman_wunsch(sequence, self.index_rc_cache[index_seq])

            max_score = max(score, score_rc)
            if max_score >= len(index_seq) - self.index_max_mismatch - self.index_max_gap:
                if max_score > best_score:
                    best_score = max_score
                    best_match = well_id

        return best_match

    def process_sequence(self, seq_record):
        sequence = str(seq_record.seq)

        # 使用列表收集日志信息，而不是直接打印
        logs = []
        logs.append("==========================================")
        logs.append(f"Processing sequence: {seq_record.id}")

        # Initialize plate_id and well_id as None at the start
        plate_id = None
        well_id = None

        # 计算最小所需的index长度（考虑允许的gap）
        min_index_length = self.index_length - self.index_max_gap

        # Forward-reverse structure
        fp_pos = -1
        rp_rc_pos = -1

        # Try exact match for forward primer
        try:
            fp_pos = sequence.index(self.forward_primer)
            logs.append(f"Forward primer: EXACT match at position {fp_pos}")
        except ValueError:
            max_score, end_pos = self.smith_waterman(sequence, self.forward_primer)
            if max_score >= len(self.forward_primer) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2):
                fp_pos = end_pos - len(self.forward_primer)
                mismatches = len(self.forward_primer) * 2 - max_score
                gaps = mismatches // 2
                mismatches = mismatches - (gaps * 2)
                logs.append(f"Forward primer: FUZZY match at position {fp_pos}")
                logs.append(f"Score: {max_score}, Mismatches: {mismatches}, Gaps: {gaps}")
            else:
                logs.append("Forward primer: NO match found")

        # if find forward primer, try to find reverse primer
        if fp_pos != -1:
            try:
                rp_rc_pos = sequence.index(self.reverse_primer_rc)
                logs.append(f"Reverse primer RC EXACT match found at position {rp_rc_pos}")
            except ValueError:
                max_score, end_pos = self.smith_waterman(sequence, self.reverse_primer_rc)
                if max_score >= len(self.reverse_primer_rc) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2):
                    rp_rc_pos = end_pos - len(self.reverse_primer_rc)
                    logs.append(f"Reverse primer RC fuzzy match found at position {rp_rc_pos}, score: {max_score}")
                else:
                    logs.append("Reverse primer RC match failed")

        # if find forward primer and reverse primer, try to find index
        if fp_pos != -1 and rp_rc_pos != -1:
            # 检查forward index区域是否有足够的长度
            if fp_pos < min_index_length:
                logs.append("Sequence too short for forward index extraction")
                return None

            # 检查reverse index区域是否有足够的长度
            if len(sequence) - (rp_rc_pos + len(self.reverse_primer_rc)) < min_index_length:
                logs.append("Sequence too short for reverse index extraction")
                return None

            forward_index_seq = sequence[fp_pos-self.index_length:fp_pos]
            reverse_index_rc_seq = sequence[rp_rc_pos+len(self.reverse_primer_rc):
                                         rp_rc_pos+len(self.reverse_primer_rc)+self.index_length]

            logs.append(f"Forward index sequence: {forward_index_seq}")
            logs.append(f"Reverse index RC sequence: {reverse_index_rc_seq}")

            # Match plate index
            plate_id = None  # Reset plate_id before matching
            if forward_index_seq in self.plate_dict:
                plate_id = self.plate_dict[forward_index_seq]
                logs.append(f"Plate index: EXACT match -> {plate_id}")
            else:
                best_score = float('-inf')
                best_plate_id = None
                best_mismatches = 0
                best_gaps = 0

                for pid, index_seq in self.plate_indices:
                    score = self.needleman_wunsch(forward_index_seq, index_seq)
                    score_rc = self.needleman_wunsch(forward_index_seq, self.index_rc_cache[index_seq])
                    max_score = max(score, score_rc)

                    if max_score >= len(index_seq) - self.index_max_mismatch - self.index_max_gap:
                        mismatches = len(index_seq) * 2 - max_score
                        gaps = mismatches // 2
                        mismatches = mismatches - (gaps * 2)

                        if max_score > best_score:
                            best_score = max_score
                            best_plate_id = pid
                            best_mismatches = mismatches
                            best_gaps = gaps

                if best_plate_id is not None:
                    plate_id = best_plate_id
                    logs.append(f"Plate index: FUZZY match -> {plate_id}")
                    logs.append(f"Score: {best_score}, Mismatches: {best_mismatches}, Gaps: {best_gaps}")
                else:
                    logs.append("Plate index: NO match found")
                    return None  # 如果plate index没有匹配到，直接返回None

            # 只有在plate_id匹配成功的情况下才进行well index匹配
            if plate_id is not None:
                # Match well index
                if reverse_index_rc_seq in self.well_dict:
                    well_id = self.well_dict[reverse_index_rc_seq]
                    logs.append(f"Well index: EXACT match -> {well_id}")
                else:
                    best_score = float('-inf')
                    best_well_id = None
                    best_mismatches = 0
                    best_gaps = 0

                    for wid, index_seq in self.well_indices:
                        score = self.needleman_wunsch(reverse_index_rc_seq, index_seq)
                        score_rc = self.needleman_wunsch(reverse_index_rc_seq, self.index_rc_cache[index_seq])
                        max_score = max(score, score_rc)

                        if max_score >= len(index_seq) - self.index_max_mismatch - self.index_max_gap:
                            mismatches = len(index_seq) * 2 - max_score
                            gaps = mismatches // 2
                            mismatches = mismatches - (gaps * 2)

                            if max_score > best_score:
                                best_score = max_score
                                best_well_id = wid
                                best_mismatches = mismatches
                                best_gaps = gaps

                    if best_well_id is not None:
                        well_id = best_well_id
                        logs.append(f"Well index: FUZZY match -> {well_id}")
                        logs.append(f"Score: {best_score}, Mismatches: {best_mismatches}, Gaps: {best_gaps}")
                    else:
                        logs.append("Well index: NO match found")

                if plate_id and well_id:
                    # 对于forward-reverse结构，返回字典格式的结果
                    return {
                        'plate_id': plate_id,
                        'well_id': well_id,
                        'sequence': sequence[fp_pos-self.index_length:rp_rc_pos+len(self.reverse_primer_rc)+self.index_length],
                        'read_id': seq_record.id,
                        'logs': logs  # 添加日志信息
                    }

            # 如果第一种情况匹配失败，才继续尝试第二种情况
            return None

        # 检查第二种情况：reverse-forward结构
        # 同样分别检查正向和反向引物
        rp_pos = -1
        fp_rc_pos = -1

        # 先尝试反向引物的精确匹配
        try:
            rp_pos = sequence.index(self.reverse_primer)
            logs.append(f"Reverse primer: EXACT match at position {rp_pos}")
        except ValueError:
            max_score, end_pos = self.smith_waterman(sequence, self.reverse_primer)
            if max_score >= len(self.reverse_primer) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2):
                rp_pos = end_pos - len(self.reverse_primer)
                mismatches = len(self.reverse_primer) * 2 - max_score
                gaps = mismatches // 2
                mismatches = mismatches - (gaps * 2)
                logs.append(f"Reverse primer: FUZZY match at position {rp_pos}")
                logs.append(f"Score: {max_score}, Mismatches: {mismatches}, Gaps: {gaps}")
            else:
                logs.append("Reverse primer: NO match found")

        # 如果找到反向引物，再尝试正向引物
        if rp_pos != -1:
            try:
                fp_rc_pos = sequence.index(self.forward_primer_rc)
                logs.append(f"Forward primer RC: EXACT match at position {fp_rc_pos}")
            except ValueError:
                max_score, end_pos = self.smith_waterman(sequence, self.forward_primer_rc)
                if max_score >= len(self.forward_primer_rc) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2):
                    fp_rc_pos = end_pos - len(self.forward_primer_rc)
                    mismatches = len(self.forward_primer_rc) * 2 - max_score
                    gaps = mismatches // 2
                    mismatches = mismatches - (gaps * 2)
                    logs.append(f"Forward primer RC: FUZZY match at position {fp_rc_pos}")
                    logs.append(f"Score: {max_score}, Mismatches: {mismatches}, Gaps: {gaps}")
                else:
                    logs.append("Forward primer RC: NO match found")

        # if find reverse primer and forward primer, try to find index
        if rp_pos != -1 and fp_rc_pos != -1:
            # 检查reverse index区域是否有足够的长度
            if rp_pos < min_index_length:
                logs.append("Sequence too short for reverse index extraction")
                return None

            # 检查forward index区域是否有足够的长度
            if len(sequence) - (fp_rc_pos + len(self.forward_primer_rc)) < min_index_length:
                logs.append("Sequence too short for forward index extraction")
                return None

            reverse_index_seq = sequence[rp_pos-self.index_length:rp_pos]
            forward_index_rc_seq = sequence[fp_rc_pos+len(self.forward_primer_rc):
                                         fp_rc_pos+len(self.forward_primer_rc)+self.index_length]

            logs.append(f"Reverse index sequence: {reverse_index_seq}")
            logs.append(f"Forward index RC sequence: {forward_index_rc_seq}")

            # Match plate index
            if forward_index_rc_seq in self.plate_dict:
                plate_id = self.plate_dict[forward_index_rc_seq]
                logs.append(f"Plate index: EXACT match -> {plate_id}")
            else:
                best_score = float('-inf')
                best_plate_id = None
                best_mismatches = 0
                best_gaps = 0

                for pid, index_seq in self.plate_indices:
                    score = self.needleman_wunsch(forward_index_rc_seq, index_seq)
                    score_rc = self.needleman_wunsch(forward_index_rc_seq, self.index_rc_cache[index_seq])
                    max_score = max(score, score_rc)

                    if max_score >= len(index_seq) - self.index_max_mismatch - self.index_max_gap:
                        mismatches = len(index_seq) * 2 - max_score
                        gaps = mismatches // 2
                        mismatches = mismatches - (gaps * 2)

                        if max_score > best_score:
                            best_score = max_score
                            best_plate_id = pid
                            best_mismatches = mismatches
                            best_gaps = gaps

                if best_plate_id is not None:
                    plate_id = best_plate_id
                    logs.append(f"Plate index: FUZZY match -> {plate_id}")
                    logs.append(f"Score: {best_score}, Mismatches: {best_mismatches}, Gaps: {best_gaps}")
                else:
                    logs.append("Plate index: NO match found")
                    return None  # 如果plate index没有匹配到，直接返回None

            # 只有在plate_id匹配成功的情况下才进行well index匹配
            if plate_id is not None:
                # Match well index
                if reverse_index_seq in self.well_dict:
                    well_id = self.well_dict[reverse_index_seq]
                    logs.append(f"Well index: EXACT match -> {well_id}")
                else:
                    best_score = float('-inf')
                    best_well_id = None
                    best_mismatches = 0
                    best_gaps = 0

                    for wid, index_seq in self.well_indices:
                        score = self.needleman_wunsch(reverse_index_seq, index_seq)
                        score_rc = self.needleman_wunsch(reverse_index_seq, self.index_rc_cache[index_seq])
                        max_score = max(score, score_rc)

                        if max_score >= len(index_seq) - self.index_max_mismatch - self.index_max_gap:
                            mismatches = len(index_seq) * 2 - max_score
                            gaps = mismatches // 2
                            mismatches = mismatches - (gaps * 2)

                            if max_score > best_score:
                                best_score = max_score
                                best_well_id = wid
                                best_mismatches = mismatches
                                best_gaps = gaps

                    if best_well_id is not None:
                        well_id = best_well_id
                        logs.append(f"Well index: FUZZY match -> {well_id}")
                        logs.append(f"Score: {best_score}, Mismatches: {best_mismatches}, Gaps: {best_gaps}")
                    else:
                        logs.append("Well index: NO match found")

                if plate_id and well_id:
                    # 返回完整的信息，包括plate_id, well_id和序列
                    return {
                        'plate_id': plate_id,
                        'well_id': well_id,
                        'sequence': sequence[rp_pos-self.index_length:fp_rc_pos+len(self.forward_primer_rc)+self.index_length],
                        'read_id': seq_record.id,
                        'logs': logs  # 添加日志信息到第二种情况的返回值
                    }

            # 如果第二种情况匹配失败，返回None
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

        # 关闭进程池
        pool.close()
        pool.join()

        # 将结果写入文件
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
    parser.add_argument('-pmm', '--primer_max_mismatch', type=int, default=3,
                      help='Maximum number of mismatches allowed for primers')
    parser.add_argument('-pmg', '--primer_max_gap', type=int, default=1,
                      help='Maximum number of gaps allowed for primers')
    parser.add_argument('-imm', '--index_max_mismatch', type=int, default=4,
                      help='Maximum number of mismatches allowed for index')
    parser.add_argument('-img', '--index_max_gap', type=int, default=1,
                      help='Maximum number of gaps allowed for index')
    parser.add_argument('-o', '--output', required=True,
                      help='Output directory path')
    parser.add_argument('-n', '--num_processes', type=int, default=None,
                      help='Number of processes for parallel processing, default uses all available CPU cores')
    parser.add_argument('-c', '--chunk_size', type=int, default=10000,
                      help='Number of sequences to process in each chunk')

    args = parser.parse_args()

    # 检查所有输入文件是否存在
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
                                    args.primer_max_mismatch, args.primer_max_gap,
                                    args.index_max_mismatch, args.index_max_gap
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
