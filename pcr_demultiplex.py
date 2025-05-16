from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import defaultdict
import argparse
import logging
import sys
import os
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
        
        # 板子index字典
        for plate_id, index_seq in plate_indices:
            self.plate_dict[index_seq] = plate_id
            self.plate_dict[self.reverse_complement(index_seq)] = plate_id
            
        # 孔index字典
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
        """
        Process single sequence, supporting both orientations of the structure
        """
        sequence = str(seq_record.seq)
        
        # 检查第一种情况：forward-reverse结构
        # 分别检查正向和反向引物
        fp_pos = -1
        rp_rc_pos = -1
        
        # 先尝试正向引物的精确匹配
        try:
            fp_pos = sequence.index(self.forward_primer)
        except ValueError:
            # 如果精确匹配失败，使用Smith-Waterman算法
            max_score, end_pos = self.smith_waterman(sequence, self.forward_primer)
            if max_score >= len(self.forward_primer) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2):
                fp_pos = end_pos - len(self.forward_primer)
        
        # 如果找到正向引物，再尝试反向引物
        if fp_pos != -1:
            # 尝试反向引物的精确匹配
            try:
                rp_rc_pos = sequence.index(self.reverse_primer_rc)
            except ValueError:
                # 如果精确匹配失败，使用Smith-Waterman算法
                max_score, end_pos = self.smith_waterman(sequence, self.reverse_primer_rc)
                if max_score >= len(self.reverse_primer_rc) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2):
                    rp_rc_pos = end_pos - len(self.reverse_primer_rc)

        if fp_pos != -1 and rp_rc_pos != -1:
            # 提取index序列
            forward_index_seq = sequence[fp_pos-self.index_length:fp_pos]
            reverse_index_rc_seq = sequence[rp_rc_pos+len(self.reverse_primer_rc):
                                         rp_rc_pos+len(self.reverse_primer_rc)+self.index_length]
            
            plate_id = self.match_plate_index(forward_index_seq)
            well_id = self.match_well_index(reverse_index_rc_seq)
            
            if plate_id and well_id:
                return (plate_id, well_id, 
                       sequence[fp_pos+len(self.forward_primer):rp_rc_pos])
                       
        # 检查第二种情况：reverse-forward结构
        # 同样分别检查正向和反向引物
        rp_pos = -1
        fp_rc_pos = -1
        
        # 先尝试反向引物的精确匹配
        try:
            rp_pos = sequence.index(self.reverse_primer)
        except ValueError:
            # 如果精确匹配失败，使用Smith-Waterman算法
            max_score, end_pos = self.smith_waterman(sequence, self.reverse_primer)
            if max_score >= len(self.reverse_primer) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2):
                rp_pos = end_pos - len(self.reverse_primer)

        # 如果找到反向引物，再尝试正向引物
        if rp_pos != -1:
            # 尝试正向引物的精确匹配
            try:
                fp_rc_pos = sequence.index(self.forward_primer_rc)
            except ValueError:
                # 如果精确匹配失败，使用Smith-Waterman算法
                max_score, end_pos = self.smith_waterman(sequence, self.forward_primer_rc)
                if max_score >= len(self.forward_primer_rc) * 2 - (self.primer_max_mismatch * 1 + self.primer_max_gap * 2):
                    fp_rc_pos = end_pos - len(self.forward_primer_rc)

        if rp_pos != -1 and fp_rc_pos != -1:
            # 提取index序列
            reverse_index_seq = sequence[rp_pos-self.index_length:rp_pos]
            forward_index_rc_seq = sequence[fp_rc_pos+len(self.forward_primer_rc):
                                         fp_rc_pos+len(self.forward_primer_rc)+self.index_length]

            # 修改这里：使用新的匹配方法
            plate_id = self.match_plate_index(forward_index_rc_seq)
            well_id = self.match_well_index(reverse_index_seq)

            if plate_id and well_id:
                return (plate_id, well_id,
                       sequence[rp_pos+len(self.reverse_primer):fp_rc_pos])

        return None

    def demultiplex(self, input_file, output_dir, num_processes=None):
        """
        Demultiplex input FASTA file with multi-processing support
        Organize output by plate ID and well ID
        """
        if num_processes is None:
            num_processes = multiprocessing.cpu_count()
            
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
            
        # 读取所有序列记录
        records = list(SeqIO.parse(input_file, "fasta"))
        
        # 创建进程池
        with Pool(processes=num_processes) as pool:
            results = pool.map(self.process_sequence, records)
            
        # 按板子编号组织结果
        plate_well_records = defaultdict(lambda: defaultdict(list))
        for record, result in zip(records, results):
            if result:
                plate_id, well_id, sequence = result
                new_record = SeqRecord(
                    Seq(sequence),
                    id=f"{record.id}_{plate_id}_{well_id}",
                    description=""
                )
                plate_well_records[plate_id][well_id].append(new_record)
        
        # 写入文件
        for plate_id, well_records in plate_well_records.items():
            # 在输出目录下创建孔的文件
            for well_id, records in well_records.items():
                output_file = os.path.join(output_dir, f"{plate_id}_{well_id}.fasta")
                SeqIO.write(records, output_file, "fasta")
            
        return plate_well_records

def setup_logger(log_file):
    """
    Set up logging configuration
    """
    logger = logging.getLogger('PCRDemultiplexer')
    logger.setLevel(logging.INFO)
    
    # 文件处理器
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # 控制台处理器
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # 格式化器
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

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
    parser.add_argument('-l', '--log', default='demultiplex.log',
                      help='Log file path')
    parser.add_argument('-n', '--num_processes', type=int, default=None,
                      help='Number of processes for parallel processing, default uses all available CPU cores')
    
    args = parser.parse_args()
    
    # Set up logger
    logger = setup_logger(args.log)
    
    # Read primer sequences
    forward_primer, reverse_primer = read_primer_file(args.primer)
    if not forward_primer or not reverse_primer:
        logger.error("Failed to read sequences from primer file")
        sys.exit(1)
    
    # Read plate and well index sequences
    plate_indices = read_plate_index_file(args.plate_index)
    well_indices = read_well_index_file(args.well_index)
    if not plate_indices or not well_indices:
        logger.error("Failed to read sequences from index files")
        sys.exit(1)
    
    logger.info(f"Processing file: {args.fa}")
    logger.info(f"Forward Primer: {forward_primer}")
    logger.info(f"Reverse Primer: {reverse_primer}")
    logger.info(f"Loaded {len(plate_indices)} plate indices and {len(well_indices)} well indices")
    
    # Create demultiplexer instance
    demultiplexer = PCRDemultiplexer(
        forward_primer, reverse_primer, plate_indices, well_indices,
        args.primer_max_mismatch, args.primer_max_gap,
        args.index_max_mismatch, args.index_max_gap
    )
    
    # Process sequences
    try:
        results = demultiplexer.demultiplex(args.fa, args.output, args.num_processes)
        logger.info(f"Sequence demultiplexing completed, processed {sum(len(records) for records in results.values())} sequences")
        for combined_id, records in results.items():
            logger.info(f"ID {combined_id}: {len(records)} sequences")
    except Exception as e:
        logger.error(f"Error occurred during processing: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
