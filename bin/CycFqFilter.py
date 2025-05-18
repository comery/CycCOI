import os
import sys
import math
import gzip
import argparse
from matplotlib import pylab
import pylab as plt
import seaborn as sns
import numpy as np
import pandas as pd
import mappy as mp
import matplotlib.patches as patches
from Bio.SeqUtils import gc_fraction

def stat_plot(df, df_filtered=None, ax=None, outpre="Cyclone"):
    # quality plot, boxplot
    qual_output = outpre + ".quality.pdf"
    plt.figure(figsize=(12,5))
    if df_filtered is not None:
        # 创建并排的箱线图，设置不同的颜色
        bp = plt.boxplot([df['quality'], df_filtered['quality']],
                        labels=['Raw Data', 'Filtered Data'],
                        patch_artist=True)

        # 设置箱体颜色
        bp['boxes'][0].set_facecolor('#74a892')
        bp['boxes'][1].set_facecolor('#008585')

        plt.ylabel('Quality Score')
    else:
        bp = plt.boxplot([df['quality']], labels=['Raw Data'],
                        patch_artist=True)
        bp['boxes'][0].set_facecolor('#74a892')
        plt.ylabel('Quality Score')
    plt.title('Sequence Quality Distribution')
    plt.savefig(qual_output, format='pdf')
    plt.close()

    # gc plot, histogram
    gc_output = outpre + ".gc.pdf"
    plt.figure(figsize=(12,5))
    if df_filtered is not None:
        # 计算合适的binwidth
        gc_binwidth = 0.02  # GC含量每2%一个bin
        gc_bins = np.arange(0, 1 + gc_binwidth, gc_binwidth)
        plt.hist(df['gc'], bins=gc_bins, alpha=0.7, color='#74a892', label='Raw Data')
        plt.hist(df_filtered['gc'], bins=gc_bins, alpha=0.7, color='#008585', label='Filtered Data')
        plt.legend()
    else:
        gc_binwidth = 0.02
        gc_bins = np.arange(0, 1 + gc_binwidth, gc_binwidth)
        plt.hist(df['gc'], bins=gc_bins, alpha=0.7, color='#74a892')
    plt.xlabel('GC Content')
    plt.ylabel('Count')
    plt.title('GC Content Distribution')
    plt.savefig(gc_output, format='pdf')
    plt.close()

    # length plot, histogram
    length_output = outpre + ".length.pdf"
    plt.figure(figsize=(12,5))
    if df_filtered is not None:
        # 计算合适的binwidth
        length_binwidth = 100  # 每100bp一个bin
        length_max = max(df['length'].max(), df_filtered['length'].max())
        length_min = min(df['length'].min(), df_filtered['length'].min())
        length_bins = np.arange(length_min, length_max + length_binwidth, length_binwidth)

        plt.hist(df['length'], bins=length_bins, alpha=0.7, color='#74a892', label='Raw Data')
        plt.hist(df_filtered['length'], bins=length_bins, alpha=0.7, color='#008585', label='Filtered Data')
        plt.legend()
    else:
        length_binwidth = 100
        length_bins = np.arange(df['length'].min(), df['length'].max() + length_binwidth, length_binwidth)
        plt.hist(df['length'], bins=length_bins, alpha=0.7, color='#74a892')
    plt.xlabel('Read Length')
    plt.ylabel('Count')
    plt.title('Read Length Distribution')
    plt.savefig(length_output, format='pdf')
    plt.close()

def main(args):
    ids = []
    quality_for_plot = []
    gc_for_plot = []
    len_for_plot = []
    count = 0
    filtered = 0
    if args.plot_only == False:
        outfile = args.outpre + ".clean.fq.gz"
        if os.path.exists(outfile) == True:
            sys.exit(f"{outfile} has existed, please check!")
        fo = gzip.open(args.outpre + ".clean.fq.gz", 'wt')

    # Create new lists for filtered data
    filtered_ids = []
    filtered_quality = []
    filtered_gc = []
    filtered_length = []

    for read in mp.fastx_read(args.fastx, read_comment=False):
        qual = read[0].split("_")[-1]
        qual = float(qual)
        rlen = len(read[1])
        gc = gc_fraction(read[1])
        count += 1

        # 判断是否通过过滤条件，增加GC含量过滤
        passed_filter = (qual >= args.quality_cutoff and
                        rlen >= args.length_cutoff and
                        (args.max_length is None or rlen <= args.max_length) and
                        (args.min_gc is None or gc >= args.min_gc) and
                        (args.max_gc is None or gc <= args.max_gc))

        if count < args.plot_limit:
            ids.append(read[0])
            quality_for_plot.append(qual)
            len_for_plot.append(rlen)
            gc_for_plot.append(gc)

            if passed_filter:
                filtered_ids.append(read[0])
                filtered_quality.append(qual)
                filtered_gc.append(gc)
                filtered_length.append(rlen)

        if args.plot_only == False:
            if passed_filter:
                print(f"@{read[0]}\n{read[1]}\n+\n{read[2]}", file=fo)
            else:
                filtered += 1

    if args.plot_only == False:
        fo.close()
        print(f"filtering done!")
        print(f"   Total Reads: {count}")
        print(f"Filtered Reads: {filtered}")
    # list 2 dict
    data = {
            'id': ids,
            'quality': quality_for_plot,
            'gc': gc_for_plot,
            'length': len_for_plot,
    }

    # dict 2 Pandas DataFrame
    df = pd.DataFrame(data)

    # Create DataFrame for filtered data
    filtered_data = {
        'id': filtered_ids,
        'quality': filtered_quality,
        'gc': filtered_gc,
        'length': filtered_length,
    }
    df_filtered = pd.DataFrame(filtered_data)

    # Print statistics
    print("\nData Statistics:")
    print("Raw Data:")
    print(f"  Average Quality Score: {df['quality'].mean():.2f}")
    print(f"  Average GC Content: {df['gc'].mean():.2f}")
    print(f"  Average Read Length: {df['length'].mean():.2f}")
    print("\nFiltered Data:")
    print(f"  Average Quality Score: {df_filtered['quality'].mean():.2f}")
    print(f"  Average GC Content: {df_filtered['gc'].mean():.2f}")
    print(f"  Average Read Length: {df_filtered['length'].mean():.2f}")

    return df, df_filtered

if __name__ == '__main__':
    if len(sys.argv) < 2:
        usage = f"""
Usage: python3 {sys.argv[0]}  -q 7 -l 1000 -g 0.2 -G 0.8 -o output *.fq.gz # filter and plot
       python3 {sys.argv[0]} --plot_only -q 7 -l 1000 -g 0.2 -G 0.8 -o output *.fq.gz # only plot
       python3 CycFqFilter.py -q 7 -l 1000 -L 5000 -g 0.2 -G 0.8 -o output input.fq.gz  # filter reads with length between 1000-5000 and GC content between 0.2-0.8
        """
        sys.exit(usage)
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument("-o", "--outpre", dest="outpre", metavar="<STR>", type=str,
                            help="prefix for outputs", required=True)
        parser.add_argument("--plot_only", dest="plot_only", action="store_true",
                            help="plot distribution for gc content, quality, read length")
        parser.add_argument("-lim", "--plot_limit", dest="plot_limit", metavar="<INT>", type=int,
                            default=50000, help="item limitation for ploting")
        parser.add_argument("-l", "--len", dest="length_cutoff", metavar="<INT>", type=int,
                            default=1000, help="filtering cutoff for read length")
        parser.add_argument("-L", "--maxlen", dest="max_length", metavar="<INT>", type=int,
                            default=None, help="maximum read length cutoff")
        parser.add_argument("-q", "--qual", dest="quality_cutoff", metavar="<FLOAT>", type=float,
                            default=7.0, help="filtering cutoff for read quality")
        parser.add_argument("-g", "--min_gc", dest="min_gc", metavar="<FLOAT>", type=float,
                            default=None, help="minimum GC content cutoff")
        parser.add_argument("-G", "--max_gc", dest="max_gc", metavar="<FLOAT>", type=float,
                            default=None, help="maximum GC content cutoff")
        parser.add_argument("fastx", type=str,
                            help="input file, fasta or fastq")
        args = parser.parse_args()

        df,df_filtered = main(args)
        print(f"ploting...")
        stat_plot(df, df_filtered=df_filtered, ax=None,  outpre=args.outpre)
        print("ploting done!")
