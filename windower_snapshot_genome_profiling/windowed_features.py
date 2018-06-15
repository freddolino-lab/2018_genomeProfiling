# Copyright (c) 2018 Michael Wolfe University of Michigan. All rights reserved.
#
#
#Developed by: Michael Wolfe
#University of Michigan
#http://freddolino-lab.med.umich.edu/
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal with
#the Software without restriction, including without limitation the rights to
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#the Software, and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:
#
#Redistributions of source code must retain the above copyright notice, this
#list of conditions and the following disclaimers.  Redistributions in binary
#form must reproduce the above copyright notice, this list of conditions and the
#following disclaimers in the documentation and/or other materials provided with
#the distribution.  Neither the names of Michael Wolfe, University of Michigan,
#nor the names of its contributors may be used to endorse or promote products
#derived from this Software without specific prior written permission.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS
#OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
from windowtools import *
import sys
import numpy as np
import argparse
import shlex
import regex

def newSplit(value):
    lex = shlex.shlex(value)
    lex.quotes = '"'
    lex.whitespace_split = True
    lex.commenters = ''
    return list(lex)

def make_comment_dict(gff_entry):
    gff_entry.comment_dict = {}
    keyvalue = gff_entry.comments.split(";")
    for pair in keyvalue:
        key_value= newSplit(pair)
        key = key_value[0]
        if len(key_value) != 2:
            value = ""
        else:
            value = key_value[1]
        gff_entry.comment_dict[key] = value

def parse_genomic_feature_type(args_list):
    if args_list is None:
        raise ValueError("--genomic_feature must be specified if --summary_stat genomic is specified")
    featype = args_list[0]
    arg = args_list[1]
    if featype == "motif":
        motif= regex.compile(arg)
        return lambda x: motif_count_window(x, motif)
    elif featype == "basecontent":
        bases = arg.split(";")
        return lambda x: base_content_window(x, bases=bases)
    else:
        raise ValueError("genomic feature %s is not a valid option"%featype)

def parse_entry_for_name(gff_entry):
    make_comment_dict(gff_entry)
    return (str(gff_entry.comment_dict["Gene"])+"\t"+str(gff_entry.comment_dict["Synonym"]))


def base_content_window(window_signal, bases=["A", "T"]):
    """ 
    Given a sequence, return the AT content of the window

    window - (start, stop) in 1-based coordinates
    genome - FastaEntry object from /home/mbwolfe/src/circ_mapper/fasta
    
    Returns fraction of window that is A or T basepair
    """
    total = 0
    divisor = len(window_signal)/len(bases[0])
    for base in bases:
        total += window_signal.count(base)
    return total/(divisor + 0.0)

def motif_count_window(window_signal, motif):
    """ 
    Given a sequence, return the count of motifs in the sequence

    window - (start, stop) in 1-based coordinates
    motif - regular expression object
    
    Returns number of motifs in the window
    """
    matches = regex.findall(motif, window_signal)
    if matches:
        return len(matches)
    else:
        return 0

def any_over_cutoff(window_signal, cutoff=1):
    return np.sum(window_signal > cutoff) > 0


def parse_bed_into_array(infile, genome_length):
    array = np.zeros(genome_length)
    array[:] = np.nan
    with open(infile, mode="r") as f:
        for line in f:
            linearr = line.rstrip().split("\t")
            start = int(linearr[1])
            end = int(linearr[2])
            array[start:end] = float(linearr[4])
    return array

def read_genome(infile):
#    sys.path.append("/home/mbwolfe/src/circ_mapper")
    import fasta
    genome = fasta.FastaFile()
    with open(infile, mode="r") as f:
        genome.read_whole_file(f)
    chrm=genome.pull_entry(genome.chrm_names()[0])
    return chrm

def pull_chrm_seq(chrm, window, genome_size, circ=True, rc=False):
    start = window[0]-1
    end = window[1]
    return chrm.pull_seq(start, end, circ=circ, rc=rc)


def parse_data_into_array(infile, genome_length):
    infile_name = infile.lower()
    if infile_name.endswith(".gr"):
        return np.loadtxt(infile), complex_window_signal_circular
    elif infile_name.endswith(".bed"):
        return parse_bed_into_array(infile, genome_length), simple_window_signal_circular
    elif infile_name.endswith(".fa") or infile_name.endswith(".fasta") or infile_name.endswith(".fna"):
        return read_genome(infile), pull_chrm_seq
    elif infile.lower().endswith(".npy"):
        data = np.load(infile)
        if len(data.shape) > 1:
            return data, complex_window_signal_circular
        else:
            return data, simple_window_signal_circular
    else:
        raise ValueError("%s file types not supported"%infile)

def sliding_window_main(args):

    data, window_signal = parse_data_into_array(args.input_data, args.genome_length)
    if args.genomic_feature is None:
        if not args.no_finite:
            data[~np.isfinite(data)] = 0
        if not args.no_truncate:
            data[data < 0 ] = 0
        if args.convert_logical:
            if len(data.shape) > 1:
                data[:,1] = data[:,1] > args.convert_logical
            else:
                data = data > args.convert_logical

    summary_stat = args.summary_stat
    direction="+"
    # determine where to write out to
    if args.o:
        outf = open(args.o, "w")
    else:
        outf = sys.stdout
    
    column_names = []
    if args.name == "start":
        column_names.append("Start")
        name_func = lambda x: str(x[0]-1)
    elif args.name == "end":
        column_names.append("End")
        name_func = lambda x: str(x[1])
    elif args.name == "startend":
        column_names.append("Start")
        column_names.append("End")
        name_func = lambda x: str(x[0]-1)+"\t"+str(x[1])
    elif args.name == "center":
        column_names.append("Center")
        name_func = lambda x: str(int(np.mean([x[0], x[1]])-1))
    else:
        raise ValueError("--name %s option not supported. See -h for details"%(args.name))

    # determine the header for the file
    if args.far_upstream > 0:
        for i in xrange(0,args.far_upstream_bins):
            column_names.append("upstream_window_%i"%i)
    for i in xrange(0,args.window_bins):
        column_names.append("window_%i"%i)
    if args.far_downstream > 0:
        for i in xrange(0, args.far_downstream_bins):
            column_names.append("downstream_window_%i"%i)
    outf.write("\t".join(column_names) + "\n")
    if args.plot_dist:
        actual_distros=[]
    for window in sliding_bounds(args.size, args.genome_length, 
                                 args.slide_by):
        outstring = ""
        windows = None
        upstream_windows = None
        downstream_windows = None
        if args.window_bins > 1:
            windows = discretize_window(window, direction, args.window_bins)
        else:
            windows = [window]
        if args.far_upstream > 0:
            upstream_window = add_window(window, args.far_upstream, direction)
            if args.far_upstream_bins > 1:
                upstream_windows = discretize_window(upstream_window, direction, 
                                                     args.far_upstream_bins)
            else:
                upstream_windows = [upstream_window]

        if args.far_downstream > 0:
            downstream_window = add_window(window, args.far_downstream, 
                                           direction,fiveprime=False)
            if args.far_downstream_bins > 1:
                downstream_windows = discretize_window(downstream_window, direction, 
                                                     args.far_downstream_bins)
            else:
                downstream_windows = [downstream_window]
        # output 0 based windows
        outstring += name_func(window)
        if upstream_windows:
            for window in upstream_windows:
                outstring += "\t%0.4e"%(summary_stat(window_signal(data, window, args.genome_length)))
        for window in windows:
            val = summary_stat(window_signal(data, window, args.genome_length))
            if args.plot_dist:
                actual_distros.append(val)
            outstring += "\t%0.4e"%(val)

        if downstream_windows:
            for window in downstream_windows:
                outstring += "\t%0.4e"%(summary_stat(window_signal(data, window, args.genome_length)))
        outf.write(outstring+"\n")
    if args.plot_dist:
        import matplotlib.pyplot as plt
        import seaborn as sns
        if args.o:
            fname = args.o.split(".")[0]
        else:
            fname = "out"

        ax = sns.distplot(np.array(actual_distros), kde=False,color="r")
        plt.savefig("%s_distributions.png"%(fname))


def gff_window_main(args):
    import gfftools

    center_func ={'median':gff_median_center, 'threeprime':gff_threeprime_center,
            'fiveprime':gff_fiveprime_center, 'identity':gff_identity_center}
    center_func = center_func[args.center_metric]
    data, window_signal = parse_data_into_array(args.input_data, args.genome_length)


    if args.genomic_feature is None:
        if not args.no_finite:
            data[~np.isfinite(data)] = 0
        if not args.no_truncate:
            data[data < 0 ] = 0
        if args.convert_logical:
            if len(data.shape) > 1:
                data[:,1] = data[:,1] > args.convert_logical
            else:
                data = data > args.convert_logical

    summary_stat = args.summary_stat
    # determine where to write out to
    if args.o:
        outf = open(args.o, "w")
    else:
        outf = sys.stdout
    
    column_names = []
    if args.name == "fiveprime":
        column_names.append("FiveprimeBase")
        name_func = lambda x: str(x.end-1) if x.direction is "-" else str(x.start-1)
    elif args.name == "threeprime":
        column_names.append("ThreeprimeBase")
        name_func = lambda x: str(x.start-1) if x.direction is "-" else str(x.end-1)
    elif args.name == "startend":
        column_names.append("Start")
        column_names.append("End")
        name_func = lambda x: str(x.start-1)+"\t"+str(x.end)
    elif args.name == "startendstrand":
        column_names.append("Start")
        column_names.append("End")
        column_names.append("Strand")
        name_func = lambda x: str(x.start-1)+"\t"+str(x.end)+"\t"+str(x.direction)
    elif args.name == "comments":
        column_names.append("Gene")
        column_names.append("bnumber")
        name_func = parse_entry_for_name
    elif args.name == "center":
        column_names.append("Center")
        name_func = lambda x: str(int(np.mean([x.start, x.end])-1))
    else:
        raise ValueError("--name %s option not supported. See -h for details"%(args.name))

    # determine the header for the file
    if args.far_upstream > 0:
        for i in xrange(0,args.far_upstream_bins):
            column_names.append("upstream_window_%i"%i)
    for i in xrange(0,args.window_bins):
        column_names.append("window_%i"%i)
    if args.far_downstream > 0:
        for i in xrange(0, args.far_downstream_bins):
            column_names.append("downstream_window_%i"%i)
    outf.write("\t".join(column_names) + "\n")
    if args.distributions:
        actual_distros=[]
        random_distros=[]
    gffs = gfftools.GffData()
    gffs.parse_gff_file(args.gff_windows)
    for entry in gffs:
        outstring = ""
        windows = None
        upstream_windows = None
        downstream_windows = None
        outstring += name_func(entry)
        window = get_gff_window(entry, args.upstream, args.downstream, center_func)
        if args.distributions:
            for x in xrange(args.distributions):
                random_distros.append(summary_stat(window_signal(data, get_random_window(window[1]-window[0]+1, args.genome_length), args.genome_length)))
        if args.window_bins > 1:
            windows = discretize_window(window, entry.direction, args.window_bins)
        else:
            windows = [window]
        if args.far_upstream > 0:
            upstream_window = add_window(window, args.far_upstream, entry.direction)
            if args.far_upstream_bins > 1:
                upstream_windows = discretize_window(upstream_window, entry.direction, 
                                                     args.far_upstream_bins)
            else:
                upstream_windows = [upstream_window]

        if args.far_downstream > 0:
            downstream_window = add_window(window, args.far_downstream, 
                                           entry.direction,fiveprime=False)
            if args.far_downstream_bins > 1:
                downstream_windows = discretize_window(downstream_window, entry.direction, 
                                                     args.far_downstream_bins)
            else:
                downstream_windows = [downstream_window]
        if upstream_windows:
            for window in upstream_windows:
                outstring += "\t%0.4e"%(summary_stat(window_signal(data, window, args.genome_length)))
        for window in windows:
            val = summary_stat(window_signal(data, window, args.genome_length))
            if args.distributions:
                actual_distros.append(val)
            outstring += "\t%0.4e"%(val)

        if downstream_windows:
            for window in downstream_windows:
                outstring += "\t%0.4e"%(summary_stat(window_signal(data, window, args.genome_length)))
        outf.write(outstring+"\n")
    if args.distributions:
        import matplotlib.pyplot as plt
        import seaborn as sns
        if args.o:
            fname = args.o.split(".")[0]
        else:
            fname = "out"

        ax = sns.distplot(np.array(actual_distros), color="r")
        ax = sns.distplot(np.array(random_distros), color="b")
        plt.savefig("%s_distributions.png"%(fname))



if __name__ == "__main__":
    # windowing control
    parent_parser = argparse.ArgumentParser(description="Program to calculate a summary statistic over specified windows in a circular genome. The options below will work on each command. To see additional help for each command, specify 'command -h'.\
            Both .gr files and .npy files are assumed to be in 0-based coordinates. All output is also in 0-based coordinates\
            .gff files are assumed to be in 1-based coordinates")
    subparsers = parent_parser.add_subparsers(title="window type commands. User must specify one of these.", description="Each command specifies a different way of defining windows", dest="command")
    parent_parser.add_argument('--window_bins', action='store', type=int, default=1,
                        help="how many bins to break the window into")
    parent_parser.add_argument('--far_upstream', action='store', type=int, default=0,
                            help='number of bases upstream of a region, default=0')
    parent_parser.add_argument('--far_upstream_bins', action='store', type=int, default=1,
                            help='number of bins to split region into, default=1')
    parent_parser.add_argument('--far_downstream', action='store', type=int, default=0,
                            help='number of bases downstream of a region, default=0')
    parent_parser.add_argument('--far_downstream_bins', action='store', type=int, default=1,
                            help='number of bins to split region into, default=1')
    parent_parser.add_argument('--genome_length', action="store", type=int,
                        default=4639676,
                        help="length of the genome for coverage file, default=\
                                4639676")
    parent_parser.add_argument('--summary_stat', action="store", type=str,
                               default="mean",
                               help="summary stat to calculate for each window.\
                                     [mean = mean of the signal in the window,\
                                      median = median of the signal in the window,\
                                      cutoff = if a number is specified then return 1\
                                      if any datapoint is over the cutoff,\
                                      density = determine data density for a window\
                                      i.e. number of data points in the window\
                                      genomic = choosing a genomic feature, must specify --genomic_feature flag]\
                                      default= mean.")
    parent_parser.add_argument("--genomic_feature", nargs=2, type=str, help="When --summary_stat genomic is specified, this flag MUST be specified\
                                                                             Choose the type of genomic feature desired:\
                                                                             [motif regex gives count of motifs in the window\
                                                                              basecontent base1;base2;...etc. \
                                                                              gives fraction of window that is/are that/those base(s). \
                                                                              Also works with dinucleotide, tri... etc. Specifiying\
                                                                              a mix of tri/di/mono bases will give meaningless answers]")

    parent_parser.add_argument('-o', action='store', type=str,
                              help='output file to put the output into (default stdout)')

    parent_parser.add_argument('--no_truncate', action="store_true",help="don't truncate\
                        values below zero to zero.")
    parent_parser.add_argument('--no_finite', action="store_true",help="don't truncate\
                        non-finite values to zero.")
    parent_parser.add_argument('--convert_logical', type=float, default=None,help="convert the data array to a logical array that is true if the data point is above a cutoff value")

    sliding_parser = subparsers.add_parser('sliding', help="do a sliding window over the genome")
    sliding_parser.add_argument("input_data", action='store', type=str,
                               help="input data. Accepted types include .gr,\
                                     .npy")
    sliding_parser.add_argument("size", type=int,help="size of the sliding window")
    sliding_parser.add_argument("slide_by", type=int,help="number of basepairs to slide by")
    sliding_parser.add_argument("--plot_dist", action="store_true", help="plot the distribution of the values in the main window")
    sliding_parser.add_argument('--name', action='store', type=str, default='startend',
                                help='[startend = put the start and end of the region,\
                                       start = put only the starting bp,\
                                       end = put only the ending bp,\
                                       center = put the center of the region median(start,end)]\
                                       default=startend')
     
    gff_parser = subparsers.add_parser('gff_window', help="define windows from a gff")
    gff_parser.add_argument("input_data", action='store', type=str,
                               help="input data. Accepted types include .gr,\
                                     .npy")
    gff_parser.add_argument('gff_windows', action='store', type=str,
                              help='gff file containing windows to scan')
    gff_parser.add_argument('--center_metric', action='store', type=str, default='identity',
                            help='How to determine the center of the gff location?\
                            [median = center is defined as the median,\
                             threeprime = center is defined as the three prime end of the region,\
                             fiveprime = center is defined as the five prime end of the region,\
                             identity = entire region is taken as a window + upstream and downstream bp]\
                             Default = identity')
    parent_parser.add_argument('--upstream', action='store', type=int, default=0,
                            help='number of bases upstream of a region, default=0')
    parent_parser.add_argument('--downstream', action='store', type=int, default=0,
                            help='number of bases downstream of a region, default=0')
    gff_parser.add_argument('--name', action='store', type=str, default='startend',
                            help='[comments= parse comments for gene name/bnumber,\
                                   fiveprime= just put the fiveprime bp,\
                                   threeprime = just put the threeprime bp,\
                                   startend = put the start and end of the region\
                                   startendstrand = put the start, end and strand of the region\
                                   center = put the center of the region median(start,end)]\
                                   default=startend')

    gff_parser.add_argument('--distributions', action="store", type=int,
                         help="plot the distribution of the data set in the ranges compared\
                               to randomly sample bins of the same length")
    args = parent_parser.parse_args()

    gffs = gfftools.GffData()
    summary_stat = {'mean': np.mean, 'median': np.median, 'density': lambda x: len(x)}
    if args.summary_stat in summary_stat.keys():
        args.summary_stat = summary_stat[args.summary_stat]
    elif args.summary_stat == "genomic":
        args.summary_stat = parse_genomic_feature_type(args.genomic_feature) 
    else:
        try:
            cutoff = float(args.summary_stat)
            args.summary_stat = lambda x : any_over_cutoff(x, cutoff)
        except:
            print "Error interpreting number for cutoff summary stat %s"%args.summary_stat
            raise
    if args.command == "sliding": 
        sliding_window_main(args)
    elif args.command == "gff_window":
        gff_window_main(args)
    else:
        raise ValueError("unsupported command")
