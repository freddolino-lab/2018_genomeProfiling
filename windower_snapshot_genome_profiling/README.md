# General Purpose Windowing script #

The intention of this piece of software is a general purpose windowing
script to be used with genomic data across a single circular genome.
The general design idea is to extract a statistic of choice from user
specified windows over input files containing values for each bp
across the genome.

There are two major modes of specifying windows

1. Sliding windows where the size and amount to slide by is specified
2. Specified windows, where the user supplies the windows of interest
   through an input .gff file

Below are a couple examples designed to show off the usefulness of
this script. For each of the examples we will use a 21 bp genome with
the following .gr file as input:

```
0 0
1 1
2 2
3 3
4 4
5 5
6 6
7 7
8 8
9 9
10 10
11 11
12 12
13 13
14 14
15 15
16 16
17 17
18 18
19 19
20 20
```

Where the left column is the bp in 0-based coordinates and the right
column are the values associated with each bp

To do a simple 3bp sliding window and calculate the average value in each
window we would do the following:

```
python windowed_features.py --genome_length 21 sliding tests/test.gr 3 1
Start   End window_0
0   3   1.0000e+00
1   4   2.0000e+00
2   5   3.0000e+00
3   6   4.0000e+00
4   7   5.0000e+00
5   8   6.0000e+00
6   9   7.0000e+00
7   10  8.0000e+00
8   11  9.0000e+00
9   12  1.0000e+01
10  13  1.1000e+01
11  14  1.2000e+01
12  15  1.3000e+01
13  16  1.4000e+01
14  17  1.5000e+01
15  18  1.6000e+01
16  19  1.7000e+01
17  20  1.8000e+01
18  21  1.9000e+01
19  22  1.3000e+01
20  23  7.0000e+00
```

Notice that the last couple of bins go over the end of the genome in a circular 
manner. We can also split the sliding window into bins by specifying the number
of bins to split into

```
python windowed_features.py --window_bins 3 --genome_length 22 sliding tests/test.gr 3 1
Start   End window_0    window_1    window_2
0   3   0.0000e+00  1.0000e+00  2.0000e+00
1   4   1.0000e+00  2.0000e+00  3.0000e+00
2   5   2.0000e+00  3.0000e+00  4.0000e+00
3   6   3.0000e+00  4.0000e+00  5.0000e+00
4   7   4.0000e+00  5.0000e+00  6.0000e+00
5   8   5.0000e+00  6.0000e+00  7.0000e+00
6   9   6.0000e+00  7.0000e+00  8.0000e+00
7   10  7.0000e+00  8.0000e+00  9.0000e+00
8   11  8.0000e+00  9.0000e+00  1.0000e+01
9   12  9.0000e+00  1.0000e+01  1.1000e+01
10  13  1.0000e+01  1.1000e+01  1.2000e+01
11  14  1.1000e+01  1.2000e+01  1.3000e+01
12  15  1.2000e+01  1.3000e+01  1.4000e+01
13  16  1.3000e+01  1.4000e+01  1.5000e+01
14  17  1.4000e+01  1.5000e+01  1.6000e+01
15  18  1.5000e+01  1.6000e+01  1.7000e+01
16  19  1.6000e+01  1.7000e+01  1.8000e+01
17  20  1.7000e+01  1.8000e+01  1.9000e+01
18  21  1.8000e+01  1.9000e+01  2.0000e+01
19  22  1.9000e+01  2.0000e+01  0.0000e+00
20  23  2.0000e+01  0.0000e+00  5.0000e-01
```

Thus giving us a column for each window. We can also add additional windows
upstream or downstream of the sliding window and split those as well.

```
python windowed_features.py --far_downstream 3 --far_downstream_bins 3 --far_upstream 3 --far_upstream_bins 3 --genome_length 21 sliding tests/test.gr 3 1
Start   End upstream_window_0   upstream_window_1   upstream_window_2   window_0    downstream_window_0 downstream_window_1 downstream_window_2
0   3   1.9000e+01  1.9500e+01  2.0000e+01  1.0000e+00  3.0000e+00  4.0000e+00  5.0000e+00
1   4   1.9500e+01  2.0000e+01  0.0000e+00  2.0000e+00  4.0000e+00  5.0000e+00  6.0000e+00
2   5   2.0000e+01  0.0000e+00  1.0000e+00  3.0000e+00  5.0000e+00  6.0000e+00  7.0000e+00
3   6   0.0000e+00  1.0000e+00  2.0000e+00  4.0000e+00  6.0000e+00  7.0000e+00  8.0000e+00
4   7   1.0000e+00  2.0000e+00  3.0000e+00  5.0000e+00  7.0000e+00  8.0000e+00  9.0000e+00
5   8   2.0000e+00  3.0000e+00  4.0000e+00  6.0000e+00  8.0000e+00  9.0000e+00  1.0000e+01
6   9   3.0000e+00  4.0000e+00  5.0000e+00  7.0000e+00  9.0000e+00  1.0000e+01  1.1000e+01
7   10  4.0000e+00  5.0000e+00  6.0000e+00  8.0000e+00  1.0000e+01  1.1000e+01  1.2000e+01
8   11  5.0000e+00  6.0000e+00  7.0000e+00  9.0000e+00  1.1000e+01  1.2000e+01  1.3000e+01
9   12  6.0000e+00  7.0000e+00  8.0000e+00  1.0000e+01  1.2000e+01  1.3000e+01  1.4000e+01
10  13  7.0000e+00  8.0000e+00  9.0000e+00  1.1000e+01  1.3000e+01  1.4000e+01  1.5000e+01
11  14  8.0000e+00  9.0000e+00  1.0000e+01  1.2000e+01  1.4000e+01  1.5000e+01  1.6000e+01
12  15  9.0000e+00  1.0000e+01  1.1000e+01  1.3000e+01  1.5000e+01  1.6000e+01  1.7000e+01
13  16  1.0000e+01  1.1000e+01  1.2000e+01  1.4000e+01  1.6000e+01  1.7000e+01  1.8000e+01
14  17  1.1000e+01  1.2000e+01  1.3000e+01  1.5000e+01  1.7000e+01  1.8000e+01  1.9000e+01
15  18  1.2000e+01  1.3000e+01  1.4000e+01  1.6000e+01  1.8000e+01  1.9000e+01  2.0000e+01
16  19  1.3000e+01  1.4000e+01  1.5000e+01  1.7000e+01  1.9000e+01  2.0000e+01  0.0000e+00
17  20  1.4000e+01  1.5000e+01  1.6000e+01  1.8000e+01  2.0000e+01  0.0000e+00  5.0000e-01
18  21  1.5000e+01  1.6000e+01  1.7000e+01  1.9000e+01  0.0000e+00  5.0000e-01  1.0000e+00
19  22  1.6000e+01  1.7000e+01  1.8000e+01  1.3000e+01  5.0000e-01  1.0000e+00  1.5000e+00
20  23  1.7000e+01  1.8000e+01  1.9000e+01  7.0000e+00  1.0000e+00  1.5000e+00  2.0000e+00
```

IMPORTANT: By default, any values below 0 or existing as a nan are converted to
a zero before being analyzed by the windower. To turn this behavior off specify 
`--no_truncate` to allow negative values and `--no_finite` to allow nans.

If the user wants to instead analyze a specific window (or set of windows)
they can specify a .gff file for the windows of interest. For example using the
following gff

```
test	test	test	10	15	.	+	.	test
```

We can use the windower to calculate the statistic of interest for this entire
window

```
python windowed_features.py --genome_length 21 gff_window tests/test.gr tests/test.gff 
Start   End window_0
9   15  1.1500e+01
```

We can also take the center of this window and pad around that center a specified
number of basepairs upstream or downstream

```
python windowed_features.py --genome_length 21 --upstream 3 --downstream 2 gff_window tests/test.gr tests/test.gff --center_metric median
Start   End window_0
9   15  1.0500e+01
```

However the window name in the example above will always be based on the original
gff entry. (Start and End are from the gff, not the calculated Start and End).

Finally we can also get information from the genomic sequence itself if we input
a fasta as input data.

```
>test                                                                               
AAAAAATGGCCCTTTTTTT 
```

To get the AT content of each window we would do the following:
```
python windowed_features.py --genome_length 21 --summary_stat genomic --genomic_feature basecontent "A;T" sliding tests/test.fasta 3 1
Start   End window_0
0   3   1.0000e+00
1   4   1.0000e+00
2   5   1.0000e+00
3   6   1.0000e+00
4   7   1.0000e+00
5   8   6.6667e-01
6   9   3.3333e-01
7   10  0.0000e+00
8   11  0.0000e+00
9   12  0.0000e+00
10  13  3.3333e-01
11  14  6.6667e-01
12  15  1.0000e+00
13  16  1.0000e+00
14  17  1.0000e+00
15  18  1.0000e+00
16  19  1.0000e+00
17  20  1.0000e+00
18  21  1.0000e+00
19  22  1.0000e+00
20  23  1.0000e+00
```

However to get the AT di-nucleotide frequency we would do the following:

```
python windowed_features.py --genome_length 21 --summary_stat genomic --genomic_feature basecontent "AT" sliding tests/test.fasta 3 1
Start   End window_0
0   3   0.0000e+00
1   4   0.0000e+00
2   5   0.0000e+00
3   6   0.0000e+00
4   7   1.0000e+00
5   8   1.0000e+00
6   9   0.0000e+00
7   10  0.0000e+00
8   11  0.0000e+00
9   12  0.0000e+00
10  13  0.0000e+00
11  14  0.0000e+00
12  15  0.0000e+00
13  16  0.0000e+00
14  17  0.0000e+00
15  18  0.0000e+00
16  19  0.0000e+00
17  20  0.0000e+00
18  21  0.0000e+00
19  22  0.0000e+00
20  23  0.0000e+00
```

Additionally, one can search with simple regular expression motifs and get the
count of motif matches with the following:

```
python windowed_features.py --genome_length 21 --summary_stat genomic --genomic_feature motif "[AC]T" sliding tests/test.fasta 3 1
Start   End window_0
0   3   0.0000e+00
1   4   0.0000e+00
2   5   0.0000e+00
3   6   0.0000e+00
4   7   1.0000e+00
5   8   1.0000e+00
6   9   0.0000e+00
7   10  0.0000e+00
8   11  0.0000e+00
9   12  0.0000e+00
10  13  1.0000e+00
11  14  1.0000e+00
12  15  0.0000e+00
13  16  0.0000e+00
14  17  0.0000e+00
15  18  0.0000e+00
16  19  0.0000e+00
17  20  0.0000e+00
18  21  0.0000e+00
19  22  0.0000e+00
20  23  0.0000e+00
```

To see a full listing of options type:

```
python windowed_features.py -h
usage: windowed_features.py [-h] [--window_bins WINDOW_BINS]
                            [--far_upstream FAR_UPSTREAM]
                            [--far_upstream_bins FAR_UPSTREAM_BINS]
                            [--far_downstream FAR_DOWNSTREAM]
                            [--far_downstream_bins FAR_DOWNSTREAM_BINS]
                            [--genome_length GENOME_LENGTH]
                            [--summary_stat SUMMARY_STAT]
                            [--genomic_feature GENOMIC_FEATURE GENOMIC_FEATURE]
                            [-o O] [--no_truncate] [--no_finite]
                            [--convert_logical CONVERT_LOGICAL]
                            [--upstream UPSTREAM] [--downstream DOWNSTREAM]
                            {sliding,gff_window} ...

Program to calculate a summary statistic over specified windows in a circular
genome. The options below will work on each command. To see additional help
for each command, specify 'command -h'. Both .gr files and .npy files are
assumed to be in 0-based coordinates. All output is also in 0-based
coordinates .gff files are assumed to be in 1-based coordinates

optional arguments:
  -h, --help            show this help message and exit
  --window_bins WINDOW_BINS
                        how many bins to break the window into
  --far_upstream FAR_UPSTREAM
                        number of bases upstream of a region, default=0
  --far_upstream_bins FAR_UPSTREAM_BINS
                        number of bins to split region into, default=1
  --far_downstream FAR_DOWNSTREAM
                        number of bases downstream of a region, default=0
  --far_downstream_bins FAR_DOWNSTREAM_BINS
                        number of bins to split region into, default=1
  --genome_length GENOME_LENGTH
                        length of the genome for coverage file, default=
                        4639676
  --summary_stat SUMMARY_STAT
                        summary stat to calculate for each window. [mean =
                        mean of the signal in the window, median = median of
                        the signal in the window, cutoff = if a number is
                        specified then return 1 if any datapoint is over the
                        cutoff, density = determine data density for a window
                        i.e. number of data points in the window genomic =
                        choosing a genomic feature, must specify
                        --genomic_feature flag] default= mean.
  --genomic_feature GENOMIC_FEATURE GENOMIC_FEATURE
                        When --summary_stat genomic is specified, this flag
                        MUST be specified Choose the type of genomic feature
                        desired: [motif regex gives count of motifs in the
                        window basecontent base1;base2;...etc. gives fraction
                        of window that is/are that/those base(s). Also works
                        with dinucleotide, tri... etc. Specifiying a mix of
                        tri/di/mono bases will give meaningless answers]
  -o O                  output file to put the output into (default stdout)
  --no_truncate         don't truncate values below zero to zero.
  --no_finite           don't truncate non-finite values to zero.
  --convert_logical CONVERT_LOGICAL
                        convert the data array to a logical array that is true
                        if the data point is above a cutoff value
  --upstream UPSTREAM   number of bases upstream of a region, default=0
  --downstream DOWNSTREAM
                        number of bases downstream of a region, default=0

window type commands. User must specify one of these.:
  Each command specifies a different way of defining windows

  {sliding,gff_window}
    sliding             do a sliding window over the genome
    gff_window          define windows from a gff
```

All optional arguments in the help above MUST be specified before specifying the
sliding or gff_window options

To see the individual options of each way of specifying windows type the 
following:

```
python windowed_features.py sliding -h
usage: windowed_features.py sliding [-h] [--plot_dist] [--name NAME]
                                    input_data size slide_by

positional arguments:
  input_data   input data. Accepted types include .gr, .npy
  size         size of the sliding window
  slide_by     number of basepairs to slide by

optional arguments:
  -h, --help   show this help message and exit
  --plot_dist  plot the distribution of the values in the main window
  --name NAME  [startend = put the start and end of the region, start = put
               only the starting bp, end = put only the ending bp, center =
               put the center of the region median(start,end)]
               default=startend
```

Or for the gff_windows:

```
python windowed_features.py gff_window -h
usage: windowed_features.py gff_window [-h] [--center_metric CENTER_METRIC]
                                       [--name NAME]
                                       [--distributions DISTRIBUTIONS]
                                       input_data gff_windows

positional arguments:
  input_data            input data. Accepted types include .gr, .npy
  gff_windows           gff file containing windows to scan

optional arguments:
  -h, --help            show this help message and exit
  --center_metric CENTER_METRIC
                        How to determine the center of the gff location?
                        [median = center is defined as the median, threeprime
                        = center is defined as the three prime end of the
                        region, fiveprime = center is defined as the five
                        prime end of the region, identity = entire region is
                        taken as a window + upstream and downstream bp]
                        Default = identity
  --name NAME           [comments= parse comments for gene name/bnumber,
                        fiveprime= just put the fiveprime bp, threeprime =
                        just put the threeprime bp, startend = put the start
                        and end of the region startendstrand = put the start,
                        end and strand of the region center = put the center
                        of the region median(start,end)] default=startend
  --distributions DISTRIBUTIONS
                        plot the distribution of the data set in the ranges
                        compared to randomly sample bins of the same length
```
