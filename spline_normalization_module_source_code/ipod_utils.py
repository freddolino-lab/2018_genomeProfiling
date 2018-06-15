#!/usr/bin/python

# A collection of functions useful in working with microarray data from IPOD

import numpy
import scipy
import scipy.stats
import pylab
import os
#import hcluster
import table
import scipy.interpolate as si
import scipy.signal
import pickle
import math
import rpy2.robjects as robjects
import growthcurves
import tempfile
import random
from tempfile import mktemp, mkstemp
import sklearn.decomposition

# path containing utilities like bpmap_bcel_join 
AFFYPATH = "/home/petefred/ST_research/ipod/from_tvora/util"

def check_file_empty(filename):
  """
  Return true if the file is empty
  """

  instr = open(filename, "r")
  firstline = instr.readline()
  if (firstline == ""):
    return True

  instr.close()
  return False

def read_peaklocs(filename):
  """
  Read a set of start,end pairs for regions of interest from a file

  Return a list of (start,end) tuples
  """

  instr = open(filename)
  peaklocs = []

  for line in instr:
    start,end = line.split()
    peaklocs.append( (int(start), int(end)) )

  instr.close()
  return peaklocs

def read_grfile(filename, skiprows=0):
  """
  Read data from a .gr file into a numpy array

  In the process we transpose the array, so that it has a[0] the indices and a[1] the values
  """

  locs = numpy.loadtxt(filename, dtype='int', usecols=(0,), skiprows=skiprows)
  vals = numpy.loadtxt(filename, usecols=(1,), skiprows=skiprows)
  return (locs,vals)

def grfiles_to_mat(filenames):
    """
    read multiple gr files into a single matrix

    the gr files must all have the same data length. the locations from the first specified file are used; the user is
      responsible for ensuring that all other given files have compatible locations

    returns a single locs array and a data matrix, where rows are locations and columns are files
    """

    if len(filenames) == 0:
        print "Must give at least one filename"
        return numpy.na
    elif len(filenames) == 1:
        return read_grfile(filenames[0])

    firstfile = filenames[0]
    locs = numpy.loadtxt(firstfile, dtype='int', usecols=(0,))
    vals = numpy.loadtxt(firstfile, usecols=(1,))

    datmat = numpy.zeros( (len(vals), len(filenames) ) )
    datmat[:,0] = vals
    nvals=len(vals)

    for i in range(1,len(filenames)):
        filename=filenames[i]
        vals = numpy.loadtxt(filename, usecols=(1,))
        if len(vals) != nvals:
            raise(ValueError("Number of entries in .gr file does not match the first one"))

        datmat[:,i] = vals


    return (locs, datmat)
    
    

def remove_nans(infile, outfile,zeros=False):
  """
  Write a new version of a .gr file with all nan entries removed.

  If zeros is True, we instead replace the nans with zeros
  """

  offsets, data = read_grfile(infile)
  goodinds = scipy.isfinite(data)

  if zeros:
    data[numpy.invert(goodinds)] = 0.0000
    write_grfile(offsets, data, outfile)

  else:
    write_grfile(offsets[goodinds], data[goodinds], outfile)


def write_grfile(indices, data, filename, header=None):
  """
  Write a gr file given the list of sequence positions and corresponding scalars
  """

  ostr = open(filename, 'w')
  if header is not None:
    ostr.write(header)
  for ind, dat in zip(indices, data):
    ostr.write("%i %f\n" % (ind, dat))
  ostr.close()

def write_int_grfile(indices, data, filename, header=None):
  """
  Write a gr file given the list of sequence positions and corresponding integers
  """

  ostr = open(filename, 'w')
  if header is not None:
    ostr.write(header)
  for ind, dat in zip(indices, data):
    ostr.write("%i %i\n" % (ind, dat))
  ostr.close()

def write_grbin(indices, data, filename):
  """
  Write a binary form of a gr file
  """

  numpy.savez(filename, indices, data)

def write_pmmmfile(indices, pms, mms, filename):
  """
  Write a set of position/pm/mm data
  """

  ostr = open(filename, 'w')
  for ind, pm, mm in zip(indices, pms, mms):
    ostr.write("%i %i %i\n" % (ind, pm, mm))
  ostr.close()

def convert_cel_to_text(celfile, bpmapfile, outfile, flags="-fieldMap"):
  """
  Use affy's bpmap_bcel_join program to combine sequence data and intensity data

  celfile and bpmapfile should be, respectively, the raw .CEL file from an
  array scan and the file containing the sequence coordinates of each spot

  Output will be written to outfile. Any extra flags for bpmap_bcel_join can be specified with the flags optional argument
  """

  os.system("%s/bpmap_bcel_join -bpmap %s -bcel %s %s > %s" % (AFFYPATH, bpmapfile, celfile, flags, outfile))

def clean_text(txtfile, outfile):
  """
  Remove unnecessary lines and fields from a converted cel->text file

  We write the tag, position, pm, and mm fields
  """

  instr = open(txtfile, 'r')
  ostr = open(outfile, 'w')
  
  for line in instr:
    if (line[0] == "#"):
      continue

    flag, tag, pos, pm, mm = line.split()
    ostr.write("%s %s %s %s\n" % (tag, pos, pm, mm))

  instr.close()
  ostr.close()


def clean_and_correct_text(txtfile, outfile, clamp=None):
  """
  Remove unnecessary lines and fields from a converted cel->text file
  
  In the process, we replace each line's PM/MM values with the single value
  PM-MM

  We also print the tags and positions

  If clamp is not None, we set any value less than clamp to clamp
  """

  instr = open(txtfile, 'r')
  ostr = open(outfile, 'w')
  
  for line in instr:
    if (line[0] == "#"):
      continue

    flag, tag, pos, pm, mm = line.split()
    newval = float(pm) - float(mm)

    if ((clamp is not None) and (newval < clamp)):
      newval = clamp

    ostr.write("%s %s %i\n" % (tag, pos, newval))

  instr.close()
  ostr.close()

def convert_txt_to_gr(txtfile, outfile):
  """
  Convert a text file including an extra first field into a gr file
  """

  instr = open(txtfile, 'r')
  ostr = open(outfile, 'w')

  for line in instr:
    linearr = line.split()
    ostr.write("%s %s\n" % (linearr[1], linearr[2]))

  instr.close()
  ostr.close()

def sort_gr_file(txtfile, outfile):
  """
  Sort a .gr file by sequence index
  """

  offsets, data = read_grfile(txtfile)
  print "Doing zip"
  valpairs = zip(offsets, data)
  print "Doing sort"
  valpairs.sort(lambda x, y : cmp(x[0], y[0]))
  print "Splitting data"
  offsets = [x[0] for x in valpairs]
  vals = [x[1] for x in valpairs]
  print "Writing data"
  write_grfile(offsets, vals, outfile)
  

  

def extract_genome_re(txtfile, outfile, genome_name="E_coli_genome"):
  """
  Write a gr file from a genome/offset/value text file

  Only entries with a genome name corresponding to the regular expression
    genome_name will be written
  """

  import re

  instr = open(txtfile, 'r')
  ostr = open(outfile, 'w')

  for line in instr:
    linearr = line.split()
    if (re.search(genome_name, linearr[0])):
      ostr.write("%s %s\n" % (linearr[1], linearr[2]))

  instr.close()
  ostr.close()

def extract_genome_flag(txtfile, outfile, genome_name="E_coli_genome"):
  """
  Write a gr file from a genome/offset/value text file

  Only entries with a genome name corresponding to genome_name 
    will be written
  """
  instr = open(txtfile, 'r')
  ostr = open(outfile, 'w')

  for line in instr:
    linearr = line.split()
    if (linearr[0] == genome_name):
      ostr.write("%s %s\n" % (linearr[1], linearr[2]))

  instr.close()
  ostr.close()

def extract_seqgroup(infile, outfile, seqgroup):
  """
  Extract all entries from a Affy text file in the specified seqgroup
  """

  instr = open(infile)
  ostr = open(outfile, 'w')

  goodblock = False

  for line in instr:
    if (line[0] == "#"):
      linearr = line.split()
      if (linearr[0] == "#seq_group_name"):
        goodblock = (linearr[1] == seqgroup)
    else:
      if (goodblock):
        ostr.write(line)

  ostr.close()

def subtract_percentile_from_genomefile(infile, bgfile, outfile, perc=50):
  """
  Identical to subtract_percentile_from_pmmm_file, but with genome information

  The input should contain genome/offset/pm/mm, and the same will be written
  """

  import scipy.stats

  mydat = numpy.loadtxt(infile, unpack=True, dtype="|S8")
  bgdat = numpy.loadtxt(bgfile, unpack=True, dtype="|S8")

  numtypes, numpoints = mydat.shape

  bg_vals = (bgdat[2:]).astype(int)
  my_vals = (mydat[2:]).astype(int)

  bg_val = scipy.stats.scoreatpercentile( scipy.reshape(bg_vals, (-1,1)) , perc )
  for i in (0,1):
    tmparr = my_vals[i,:]
    tmparr -= bg_val
    tmparr += (-1 * tmparr * scipy.less(tmparr,0))
    my_vals[i,:] = tmparr

  ostr = open(outfile, "w")

  for i in range(numpoints):
    ostr.write("%s %s %i %i\n" % (mydat[0,i], mydat[1,i], my_vals[0,i], my_vals[1,i]))

  ostr.close()


def subtract_percentile_from_pmmm_file(infile, bgfile, outfile, perc=50):
  """
  Subtract the perc-percentile of the combined pm/mm probes in bgfile from infile

  Values will be clamped on the bottom end at zero

  We assume that the output contains offset/pm/mm, and write the same
  """

  import scipy.stats

  mydat = numpy.loadtxt(infile, unpack=True)
  bgdat = numpy.loadtxt(bgfile, unpack=True)

  bg_val = scipy.stats.scoreatpercentile( scipy.reshape(bgdat[1:], (-1,1)) , perc )
  for i in (1,2):
    tmparr = mydat[i,:]
    tmparr -= bg_val
    tmparr += (-1 * tmparr * scipy.less(tmparr,0))
    mydat[i,:] = tmparr

  write_pmmmfile(mydat[0], mydat[1], mydat[2], outfile)

def subtract_mm_from_pm_clamped(infile, outfile):
  """
  Subtract the mm values from the pm values in infile and write to outfile

  Negative values will be replaced with zeros
  """

  mydat = numpy.loadtxt(infile, unpack=True)

  offsets, pm, mm = mydat
  pmmm = pm - mm
  pmmm += (-1 * pmmm * scipy.less(pmmm, 0))

  write_grfile(offsets, pmmm, outfile)

def remove_mm_probes(infile, outfile):
  """
  Remove all the mm values from offset/pm/mm files
  """

  mystr = open(infile, 'r')
  outstr = open(outfile, 'w')

  for line in mystr:
    linearr = line.split()
    outstr.write("%s %s\n" % (linearr[0], linearr[1]))

  mystr.close()
  outstr.close()

def calc_ratio(firstfile, secondfile, outfile, assume_common=True,minval=0.0):
  """
  Calculate and write the ratio at each probe position

  We ignore any offsets not present in both files

  If assume_common is true, we don't bother checking that the offsets are the same

  minval should be small relative to the data, and is added to all
    data points to avoid division by zero
  """

  offsets1, pm1 = numpy.loadtxt(firstfile, unpack=True)

  offsets2, pm2 = numpy.loadtxt(secondfile, unpack=True)

  pm1 += minval
  pm2 += minval

  if (assume_common):
    pmratio = pm1 / pm2
    write_grfile(offsets1, pmratio, outfile)
    return

  all_offsets = set.intersection(set(offsets1), set(offsets2))
  all_offsets = list(all_offsets)
  all_offsets.sort()

  offsets1l = list(offsets1)
  offsets2l = list(offsets2)

  # hold the probe values corresponding to common offsets
  goodvals1 = []
  goodvals2 = []

  for off in all_offsets:
    goodvals1.append( pm1[ offsets1l.index(off) ] )
    goodvals2.append( pm2[ offsets2l.index(off) ] )

  pm1_good = scipy.array(goodvals1)
  pm2_good = scipy.array(goodvals2)

  write_grfile(all_offsets, pm1_good / pm2_good, outfile)

def calc_logratio(firstfile, secondfile, outfile, assume_common=True,minval=None):
  """
  Calculate and write the log2ratio at each probe position

  We ignore any offsets not present in both files

  If assume_common is true, we don't bother checking that the offsets are the same

  If minval is specified, it is added to all values prior to taking the log (to avoid nans)
  """

  offsets1, pm1 = numpy.loadtxt(firstfile, unpack=True)

  offsets2, pm2 = numpy.loadtxt(secondfile, unpack=True)

  if (minval is not None):
    pm1 += minval
    pm2 += minval

  if (assume_common):
    pmratio = pm1 / pm2
    write_grfile(offsets1, scipy.log2(pmratio), outfile)
    return

  all_offsets = set.intersection(set(offsets1), set(offsets2))
  all_offsets = list(all_offsets)
  all_offsets.sort()

  offsets1l = list(offsets1)
  offsets2l = list(offsets2)

  # hold the probe values corresponding to common offsets
  goodvals1 = []
  goodvals2 = []

  for off in all_offsets:
    goodvals1.append( pm1[ offsets1l.index(off) ] )
    goodvals2.append( pm2[ offsets2l.index(off) ] )

  pm1_good = scipy.array(goodvals1)
  pm2_good = scipy.array(goodvals2)

  write_grfile(all_offsets, scipy.log2( pm1_good / pm2_good), outfile)


def subtract_mm_from_pm_withgenome(infile, outfile):
  """
  Subtract the mm values from the pm values in infile and write to outfile

  Negative values will be replaced with zeros

  Here we assume that the file has genome/offset/pm/mm values

  We then write genome/offset/pm-mm
  """

  mydat = numpy.loadtxt(infile, unpack=True, dtype="|S8")
  numblocks, numentries = mydat.shape

  offsets, pm, mm = (mydat[1:].astype(int))
  pmmm = pm - mm
  pmmm += (-1 * pmmm * scipy.less(pmmm, 0))

  ostr = open(outfile, "w")

  for i in range(numentries):
    ostr.write("%s %s %i\n" % (mydat[0,i], mydat[1,i], pmmm[i]))

  ostr.close()


def extract_genome_flag_pmmm(txtfile, outfile, genome_re = "E_coli_genome"):
  """
  Write a file containing offset/pm/mm from a file containing genome/offset/pm/mm

  Only lines where the genome matches genome_re will be written
  """

  import re

  myre = re.compile(genome_re)

  instr = open(txtfile, 'r')
  ostr = open(outfile, 'w')

  for line in instr:
    linearr = line.split()
    if (myre.search(linearr[0])):
      ostr.write("%s %s %s\n" % (linearr[1], linearr[2], linearr[3]))

  instr.close()
  ostr.close()

def gauss_file_ipod(txtfile, outfile, gausswidth=1, width=11, dropvals=0):
  """
  Apply a gaussian convolution of specified width to input text data
  
  We assume .gr input and output

  Note that width is the total size of the window that will be considered in
  the convolution, and should be wide enough so that the values beyond it
  are negligable for the chosen value of gausswidth

  Note that all widths are in indices, **not** bp. This function
  should not be used unless the spacings between probes are uniform.

  If dropvals is nonzero, the DROPVALS highest and lowest points in
  each window are ignored, and the resulting sum will be scaled
  according to the sum of the remaining weights after masking these
  values
  """

  import scipy.signal
  import scipy.stats
  myxvals = scipy.arange(-(width/2),(width+1)/2)
  kernel = [scipy.stats.norm.pdf(x,scale=gausswidth) for x in myxvals]

  # test to make sure out window size is reasonable
  print "Performing gaussian convolution with a %i-element kernel" % len(kernel)
  print "Integral of the kernel is %f" % scipy.integrate.trapz(kernel, myxvals)

  if (width % 2 == 0):
    print "WARNING! Convolutions with an even width are not a good idea"

  instr1 = open(txtfile, 'r')
  instr2 = open(txtfile, 'r')
  ostr = open(outfile, 'w')

  # first figure out the total number of lines
  start = instr1.tell()
  nr = 0
  currline = " "
  while(currline != ""):
    nr += 1
    currline = instr1.readline()

  instr1.seek(start)

  vallist = []; # list of the current values for the average
  # read the values for the wrapped end of the chromosome
  for i in range(nr - ( (width+1)/2)):
    instr1.readline()

  currline = instr1.readline()
  while (currline != ""):
    index, val = currline.split()
    vallist.append(float(val))
    currline = instr1.readline()

  instr1.seek(start)

  print "Back half of the array has %i elements" % len(vallist)


  # use instr1 to load up the sets of values up to the end of the first
  #   moving average
  for i in range((width+1)/2):
    currline = instr1.readline()
    index, val = currline.split()
    vallist.append(float(val))

  print "Calculating running average with %i elements" % len(vallist)
  #ostr.write("Current list of elements for vallist: \n")
  #for elem in vallist:
  #  ostr.write("%s\n" % elem)
  #ostr.write("----------\n")
  if (dropvals == 0):
    usekern = kernel
  else:
    numdropped = 0
    while (numdropped < dropvals):
      min_ind, max_ind = find_min_max_ind(vallist)
      usekern = numpy.array(kernel, copy=False)
      usekern[min_ind] = 0
      usekern[max_ind] = 0
      vallist[min_ind] = numpy.median(vallist)
      vallist[max_ind] = numpy.median(vallist)

  curravg = numpy.dot(vallist, usekern) / numpy.sum(usekern)

  # now start walking through the array and writing the correct values to ostr
  currline = instr2.readline()
  lineind = 0
  while (currline != ""):
    index, val = currline.split()
    ostr.write("%i %f\n" % ( int(index), curravg ) )
    newline = instr1.readline()

    # loop around the chromosome if needed
    if (newline == ""):
      instr1.seek(start)
      newline = instr1.readline()
    index, val = newline.split()
    newval = float(val)

    vallist.pop(0)
    vallist.append(newval)

    if (dropvals == 0):
      usekern = kernel
    else:
      numdropped = 0
      usekern = numpy.array(kernel, copy=False)
      while (numdropped < dropvals):
        min_ind, max_ind = find_min_max_ind(vallist)
        usekern[min_ind] = 0
        usekern[max_ind] = 0
        vallist[min_ind] = numpy.median(vallist)
        vallist[max_ind] = numpy.median(vallist)

    curravg = numpy.dot(vallist, usekern) / numpy.sum(usekern)

    currline = instr2.readline()
    lineind += 1

  instr1.close()
  instr2.close()
  ostr.close()


def do_runningavg_ipod(txtfile, outfile, width=20001):
  """
  Calculate the running average over a given length on the input text data
  
  We assume the input format specified by extract_genome_flag
  Note that a .gr file is written -- ie, the "genome" field is ignored

  This function is **not** safe if there are any gaps in bp numbering
  """

  # for anything more complicated than a running average,
  # it will pay to use a convolution
  #import scipy.signal
  #kernel = scipy.array([1.0/width for i in range(width)])

  if (width % 2 == 0):
    print "WARNING! Convolutions with an even width are not a good idea"

  instr1 = open(txtfile, 'r')
  instr2 = open(txtfile, 'r')
  ostr = open(outfile, 'w')

  # first figure out the total number of lines
  start = instr1.tell()
  nr = 0
  currline = " "
  while(currline != ""):
    nr += 1
    currline = instr1.readline()

  instr1.seek(start)

  vallist = []; # list of the current values for the average
  # read the values for the wrapped end of the chromosome
  for i in range(nr - ( (width+1)/2)):
    instr1.readline()

  currline = instr1.readline()
  while (currline != ""):
    index, val = currline.split()
    vallist.append(float(val))
    currline = instr1.readline()

  instr1.seek(start)

  print "Back half of the array has %i elements" % len(vallist)


  # use instr1 to load up the sets of values up to the end of the first
  #   moving average
  for i in range((width+1)/2):
    currline = instr1.readline()
    index, val = currline.split()
    vallist.append(float(val))

  print "Calculating running average with %i elements" % len(vallist)
  #ostr.write("Current list of elements for vallist: \n")
  #for elem in vallist:
  #  ostr.write("%s\n" % elem)
  #ostr.write("----------\n")
  curravg = scipy.mean(vallist)

  # now start walking through the array and writing the correct values to ostr
  currline = instr2.readline()
  lineind = 0
  while (currline != ""):
    index, val = currline.split()
    ostr.write("%i %f\n" % ( int(index), curravg ) )
    curravg -= vallist[0] / width
    newline = instr1.readline()

    # loop around the chromosome if needed
    if (newline == ""):
      instr1.seek(start)
      newline = instr1.readline()
    index, val = newline.split()
    newval = float(val)
    curravg += newval / width

    vallist.pop(0)
    vallist.append(newval)

    currline = instr2.readline()
    lineind += 1

  instr1.close()
  instr2.close()
  ostr.close()

def do_runningmedian_opt(txtfile, outfile, width=49, genomelength=None):
  """
  Do a running median over a given width (in bp)

  This function is safe to use with jumps in bp numbering
  """

  import bisect

  if (width % 2 == 0):
    print "WARNING! Convolutions with an even width are not a good idea"

  halfwidth = width / 2

  instr1 = open(txtfile, 'r')
  instr2 = open(txtfile, 'r')
  ostr = open(outfile, 'w')
  # first figure out the total number of lines
  start = instr1.tell()
  nr = 0
  currline = instr1.readline()
  firstind, firstval = currline.split()
  firstind = int(firstind)
  finalind, finalval = (0,0)
  while(currline != ""):
    nr += 1
    if (currline != " "):
      lastind, lastval = (finalind, finalval)
      finalind, finalval = currline.split()
    currline = instr1.readline()

  if not genomelength:
    genomelength = int (finalind) + int(finalind) - int(lastind)
    print "Guessing genome length of %i" % genomelength

  instr1.seek(start)

  vallist = []; # list of the current values for the average
  offsetlist = []; # corresponding offsets

  # read the values for the wrapped end of the chromosome
  # Note that we search out the appropriate entry based on the difference in bp
  for i in range(nr - ( (width+1)/2)):
    instr1.readline()

  firstbp = firstind - halfwidth
  lastbp = firstind + halfwidth

  currline = instr1.readline()
  while (currline != ""):
    index, val = currline.split()
    index = int(index)
    val = float(val)
    if index > firstbp:
      target_ind = bisect.bisect_right(vallist, val)
      vallist.insert(target_ind, val)
      offsetlist.insert(target_ind, index)
    currline = instr1.readline()

  instr1.seek(start)

  # use instr1 to load up the sets of values up to the end of the first
  #   moving average
  currline = instr1.readline()
  index, val = currline.split()
  index = int(index)
  val = float(val)
  while (index < lastbp):
    target_ind = bisect.bisect_right(vallist, val)
    vallist.insert(target_ind, val)
    offsetlist.insert(target_ind, index)
    currline = instr1.readline()
    index, val = currline.split()
    index = int(index)
    val = float(val)

  #print "Calculating running average with %i elements" % len(vallist)
  #ostr.write("Current list of elements for vallist: \n")
  #for offset, elem in zip(offsetlist, vallist):
  #  ostr.write("%s %s\n" % (offset,elem))
  #ostr.write("----------\n")
  curravg = sorted_list_median(vallist)

  # now start walking through the array and writing the correct values to ostr
  curr_center_line = instr2.readline()

  while (curr_center_line != ""):
    center_offset, center_val = curr_center_line.split()
    center_offset = int(center_offset)
    index, val = currline.split()
    index = int(index)
    val = float(val)

    goodrange = [center_offset - halfwidth, center_offset + halfwidth]
    #print "Centered on %s; range is %s - %s" % (curr_center_line, goodrange[0], goodrange[1])

    if goodrange[0] < 0:
      prevrange = [genomelength + goodrange[0], genomelength]
      goodrange[0] = 0
    elif goodrange[1] > genomelength:
      prevrange = [0, goodrange[1] - genomelength]
      goodrange[1] = genomelength
    else:
      prevrange = [-1,-1]

    #print "Acceptable ranges: %s %s" % (goodrange, prevrange)

    # remove unneeded elements
    i=0
    while (i < len(offsetlist)):
      curroffset = offsetlist[i]
      if (curroffset <= goodrange[1] and curroffset >= goodrange[0]):
        i += 1
        continue
      if (curroffset <= prevrange[1] and curroffset >= prevrange[0]):
        i += 1 
        continue
      offsetlist.pop(i)
      vallist.pop(i)
      
    if (index > genomelength):
        instr1.seek(start)
      
    #print "New index start: %i" % index
    while ( (index >= goodrange[0] and index <= goodrange[1]) or (index >= prevrange[0] and index <= prevrange[1])):


      newloc = bisect.bisect_right(vallist, val)
      vallist.insert(newloc, val)
      offsetlist.insert(newloc, index)

      currline = instr1.readline()
      #print "Current line: %s" % currline
      if (currline == ""):
        instr1.seek(start)
        currline = instr1.readline()
        #print "Wrapping to beginning"
      index, val = currline.split()
      index = int(index)
      val = float(val)

    #print "For center offset %i, lists are %s || %s" % (center_offset, vallist, offsetlist)
    #curravg = scipy.mean(vallist)
    ostr.write("%i %f\n" % ( int(center_offset), sorted_list_median(vallist) ) )

    curr_center_line = instr2.readline()

  instr1.close()
  instr2.close()
  ostr.close()


def do_runningavg_opt(txtfile, outfile, width=49, genomelength=None):
  """
  Do a running average over a given width (in bp)

  This is similar to running convolve_file with kerneltype=runavg, but
    is faster due to optimizations specific for doing running averages

  This function **is** safe to jumps in bp numbering
  """

  
  # for anything more complicated than a running average,
  # it will pay to use a convolution
  #import scipy.signal
  #kernel = scipy.array([1.0/width for i in range(width)])

  if (width % 2 == 0):
    print "WARNING! Convolutions with an even width are not a good idea"

  halfwidth = width / 2

  instr1 = open(txtfile, 'r')
  instr2 = open(txtfile, 'r')
  ostr = open(outfile, 'w')

  # first figure out the total number of lines
  start = instr1.tell()
  nr = 0
  currline = instr1.readline()
  firstind, firstval = currline.split()
  firstind = int(firstind)
  finalind, finalval = (0,0)
  while(currline != ""):
    nr += 1
    if (currline != " "):
      lastind, lastval = (finalind, finalval)
      finalind, finalval = currline.split()
    currline = instr1.readline()


  if not genomelength:
    genomelength = int (finalind) + int(finalind) - int(lastind)
    print "Guessing genome length of %i" % genomelength

  instr1.seek(start)

  vallist = []; # list of the current values for the average
  offsetlist = []; # corresponding offsets

  # read the values for the wrapped end of the chromosome
  # Note that we search out the appropriate entry based on the difference in bp
  for i in range(nr - ( (width+1)/2)):
    instr1.readline()

  firstbp = firstind - halfwidth
  lastbp = firstind + halfwidth

  currline = instr1.readline()
  while (currline != ""):
    index, val = currline.split()
    index = int(index)
    if index > firstbp:
      vallist.append(float(val))
      offsetlist.append(int(index))
    currline = instr1.readline()

  instr1.seek(start)

  # use instr1 to load up the sets of values up to the end of the first
  #   moving average
  currline = instr1.readline()
  index, val = currline.split()
  index = int(index)
  while (index < lastbp):
    vallist.append(float(val))
    offsetlist.append(index)
    currline = instr1.readline()
    index, val = currline.split()
    index = int(index)

  #print "Calculating running average with %i elements" % len(vallist)
  #ostr.write("Current list of elements for vallist: \n")
  #for offset, elem in zip(offsetlist, vallist):
  #  ostr.write("%s %s\n" % (offset,elem))
  #ostr.write("----------\n")
  curravg = scipy.mean(vallist)
  currlen = len(vallist)

  # now start walking through the array and writing the correct values to ostr
  curr_center_line = instr2.readline()

  while (curr_center_line != ""):
    center_offset, center_val = curr_center_line.split()
    center_offset = int(center_offset)
    index, val = currline.split()
    index = int(index)
    val = float(val)

    goodrange = [center_offset - halfwidth, center_offset + halfwidth]
    #print "Centered on %s; range is %s - %s" % (curr_center_line, goodrange[0], goodrange[1])

    if goodrange[0] < 0:
      prevrange = [genomelength + goodrange[0], genomelength]
      goodrange[0] = 0
    elif goodrange[1] > genomelength:
      prevrange = [0, goodrange[1] - genomelength]
      goodrange[1] = genomelength
    else:
      prevrange = [-1,-1]

    #print "Acceptable ranges: %s %s" % (goodrange, prevrange)

    # find how big a slice of the list we need to cut
    removeind = -1
    for i in range(len(offsetlist)):
      curroffset = offsetlist[i]
      if (curroffset <= goodrange[1] and curroffset >= goodrange[0]):
        break
      if (curroffset <= prevrange[1] and curroffset >= prevrange[0]):
        break
      removeind = i
      
    #print removeind
    #print currlen
    #print vallist
    #print offsetlist
    #print "---"
    if (removeind >= 0):
      removeind += 1
      removedvals = vallist[:removeind]
      vallist = vallist[removeind:]
      offsetlist = offsetlist[removeind:]

      newlen = currlen - removeind
      if (newlen == 0):
        currlen = 0
        curravg = 0

      else:
        curravg = (currlen/float(newlen)) * curravg - numpy.sum(removedvals) / float(newlen)
        currlen=newlen

    # Add in new values as needed
    #if (index > genomelength):
    #  index -= genomelength

    if (index > genomelength):
        instr1.seek(start)
      
    #print "New index start: %i" % index
    while ( (index >= goodrange[0] and index <= goodrange[1]) or (index >= prevrange[0] and index <= prevrange[1])):
      vallist.append(val)
      newlen = currlen + 1
      curravg = curravg * (currlen / float(newlen)) + val / float(newlen)
      currlen = newlen
      offsetlist.append(index)
      currline = instr1.readline()
      #print "Current line: %s" % currline
      if (currline == ""):
        instr1.seek(start)
        currline = instr1.readline()
        #print "Wrapping to beginning"
      index, val = currline.split()
      index = int(index)
      val = float(val)

    #print "For center offset %i, lists are %s || %s" % (center_offset, vallist, offsetlist)
    #curravg = scipy.mean(vallist)
    ostr.write("%i %f\n" % ( int(center_offset), curravg ) )

    curr_center_line = instr2.readline()

  instr1.close()
  instr2.close()
  ostr.close()


def convolve_core(offsets, data, gausswidth=4, width=49, genomelength=None, kerneltype="gauss", kernel=None, wrap=True):
  """
  Apply a normalized discrete convolution of specified width to a set of locations and values
  
  Note that width is the total size of the window that will be considered in
  the convolution, and should be wide enough so that the values beyond it
  are negligable for the chosen value of gausswidth

  Note that all widths are in terms of bp; there may be relatively few
  entries in a given bp range. This function is safe to use with uneven
  spacing between probes.

  By default, a gaussian convolution is used. If kerneltype is instead 
    "runavg", a running average over a window of size width will
    be performed instead, and gausswidth will be ignored

  Wrapping will only be performed if "wrap" is true; otherwise we only
    use whatever valid values are available near the ends

  If kernel is defined, it supercedes gausswidth, width, and kerneltype. It must be a single list of floats containing the weights,
    going bp by bp, and centered on the center of the distribution. 

  returns the convolved array
  """

  import scipy.signal
  import scipy.stats
  myxvals = scipy.arange(-(width/2),(width+1)/2)

  output = numpy.zeros_like(data)


  if (width % 2 == 0):
    print "ERROR! Convolution width must be odd"
    raise(ValueError)

  if (kernel is not None):
    width = len(kernel)
    if width % 2 == 0:
      raise(ValueError("Kernel must have an odd number of elements"))

  #  raise("explicit kernel not yet implemented")
    kernel = numpy.array(kernel)

  elif (kerneltype == "gauss"):
    kernel = numpy.array([scipy.stats.norm.pdf(x,scale=gausswidth) for x in myxvals])

    # test to make sure out window size is reasonable
    print "Performing gaussian convolution with a %i-element kernel" % len(kernel)
    print "Integral of the kernel is %f" % scipy.integrate.trapz(kernel, myxvals)
  elif (kerneltype == "runavg"):
    print "Doing running average of width %i" % width
    kernel = scipy.array([1.0 for x in myxvals])
  else:
    print "ERROR: Unknown kernel type specified"
    raise(ValueError)

  halfwidth = width / 2

  if not(genomelength):
    genomelength = offsets[-1] + (offsets[-1] - offsets[-2])

  guess_start = 0
  guess_end = 1
  for off_ind,offset in enumerate(offsets):

    if (offset % 1000 == 0):
      print "Working on offset %i" % offset

    # Get the region of the data that we need to work on
    offset_slice, data_slice, guess_start, guess_end = circular_range_bps(offsets, data, offset-halfwidth, offset+halfwidth, genomelength=genomelength, guess_start=guess_start, guess_end=guess_end, return_inds=True)
    #print " Correct guesses for center %i would have been %i/%i" % (offset, guess_start, guess_end)
    #print "  Offsets: %i -- %i -- %i" % (offsets[guess_start], offset, offsets[guess_end])
    #print "Working on slice centered at offset %i (%i -- %i). Our offsets and data are %s and %s" % (offset, offset-halfwidth, offset+halfwidth, offset_slice, data_slice)

    if (wrap):
      if (offset < (genomelength/2)):
        offset_slice -= genomelength * scipy.greater(offset_slice, genomelength-width)
      else:
        offset_slice += genomelength * scipy.less(offset_slice, width)
    else:
      goodinds = scipy.all(scipy.vstack((scipy.greater(offset_slice, -1), scipy.less(offset_slice, genomelength))), axis=0)
      offset_slice = offset_slice[goodinds]
      data_slice = data_slice[goodinds]

    #print offset_slice
    offset_slice -= offset
    offset_slice += len(kernel)/2

    #print "Our new offsets corresponding to the kernel are %s" % (offset_slice)
    #print len(offset_slice)


    total = 0
    total_weights = 0
    #print offset_slice
    #print data_slice
    #print "----"
    #for (myoff, mydat) in zip(offset_slice, data_slice):
    #  total += kernel[myoff] * mydat
    #  total_weights += kernel[myoff]
    total = numpy.dot(kernel[offset_slice],data_slice)
    total_weights = numpy.sum(kernel[offset_slice])


    output[off_ind] = (total / total_weights)
    #raise("")


  return output


def convolve_file(txtfile, outfile, **kwargs):
  """
  Apply a normalized discrete convolution of specified width to input text data
  
  arguments are almost entirely passed to convolve_core
  """


  offsets, data = read_grfile(txtfile)
  outdat = convolve_core(offsets, data, **kwargs)
  write_grfile(offsets, outdat, outfile)

def write_ref_mean_std(filenames, outfile, outmeans=False):
  """
  Calculate and write the mean and stdev of a set of hyb data

  filenames should contain one or more files containing a .gr-formatted set
  of hyb data

  If outmeans is not False, write another file contains only the means

  If outfile is False, the main output is not written
  """

  offsets, vals = numpy.loadtxt(filenames[0], unpack=True)

  for name in filenames[1:]:
    newvals = numpy.loadtxt(name, usecols=(1,), unpack=True)
    vals = scipy.vstack((vals,newvals))

  vals = numpy.ma.masked_array(vals, numpy.invert(numpy.isfinite(vals)))
  #print scipy.shape(vals)
  means = numpy.ma.mean(vals, axis=0)
  stdevs = numpy.ma.std(vals, axis=0, ddof=1)

  # Fix cases with only one sample in each case

  # Fix cases with no samples
  #sortedvals = numpy.ma.sort(means)
  #numgood = numpy.ma.count(sortedvals)
  #second_perc = sortedvals[int(numgood * 0.02)]
  #fiftieth_perc = sortedvals[int(numgood*0.5)]
  #means[goodcounts == 0] = second_perc
  #stdevs[goodcounts == 0] = fiftieth_perc
  #stdevs[means < 0.001] = fiftieth_perc


  if (outfile):
    ostr = open(outfile, "w")
    for offset, mean, std in zip(offsets, means, stdevs):
      ostr.write("%i %f %f\n" % (offset, mean, std))
    ostr.close()


  if (outmeans):
    ostr = open(outmeans, "w")
    for offset, mean in zip(offsets, means):
      ostr.write("%i %f\n" % (offset, mean))
    ostr.close()

def write_ref_log_mean_std(filenames, outfile, outmeans=False):
  """
  Calculate and write the mean and stdev of a log-transformed set of hyb data

  filenames should contain one or more files containing a .gr-formatted set
  of hyb data

  If outmeans is not False, write another file contains only the means

  If outfile is False, the main output is not written

  All operations will be performed on the log2 of the input data

  We make two additional corrections:
    -if there is only one non-nan value at a given probe, we set sigma to
      the value at that probe
    -if there are *no* non-nan values at a probe, we set the value to the 2nd percentile
      of all data, and the sigma to the 50th percentile
  """

  offsets, vals = numpy.loadtxt(filenames[0], unpack=True)
  vals = numpy.log2(vals)

  for name in filenames[1:]:
    newvals = numpy.log2(numpy.loadtxt(name, usecols=(1,), unpack=True))
    vals = scipy.vstack((vals,newvals))

  vals = numpy.ma.masked_array(vals, numpy.invert(numpy.isfinite(vals)))
  #print scipy.shape(vals)
  means = numpy.ma.mean(vals, axis=0)
  stdevs = numpy.ma.std(vals, axis=0, ddof=1)

  # Fix cases with only one sample in each case
  #goodcounts = numpy.ma.count(vals, axis=0)
  #stdevs[goodcounts == 1] = means[goodcounts == 1]

  # Fix cases with no samples
  #sortedvals = numpy.ma.sort(means)
  #numgood = numpy.ma.count(sortedvals)
  #second_perc = sortedvals[int(numgood * 0.02)]
  #fiftieth_perc = sortedvals[int(numgood*0.5)]
  #means[goodcounts == 0] = second_perc
  #stdevs[goodcounts == 0] = fiftieth_perc

  #print means[32550]
  #print stdevs[32550]

  if (outfile):
    ostr = open(outfile, "w")
    for offset, mean, std in zip(offsets, means, stdevs):
      ostr.write("%i %f %f\n" % (offset, mean, std))
    ostr.close()


  if (outmeans):
    ostr = open(outmeans, "w")
    for offset, mean in zip(offsets, means):
      ostr.write("%i %f\n" % (offset, mean))
    ostr.close()

def calc_log_foldchange(indat, refdat, outdat):
  """
  Write a new array containing, at each point, log2(indat/refdat)
  """

  offsets, expvals = numpy.loadtxt(indat, unpack=True)
  dummy, refvals = numpy.loadtxt(refdat, unpack=True)

  newvals = scipy.log2(expvals / refvals)
  write_grfile(offsets, newvals, outdat)


def normalize_vs_ref(indat, indat80, refdatmeanstd, refdat80, outraw, outz, oricorr = True, gaps=False):
  """
  Normalize hyb data following the procedure of Vora2009

  In brief, each data point in the output array (outraw) is:
    exp * (ref80/exp80)
    Where exp is the experimental raw data and ref80/exp80 are the 
      experimental and reference hybs smoothed over 80 kb

  Then, we calculate a z score at each point, z = (expnorm - refmean) / refstd
    Where refmean and refstd are of the *unsmoothed* ref data

  refdatmeanstd is a special case, and must contain a set of rows containing
    offsets, means, and stdevs of the reference sets
    It can be written with write_ref_mean_std

  If oricorr is False, then the 80 kb window normalization will be
    skipped, and z scores calculated just based on the significance
    of the experimental vs. reference data
    In this case refdat80 is not used, and may be None

  If gaps is true, we only use data from the reference hybs that matches offsets from the input data
    By default we assume that they are already equivalent
  """


  if (gaps):

    tmp80 = mktemp()
    tmp1 = mktemp()
    filter_grfile_by_bpranges(indat80, refdat80, tmp80)
    filter_grfile_by_bpranges(indat, refdatmeanstd, tmp1)

    os.rename(tmp80, refdat80)
    os.rename(tmp1, refdatmeanstd)

  offsets, expvals = numpy.loadtxt(indat, unpack=True)
  if (oricorr):
    expvals80 = numpy.loadtxt(indat80, usecols=(1,), unpack=True)
    refvals80 = numpy.loadtxt(refdat80, usecols=(1,), unpack=True)

    # Get and write the normalized experimental data
    expvals_normed = expvals * (refvals80 / expvals80)
    write_grfile(offsets, expvals_normed, outraw)
  else:
    expvals_normed = expvals
    write_grfile(offsets, expvals, outraw)


  # Now calculate and write the z scores
  refmean, refstd = numpy.loadtxt(refdatmeanstd, usecols=(1,2), unpack=True)
  expvals_zscores = (expvals_normed - refmean) / refstd
  write_grfile(offsets, expvals_zscores, outz)

def plot_zscore_hist(zscoresfile, histfigfile, flagfile=False, zerolabel="0", onelabel="1", min=-10, max=10, bins=500):
  """
  Plot a histogram showing histograms of z scores 

  If flagfile exists, it must give at each position either a 0 or 1, in which case
    separate histograms will be constructed for the two flags

  Zerolabel and onelabel are used to label the two histograms in that case

  Flagfile must have the same number of entries as zscoresfile
  """

  import pylab

  indices, data = numpy.loadtxt(zscoresfile, unpack=True)

  if (flagfile):
    #flagfile_stream, flagfile_use = mkstemp()
    #filter_grfile_by_bpranges(zscoresfile, flagfile, flagfile_use)
    flagdata = numpy.loadtxt(flagfile, usecols=(1,), unpack=True)
    zerodata = data[scipy.equal(flagdata, 0)]
    onedata = data[scipy.equal(flagdata, 1)]
    pylab.figure()
    pylab.hist(zerodata, bins=bins, normed=True, alpha=0.5, facecolor='blue', label=zerolabel, range=(min,max), edgecolor='none')
    pylab.hist(onedata, bins=bins, normed=True, alpha=0.5, facecolor='red', label=onelabel, range=(min,max), edgecolor='none')
    pylab.legend()
    pylab.savefig(histfigfile)
    
  else:
    pylab.figure()
    pylab.hist(data, bins=50, normed=True)
    pylab.savefig(histfigfile)

    
def plot_zscore_cum_hist(zscoresfile, histfigfile, flagfile=False, zerolabel="0", onelabel="1", hyptest=True, min=-10, max=10):
  """
  Plot a cumulative histogram showing histograms of z scores 

  If flagfile exists, it must give at each position either a 0 or 1, in which case
    separate histograms will be constructed for the two flags

  Zerolabel and onelabel are used to label the two histograms in that case

  If hyptest is true and two distributions are given, we also do a rank
    sum test to see if their means are equal

  Flagfile must have the same number of entries as zscoresfile
  """

  import pylab

  indices, data = numpy.loadtxt(zscoresfile, unpack=True)

  if (flagfile):
    #flagfile_stream, flagfile_use = mkstemp()
    #filter_grfile_by_bpranges(zscoresfile, flagfile, flagfile_use)
    flagdata = numpy.loadtxt(flagfile, usecols=(1,), unpack=True)
    zerodata = data[scipy.equal(flagdata, 0)]
    onedata = data[scipy.equal(flagdata, 1)]
    pylab.figure()
    pylab.hist(zerodata, bins=5000, normed=True, cumulative=True, histtype='step', label=zerolabel, range=(min,max))
    pylab.hist(onedata, bins=5000, normed=True, cumulative=True, histtype='step', label=onelabel, range=(min,max))
    pylab.legend()
    pylab.axis('tight')
    pylab.savefig(histfigfile)
    print "NOW SAVING TO %s in plot_zscore_cum_hist A" % histfigfile
    #if (hyptest):
    #  print "%10e" % rpy.r.t_test(zerodata, onedata)['p.value']
    
  else:
    pylab.figure()
    pylab.hist(data, bins=500, normed=True, cumulative=True, histtype='step')
    print "NOW SAVING TO %s in plot_zscore_cum_hist B" % histfigfile
    pylab.savefig(histfigfile)

def plot_multiway_cum_hist(zscoresfile, histfigfile, flagfile, keyvals, labels, minval=-10, maxval=10, xlabel=None, ylabel=None, searchkeys=False):
  """
  Plot a cumulative histogram showing histograms of z scores 

  Keyvals is a list of keys from flagfile corresponding to each histogram

  labels should be a similar length list of labels

  If searchkeys is false (the default) we treat the flags in flagfile as integers
  If it is true, we instead do a string search for each key in each flag
    Note that with searchkeys=True, one can plot histograms for nondisjoint sets
  """

  import pylab

  indices, data = numpy.loadtxt(zscoresfile, unpack=True)

  collist = ["black", "black", "red", "red"]
  linelist = ["dashed", "solid", "dashed", "solid"]

  flagdata = numpy.loadtxt(flagfile, usecols=(1,), unpack=True)
  datarrs = []
  print data.shape
  print flagdata.shape
  for key in keyvals:
    if searchkeys:
      goodflags = [mystr.find(key) for mystr in data]
      goodflags = scipy.array(goodflags, type=bool)
      datarrs.append(goodflags)
    else:
      datarrs.append(data[scipy.equal(flagdata, key)])
  print len(datarrs)
  pylab.figure()
  for i in range(len(datarrs)):
    pylab.hist(datarrs[i], bins=5000, normed=True, cumulative=True, histtype='step', range=(minval,maxval), ls=linelist[i], ec=collist[i], fc=collist[i], lw=2)
    pylab.plot( [-11,-10], [-11,-10], ls=linelist[i], c=collist[i], lw=2, label=labels[i])
  pylab.legend(loc=4)
  pylab.xlabel(xlabel)
  pylab.ylabel(ylabel)
  pylab.axis('tight')
  pylab.xlim((minval,maxval))
  pylab.ylim((0,1))
  print "NOW SAVING TO %s in plot_multiway_cum_hist" % histfigfile
  pylab.savefig(histfigfile)

  #return { 'mean' : , 'std': , 'median': }

def make_refs_for_09032010_mg1655(prefixes, refprefix, bpmapfile):
  """
  Make the reference files necessary for do_09032010_analysis_mg1655

  We do the following:
    -extract the pm values (ignore mm)
    -take a 80000 bp running average
    -Fit the 80kbp running average with a spline
    -At each point, scale the data by 2500.0/runavg_spline
    -rescale the data to yield a median value of 2500.0

    We then write on gr file with the means of the 6 references,
      and one with the means and standard deviations
      We also write equivalent files for the log2-transformed data
  """

  GENOME_LENGTH = 4641652

  ref_meanstdfile = "%s_refmeanstd_09032010.txt" % refprefix
  ref_meanfile = "%s_refmeans_09032010.gr" % refprefix
  ref_logmeanstdfile = "%s_reflogmeanstd_09032010.txt" % refprefix
  ref_logmeanfile = "%s_reflogmeans_09032010.gr" % refprefix
  scaled_files = []


  # We need to, for each reference file, go through all the stages up through
  #  normalization of the mean value, and the calculate the mean and std of the six
  #  reference hybs and do an 80 kb  running average over them

  for prefix in prefixes:
    celfile = "%s.CEL" % prefix
    txtfile = "%s_rawval_09032010.txt" % prefix
    cleantxtfile = "%s_cleaned_09032010.txt" % prefix
    ecfile = "%s_Ecoli_probes_09032010.txt" % prefix
    grfile = "%s_ecolionly_09032010.gr" % prefix
    gr80file = "%s_ecolionly_avg80kb_09032010.gr" % prefix
    gr80splinefile = "%s_ecolionly_avg80kb_splinefit_09032010.gr" % prefix
    splinefile_all = "%s_ecolionly_avg80kb_splinefit_alllocs_09032010.gr" % prefix
    ci_file = "%s_ecolionly_copy_indices_09032010.gr" % prefix
    grfile_localnorm = "%s_ecolionly_localavg_09032010.gr" % prefix
    normalized_grfile = "%s_ecolionly_normalized_09032010.gr" % prefix

    # Get a pm gr file containing only E. coli probes
    convert_cel_to_text(celfile, bpmapfile, txtfile)
    clean_text(txtfile, cleantxtfile)
    extract_genome_flag(cleantxtfile, grfile)
    #subtract_mm_from_pm_clamped(ecfile, grfile)

    # Normalization: scale to have median pm of 2500
    do_runningavg_opt(grfile, gr80file, width=80001, genomelength=GENOME_LENGTH)
    spline_smooth_genome(gr80file, gr80splinefile, plot=False, outfile_all = splinefile_all)
    make_ci_file(splinefile_all, ci_file)
    normalize_by_runavg(grfile, gr80splinefile, grfile_localnorm, 2500.0)
    normalize_median_value(grfile_localnorm, normalized_grfile, 2500.0)
    scaled_files.append(normalized_grfile)

  # Now that we have all the files, write the mean and std
  write_ref_mean_std(scaled_files, outfile=ref_meanstdfile, outmeans=ref_meanfile)
  write_ref_log_mean_std(scaled_files, outfile=ref_logmeanstdfile, outmeans=ref_logmeanfile)
    

def make_refs_for_plf_mg1655(prefixes, refprefix, bpmapfile):
  """
  Make the necessary mean/std and 80kb averaged reference files for 11/21/2009 analysis

  """

  GENOME_LENGTH = 4641652

  ref_meanstdfile = "%s_refmeanstd.txt" % refprefix
  ref_meanfile = "%s_refmeans.gr" % refprefix
  ref_80kbfile = "%s_ref80kbavg.gr" % refprefix

  # We need to, for each reference file, go through all the stages up through
  #  normalization of the mean value, and the calculate the mean and std of the six
  #  reference hybs and do an 80 kb  running average over them

  for prefix in prefixes:
    celfile = "%s.CEL" % prefix
    txtfile = "%s_rawvel.txt" % prefix
    cleantxtfile = "%s_allprobes.txt" % prefix
    atfile = "%s_Athaliana_probes.txt" % prefix
    ecfile = "%s_Ecoli_probes.txt" % prefix
    ecfile_bgsub = "%s_Ecoli_probes_bgcorrected.txt" % prefix
    ecfile_corr = "%s_Ecoli_probecorr.gr" % prefix
    ecfile_norm = "%s_Ecoli_normalized.gr" % prefix

    # Obtain a text file that has only the genome name, position, pm, and mm values
    convert_cel_to_text(celfile, bpmapfile, txtfile)
    clean_text(txtfile, cleantxtfile)

    # background correction: subtract median of A. thaliana hybs
    extract_genome_flag_pmmm(cleantxtfile, atfile, genome_re = "chloroplast|chr.")
    extract_genome_flag_pmmm(cleantxtfile, ecfile)
    subtract_percentile_from_pmmm_file(ecfile, atfile, ecfile_bgsub, perc=50)

    # Probe specific correction: Take pm-mm, clamped at zero
    subtract_mm_from_pm_clamped(ecfile_bgsub, ecfile_corr)

    # Normalization: scale to have mean pm-mm of 2500
    normalize_mean_value(ecfile_corr, ecfile_norm)

  # Now that we have all the files, write the mean and std, and then
  #  do an 80 kb running averagr
  write_ref_mean_std(["%s_Ecoli_normalized.gr" % prefix for prefix in prefixes], outfile=ref_meanstdfile, outmeans=ref_meanfile)
  do_runningavg_opt(ref_meanfile, ref_80kbfile, width=80001, genomelength=GENOME_LENGTH)

def smooth_plf_analysis_human(intprefix):
  """
  Do some smoothing of the data from do_plf_analysis_human
  """

  import threading

  CHRLIST = range(1,23) + ["X"]
  #CHRLIST = [1]

  class ipod_slave(threading.Thread):
    def __init__(self, prefix, chr):
      threading.Thread.__init__(self)
      self.prefix = prefix
      self.chr = chr
    def run(self):
      infile = "%s_chr%s_logfoldchange.gr" % (self.prefix, self.chr)
      outfile = "%s_chr%s_logfoldchange_gauss50bp.gr" % (self.prefix,self.chr)
      convolve_file(infile, outfile, 50, 501, wrap=False)

  threads = []
  for chr in CHRLIST:
    print "Starting thread for chromosome %s" % chr
    current = ipod_slave(intprefix, chr)
    threads.append(current)
    current.start()

  for thread in threads:
    thread.join()

def do_plf_analysis_human(intprefix, aqprefix):
  """
  Do analysis on a microarray file as defined in my notebook from 11/28/2009

  We work throgh each chromosome separately
  """


  BPMAPFILE = "/home/petefred/ST_research/ipod/human_genome_data/CD_ENCODE_2/BPMAP/ENCODE_2F_v01-2_NCBIv35.bpmap"

  CHRLIST = range(1,23) + ["X"]
  for chr in CHRLIST:
    print "Working on chromosome %s" % chr
    intfile_tmp = "%s_chr%s.gr" % (intprefix,chr)
    aqfile_tmp = "%s_chr%s.gr" % (aqprefix, chr)
    outfile = "%s_chr%s_logfoldchange.gr" % (intprefix, chr)

    extract_genome_flag(intfile, intfile_tmp, "chr%s" % chr)
    extract_genome_flag(aqfile, aqfile_tmp, "chr%s" % chr)

    try:
      sort_gr_file(intfile_tmp, intfile_tmp)
    except:
      print "  Failed on current chromosome... skipping"
      continue
    sort_gr_file(aqfile_tmp, aqfile_tmp)

    calc_log_foldchange(intfile_tmp, aqfile_tmp, outfile)

  # first, convert the cel file to text and do normalization
  
  for prefix in [intprefix, aqprefix]:
    celfile = "%s.CEL" % prefix 
    txtfile = "%s_rawvel.txt" % prefix

    atfile = "%s_Athaliana_probes.txt" % prefix
    humanfile = "%s_humanprobes.txt" % prefix
    atfile_clean = "%s_Athaliana_probes_clean.txt" % prefix
    humanfile_clean = "%s_humanprobes_clean.txt" % prefix
    humanfile_bgsub = "%s_humanprobes_bgcorrected.txt" % prefix
    humanfile_corr = "%s_humanprobes_probecorr.txt" % prefix
    humanfile_norm = "%s_humanprobes_norm.txt" % prefix

    convert_cel_to_text(celfile, BPMAPFILE, txtfile)

    extract_seqgroup(txtfile, atfile, "At")
    extract_seqgroup(txtfile, humanfile, "Hs")
    clean_text(atfile, atfile_clean)
    clean_text(humanfile, humanfile_clean)
    subtract_percentile_from_genomefile(humanfile_clean, atfile_clean, humanfile_bgsub)
    subtract_mm_from_pm_withgenome(humanfile_bgsub, humanfile_corr)
    normalize_mean_value_withgenome(humanfile_corr, humanfile_norm)

  # now, generate summary scores for each chromosome
  intfile = "%s_humanprobes_norm.txt" % intprefix
  aqfile = "%s_humanprobes_norm.txt" % aqprefix

  for chr in CHRLIST:
    print "Working on chromosome %s" % chr
    intfile_tmp = "%s_chr%s.gr" % (intprefix,chr)
    aqfile_tmp = "%s_chr%s.gr" % (aqprefix, chr)
    outfile = "%s_chr%s_logfoldchange.gr" % (intprefix, chr)

    extract_genome_flag(intfile, intfile_tmp, "chr%s" % chr)
    extract_genome_flag(aqfile, aqfile_tmp, "chr%s" % chr)

    try:
      sort_gr_file(intfile_tmp, intfile_tmp)
    except:
      print "  Failed on current chromosome... skipping"
      continue
    sort_gr_file(aqfile_tmp, aqfile_tmp)

    calc_log_foldchange(intfile_tmp, aqfile_tmp, outfile)

def calc_normed_score_deltas(prefix, pathdict):
  """
  Calculate the difference in normed scores between a given prefix and both tjv's and plf's genomic reference hybs

  We assume that get_normed_corrected_intensities has already been run for both prefix and the reference hybs

  Each position is then assigned a scaled score based on (val - mean(refvals)) / 2500.0

  Note that these are *not* z-scores, and have nothing to do with the uncertainty at a particular position. We're just taking an average
    of the references and then dividing by the median score of both arrays
  """

  ingrfile = "%s_Ecoli_normed_and_corrected.gr" % prefix

  reffiles = pathdict['REFGRFILES']
  refnames = pathdict['REFNAMES']

  origpos,origval = read_grfile(ingrfile)

  for i in range(len(reffiles)):
    reffilelist = reffiles[i]
    refname = refnames[i]
    outgrfile = "%s_delta_from_%s_refs.gr" % (prefix, refname)
    gaussedgrfile = "%s_delta_from_%s_refs_gauss4bp.gr" % (prefix, refname)


    refhybs_list=[]

    for hyb in reffilelist:
      pos,val = read_grfile(hyb)
      refhybs_list.append(val)

    refhybs_arr = scipy.array(refhybs_list)
    refhybs_list = []

    delta_vals = (origval - scipy.mean(refhybs_arr, axis=0)) / 2500.0

    write_grfile(origpos, delta_vals, outgrfile)
    convolve_file(outgrfile, gaussedgrfile)

def get_normed_corrected_intensities_pmonly(prefix, pathdict):
  """
  Obtain PM intensities with background correction and normalization to a median value of 2500.0 by scaling

  This analysis assumes that we're working on mg1655

  Because we ignore the mm values, this should only be used when
    comparing to a control hyb
  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']

  celfile = "%s.CEL" % prefix
  txtfile = "%s_rawvel.txt" % prefix
  cleantxtfile = "%s_allprobes.txt" % prefix

  atfile = "%s_Athaliana_probes.txt" % prefix
  ecfile = "%s_Ecoli_probes.txt" % prefix
  ecfile_bgsub = "%s_Ecoli_probes_bgcorrected.txt" % prefix
  ecfile_bgsub_nomm = "%s_Ecoli_probes_bgcorrected_pmonly.txt" % prefix
  ecfile_normalized = "%s_Ecoli_normed_and_corrected_pmonly.gr" % prefix

  # Obtain a text file that has only the genome name, position, pm, and mm values
  #convert_cel_to_text(celfile, BPMAPFILE, txtfile)
  #clean_text(txtfile, cleantxtfile)

  # background correction: subtract median of A. thaliana hybs
  #extract_genome_flag_pmmm(cleantxtfile, atfile, genome_re = "chloroplast|chr.")
  #extract_genome_flag_pmmm(cleantxtfile, ecfile)
  #subtract_percentile_from_pmmm_file(ecfile, atfile, ecfile_bgsub, perc=50)

  remove_mm_probes(ecfile_bgsub, ecfile_bgsub_nomm)

  # Normalization of median corrected intensity
  normalize_median_value(ecfile_bgsub_nomm, ecfile_normalized)


def get_normed_corrected_intensities(prefix, pathdict):
  """
  Obtain PM-MM intensities with background correction and normalization to a median value of 2500.0 by scaling

  This analysis assumes that we're working on mg1655
  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']

  celfile = "%s.CEL" % prefix
  txtfile = "%s_rawvel.txt" % prefix
  cleantxtfile = "%s_allprobes.txt" % prefix

  atfile = "%s_Athaliana_probes.txt" % prefix
  ecfile = "%s_Ecoli_probes.txt" % prefix
  ecfile_bgsub = "%s_Ecoli_probes_bgcorrected.txt" % prefix
  ecfile_corr = "%s_Ecoli_probecorr.gr" % prefix
  ecfile_normalized = "%s_Ecoli_normed_and_corrected.gr" % prefix

  # Obtain a text file that has only the genome name, position, pm, and mm values
  convert_cel_to_text(celfile, BPMAPFILE, txtfile)
  clean_text(txtfile, cleantxtfile)

  # background correction: subtract median of A. thaliana hybs
  extract_genome_flag_pmmm(cleantxtfile, atfile, genome_re = "chloroplast|chr.")
  extract_genome_flag_pmmm(cleantxtfile, ecfile)
  subtract_percentile_from_pmmm_file(ecfile, atfile, ecfile_bgsub, perc=50)

  # Probe specific correction: Take pm-mm, clamped at zero
  subtract_mm_from_pm_clamped(ecfile_bgsub, ecfile_corr)

  # Normalization of median corrected intensity
  normalize_median_value(ecfile_corr, ecfile_normalized)

def do_plf_analysis_mg1655(prefix, pathdict):
  """
  Do analysis on a microarray file as defined in my notebook from 11/21/2009

  This analysis assumes that we're working on mg1655
  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']

  celfile = "%s.CEL" % prefix
  txtfile = "%s_rawvel.txt" % prefix
  cleantxtfile = "%s_allprobes.txt" % prefix

  atfile = "%s_Athaliana_probes.txt" % prefix
  ecfile = "%s_Ecoli_probes.txt" % prefix
  ecfile_bgsub = "%s_Ecoli_probes_bgcorrected.txt" % prefix
  ecfile_corr = "%s_Ecoli_probecorr.gr" % prefix
  ecfile_80kbavg = "%s_Ecoli_80kbavg.gr" % prefix
  ecscores_raw = "%s_plfanalysis_fullnorm.gr" % prefix
  ecscores_zscore = "%s_plfanalysis_zscores.gr" % prefix
  

  print "Performing PLF (11/21/2009) analysis on cel file %s" % celfile

  # Obtain a text file that has only the genome name, position, pm, and mm values
  convert_cel_to_text(celfile, BPMAPFILE, txtfile)
  clean_text(txtfile, cleantxtfile)

  # background correction: subtract median of A. thaliana hybs
  extract_genome_flag_pmmm(cleantxtfile, atfile, genome_re = "chloroplast|chr.")
  extract_genome_flag_pmmm(cleantxtfile, ecfile)
  subtract_percentile_from_pmmm_file(ecfile, atfile, ecfile_bgsub, perc=2)

  # Probe specific correction: Take pm-mm, clamped at zero
  subtract_mm_from_pm_clamped(ecfile_bgsub, ecfile_corr)

  # Perform replication-dependent normalization and produce summary statistic
  do_runningavg_opt(ecfile_corr, ecfile_80kbavg, width=80001, genomelength=GENOME_LENGTH)
  normalize_vs_ref(ecfile_corr, ecfile_80kbavg, REFS_MEANSTD, REFS_80KB, ecscores_raw, ecscores_zscore)

def do_plf_analysis_mg1655_splitstrand(prefix, pathdict):
  """
  Do split-strand analysis on a microarray file

  This analysis assumes that we're working on mg1655

  We do the same analysis as in do_plf_analysis_mg1655, but on the t and f
    strands separately

  This function assumes that do_plf_analysis_mg1655 has already been run; it just
    produces separate z score files

  We also produce smoothed data that is ready for viewing/analysis

  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']
  forward_flags = pathdict['FORWARD_FLAG_FILE']
  reverse_flags = pathdict['REVERSE_FLAG_FILE']

  ecscores_zscore = "%s_plfanalysis_zscores.gr" % prefix

  forward_zscores = "%s_plfanalysis_zscores_plusstrand.gr" % prefix
  reverse_zscores = "%s_plfanalysis_zscores_minusstrand.gr" % prefix

  forward_zscore_filtered = "%s_plfanalysis_zscores_plusstrand_filtered.gr" % prefix
  forward_gaussed = "%s_plfanalysis_zscores_plusstrand_gaussed4bp.gr" % prefix

  reverse_zscore_filtered = "%s_plfanalysis_zscores_minusstrand_filtered.gr" % prefix
  reverse_gaussed = "%s_plfanalysis_zscores_minusstrand_gaussed4bp.gr" % prefix

  filter_grfile_by_bpranges(forward_flags, ecscores_zscore, forward_zscores)
  filter_grfile_by_bpranges(reverse_flags, ecscores_zscore, reverse_zscores)

  remove_nans(forward_zscores, forward_zscore_filtered)
  convolve_file(forward_zscore_filtered, forward_gaussed, genomelength=GENOME_LENGTH, gausswidth=8, width=97)
  remove_nans(reverse_zscores, reverse_zscore_filtered)
  convolve_file(reverse_zscore_filtered, reverse_gaussed, genomelength=GENOME_LENGTH, gausswidth=8, width=97)


def prep_for_viewing_mg1655(prefix, pathdict):
  """
  Filter any inf/nan values from a zscore file and gauss it
  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']

  ecscores_zscore = "%s_plfanalysis_zscores.gr" % prefix
  ecscores_zscore_filtered = "%s_plfanalysis_zscores_filtered.gr" % prefix
  ecscores_gaussed = "%s_plfanalysis_zscores_gaussed4bp.gr" % prefix

  
  remove_nans(ecscores_zscore, ecscores_zscore_filtered)
  convolve_file(ecscores_zscore_filtered, ecscores_gaussed, genomelength=GENOME_LENGTH)

def prep_for_viewing_tjvscores(prefix, pathdict):
  """
  Filter any inf/nan values from a zscore file and gauss it
  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']

  ecscores_zscore = "%s_standard-analysis_ecolionly_zscores.gr" % prefix
  ecscores_zscore_filtered = "%s_standard-analysis_zscores_filtered.gr" % prefix
  ecscores_gaussed = "%s_standard-analysis_gaussed_zscores_4bp.gr" % prefix

  
  remove_nans(ecscores_zscore, ecscores_zscore_filtered)
  convolve_file(ecscores_zscore_filtered, ecscores_gaussed, genomelength=GENOME_LENGTH)

def do_09032010_analysis_mg1655(prefix, pathdict):
  """
  Run through analysis of ipod data on the file prefix.CEL

  We do the following:
    -extract the pm values (ignore the mm)
    -take a 80000 bp running average
    -Smooth the 80kbp running average with a spline
    -At each point, scale the data by 2500.0/runavg_spline
    -rescale the data to yield a median value of 2500.0
    -score the probe as log2(val) - mean(log2(refs)) / sigma(log2(refs))
    -z-score the probe as val - mean(refs) / sigma(refs)
    -also calculate a t-score as defined in 9/10/2010 notes
    -change the probe indices by CENTER_OFFSET to correct for the fact
      that the affy probes are labeled by their side, not center

      Here refs are Tiffany's 6 reference hybs, prepared by
        make_refs_09032010

  This analysis assumes that we're working on E. coli
  """

  # Distance by which we need to offset the probe index to reach the 
  #   probe center
  CENTER_OFFSET = 12

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_MEANSTD_LOG = pathdict['REFS_MEANSTD_LOG']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']

  celfile = "%s.CEL" % prefix
  txtfile = "%s_09032010.txt" % prefix
  cleantxtfile = "%s_cleaned_09032010.txt" % prefix
  ecfile = "%s_Ecoli_probes_09032010.txt" % prefix
  grfile = "%s_ecolionly_09032010.gr" % prefix
  gr80file = "%s_ecolionly_avg80kb_09032010.gr" % prefix
  gr80splinefile = "%s_ecolionly_avg80kb_splinefit_09032010.gr" % prefix
  splinefile_all = "%s_ecolionly_avg80kb_splinefit_alllocs_09032010.gr" % prefix
  grfile_localnorm = "%s_ecolionly_localavg_09032010.gr" % prefix
  ci_file = "%s_ecolionly_copy_index_09032010.gr" % prefix
  normdat_raw = "%s_ecolionly_normalized_09032010.gr" % prefix
  normdat_zscore = "%s_standard-analysis_ecolionly_zscores_09032010.gr" % prefix
  normdat_zscore_offset = "%s_standard-analysis_ecolionly_zscores_offset_09032010.gr" % prefix
  normdat_zscore_log = "%s_standard-analysis_ecolionly_zscores_logtransformed_09032010.gr" % prefix
  normdat_zscore_offset_log = "%s_standard-analysis_ecolionly_zscores_logtransformed_offset_09032010.gr" % prefix
  normdat_tscore = "%s_standard-analysis_ecolionly_tscores_09032010.gr" % prefix
  normdat_tscore_offset = "%s_standard-analysis_ecolionly_tscores_offset_09032010.gr" % prefix
  normdat_tscore_log = "%s_standard-analysis_ecolionly_tscores_logtransformed_09032010.gr" % prefix
  normdat_tscore_offset_log = "%s_standard-analysis_ecolionly_tscores_logtransformed_offset_09032010.gr" % prefix

  convert_cel_to_text(celfile, BPMAPFILE, txtfile)
  clean_text(txtfile, cleantxtfile)
  extract_genome_flag(cleantxtfile, grfile)
  do_runningavg_opt(grfile, gr80file, width=80001, genomelength=GENOME_LENGTH)
  spline_smooth_genome(gr80file, gr80splinefile, plot=False, outfile_all =  splinefile_all)
  make_ci_file(splinefile_all, ci_file)
  normalize_by_runavg(grfile, gr80splinefile, grfile_localnorm, 2500.0)
  normalize_median_value(grfile_localnorm, normdat_raw, 2500.0)
  score_normed_vals_nolog(normdat_raw, REFS_MEANSTD,  normdat_zscore)
  apply_index_offset(normdat_zscore, normdat_zscore_offset, CENTER_OFFSET)
  score_normed_vals(normdat_raw, REFS_MEANSTD_LOG,  normdat_zscore_log)
  apply_index_offset(normdat_zscore_log, normdat_zscore_offset_log, CENTER_OFFSET)


def do_standard_analysis_mg1655(prefix, pathdict):
  """
  Run through the standard set of analysis for ipod data on the file prefix.CEL

  This follows, as closely as possible, the procedure used in tvora's scripts

  This analysis assumes that we're working on E. coli
  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']

  celfile = "%s.CEL" % prefix
  txtfile = "%s.txt" % prefix
  cleantxtfile = "%s_cleaned.txt" % prefix
  grfile = "%s_ecolionly.gr" % prefix
  grfile_sc2500 = "%s_ecolionly_sc2500.gr" % prefix
  gr80file = "%s_ecolionly_sc2500_avg80kb.gr" % prefix
  normdat_raw = "%s_ecolionly_normalized.gr" % prefix
  normdat_zscore = "%s_standard-analysis_ecolionly_zscores.gr" % prefix
  hist_fig_file = "%s_standard-analysis_coding-noncoding_vs_zscore.pdf" % prefix
  chist_fig_file = "%s_standard-analysis_coding-noncoding_vs_zscore_chist.pdf" % prefix
  smoothed_zscore_file = "%s_standard-analysis_smoothed_zscores_48bp.gr" % prefix
  gaussed_zscore_file = "%s_standard-analysis_gaussed_zscores_4bp.gr" % prefix
  smoothed_prefix = "%s_standard-analysis_smoothed" % prefix
  gaussed_prefix = "%s_standard-analysis_gaussed" % prefix

  convert_cel_to_text(celfile, BPMAPFILE, txtfile)
  clean_and_correct_text(txtfile, cleantxtfile)
  extract_genome_flag(cleantxtfile, grfile)
  offset_median_value(grfile, grfile_sc2500, 2500.0)
  do_runningavg_opt(grfile_sc2500, gr80file, width=80001, genomelength=GENOME_LENGTH)
  normalize_vs_ref(grfile_sc2500, gr80file, REFS_MEANSTD, REFS_80KB, normdat_raw, normdat_zscore)
  do_runningavg_opt(normdat_zscore, "%s_runavg128.gr" % smoothed_prefix, width=129, genomelength=GENOME_LENGTH)
  plot_zscore_hist("%s_runavg128.gr" % smoothed_prefix, hist_fig_file, CODINGFLAGFILE, zerolabel="Coding", onelabel="Noncoding")
  plot_zscore_cum_hist("%s_runavg128.gr" % smoothed_prefix, chist_fig_file, CODINGFLAGFILE, zerolabel="Coding", onelabel="Noncoding")
  convolve_file(normdat_zscore, gaussed_zscore_file, genomelength=GENOME_LENGTH)

def do_standard_analysis_mg1655_splitstrand(prefix, pathdict):
  """
  Run through the standard set of analysis for ipod data on the file prefix.CEL

  We produce the same results as do_standard_analysis_mg1655, but produce
    separate z score files for each strand

  We assume do_standard_analysis_mg1655 has already been run

  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']
  forward_flags = pathdict['FORWARD_FLAG_FILE']
  reverse_flags = pathdict['REVERSE_FLAG_FILE']

  normdat_zscore = "%s_standard-analysis_ecolionly_zscores.gr" % prefix

  forward_zscores = "%s_standard-analysis_zscores_plusstrand.gr" % prefix
  reverse_zscores = "%s_standard-analysis_zscores_minusstrand.gr" % prefix

  forward_zscore_filtered = "%s_standard-analysis_zscores_plusstrand_filtered.gr" % prefix
  forward_gaussed = "%s_standard-analysis_zscores_plusstrand_gaussed4bp.gr" % prefix

  reverse_zscore_filtered = "%s_standard-analysis_zscores_minusstrand_filtered.gr" % prefix
  reverse_gaussed = "%s_standard-analysis_zscores_minusstrand_gaussed4bp.gr" % prefix
  
  filter_grfile_by_bpranges(forward_flags, normdat_zscore, forward_zscores)
  filter_grfile_by_bpranges(reverse_flags, normdat_zscore, reverse_zscores)

  remove_nans(forward_zscores, forward_zscore_filtered)
  convolve_file(forward_zscore_filtered, forward_gaussed, genomelength=GENOME_LENGTH, gausswidth=8, width=97)
  remove_nans(reverse_zscores, reverse_zscore_filtered)
  convolve_file(reverse_zscore_filtered, reverse_gaussed, genomelength=GENOME_LENGTH, gausswidth=8, width=97)


def analyze_files(prefixes, pathdict, analysis_function):
  """
  Analyze a set of files with given prefixes in parallel

  This runs the specified analysis function on each prefix

  pathdict must contain a set of paths labeled
    REFS_MEANSTD (mean and stdev of the reference hybs)
    REFS_80KB (running averages of the reference hyb means over 80 kb)
    BPMAPFILE (path to the relevant .bpmap file)
    CODINGFLAGFILE (path containing flags for whether or not each probe is
      in a coding sequence)
  """

  import threading 

  class ipod_slave(threading.Thread):
    def __init__(self, prefix, pathdict):
      threading.Thread.__init__(self)
      self.prefix = prefix
      self.pathdict = pathdict
    def run(self):
      analysis_function(self.prefix, self.pathdict)

  threads = []
  for prefix in prefixes:
    print "Starting thread for prefix %s" % prefix
    current = ipod_slave(prefix, pathdict)
    threads.append(current)
    current.start()

  for thread in threads:
    thread.join()

def calc_corrcoef_matrix(filenames, displaynames, printmat=True, clampvals = None):
  """
  Calculate the correlation coefficients between the given files

  The matrix will be returned if printmat is false, and printed if
    printmat is true
  """

  from sys import stdout

  corrs = scipy.ones((len(filenames), len(filenames)))

  for i in range(len(filenames)):
    file1 = filenames[i]
    for j in range(i):
      file2 = filenames[j]
      corrs[i,j] = calc_corrcoef(file1, file2, clampvals)
      corrs[j,i] = corrs[i,j]

  if (printmat):
    stdout.write("           |")
    for name in displaynames:
      stdout.write(" %10s |" % name) 
    stdout.write("\n")
    for i in range(len(filenames)):
      file1 = filenames[i]
      stdout.write("%10s |" % displaynames[i])
      for j in range(i):
        stdout.write("    %5f    |" % corrs[i,j])
      print ""
  else:
    return corrs


def calc_corrcoef(file1, file2, clampvals):
  """
  Calculate the correlation coefficient for ipod data in file1 and file2
  """

  data1 = numpy.loadtxt(file1, usecols=(1,), unpack=True)
  data2 = numpy.loadtxt(file2, usecols=(1,), unpack=True)

  if (clampvals):
    print "Clamping data to be between %f and %f" % (minval, maxval)
    minval, maxval = clampvals
    data1[scipy.less(data1, minval)] = minval
    data2[scipy.less(data2, minval)] = minval
    data1[scipy.greater(data1, maxval)] = maxval
    data2[scipy.greater(data2, maxval)] = maxval
  cc = scipy.stats.spearmanr(data1, data2)
  print "Correlation between %s and %s is %s" % (file1, file2, cc)
  return cc[0]

def generate_averages(indata, avgprefix, windowsizes = []):
  """
  For each specified window size, generate running averages of that width
  """

  for window in windowsizes:
    convolve_file(indata, "%s_runavg%s.gr" % (avgprefix, window), window+1, kerneltype="gauss")
    convolve_file(indata, "%s_gauss%s.gr" % (avgprefix, window), window/10.0, window+1)

def circular_range(datvec, startindex, endindex):
  """
  Return the range of the data array corresponding to the start and end indices

  Unlike normal array indexing, we transparently assume a circular array
    and retrieve the corresponding data points

  Also, unlike normal python array slicing, the last element **is included**

  This is done to make working with ranges centered on a specific value
    more transparent, as this occurs frequently in ipod_utils
  """

  arrlen = len(datvec)
  
  startvec = []
  endvec = []

  if (startindex < 0):
    beforeindex = arrlen + startindex
    startindex = 0
    startvec = datvec[beforeindex:]
  
  if (endindex >= arrlen):
    afterindex = endindex - arrlen
    endindex = len(datvec) - 1
    endvec = datvec[:(afterindex + 1)]

  bodyslice = datvec[startindex:(endindex + 1)]

  if (len(startvec)):
    bodyslice = scipy.hstack((startvec, bodyslice))
  if (len(endvec)):
    bodyslice = scipy.hstack((bodyslice, endvec))

  return bodyslice

def circular_range_bps(offsetvec_orig, datvec_orig, startbp, endbp, genomelength=None, guess_start=None, guess_end = None, return_inds = False):

  """
  Similar to circular_range, but using bp units

  This means that the number of elements returned is uncertain, since it
    depends on the spacing between entries

  We return the set of all entries completely contained with startbp
    and endbp, and the corresponding slice of offsetvec

  If genomelength is given, it is used as the total length of the genome
  Otherwise, we estimate it based on the last two entries

  guess_start and guess_end are guesses at where we should start looking,
    to avoid having to search through the entire array

  If return_inds is true, we return the indices passed to circular_range
    as well as the usual values
  """

  import bisect


  offsetvec = numpy.array(offsetvec_orig, copy=True)
  datvec = datvec_orig

  startbp=int(startbp)
  endbp=int(endbp)

  arrlen = len(datvec)
  if genomelength:
    lastbp = genomelength ;# last coordinate in the genome
  else:
    lastbp = offsetvec[-1] + (offsetvec[-1] - offsetvec[-2])

  if (startbp > endbp):
    tmpbp = endbp
    endbp=startbp
    startbp=tmpbp

    guesstmp = guess_start
    guess_start = guess_end
    guess_end = guesstmp

  if (startbp < 1):
    startbp += lastbp

  if endbp > offsetvec[-1]:
    endbp -= lastbp

  # limit the bounds of the search if we are able to
  searchleft_start = 0
  searchleft_end = len(offsetvec)
  searchright_start = 0
  searchright_end = len(offsetvec)
  width=2*abs(startbp-endbp)


  if (guess_start is not None) and (guess_start%arrlen > width) and (guess_start%arrlen < (arrlen-width)):
    searchleft_start = guess_start - width
    searchleft_end = guess_start + width
    #print "Using guessed start"
  if (guess_end is not None) and (guess_end%arrlen > width) and (guess_end%arrlen < (arrlen-width)):
    searchright_start = guess_end - width
    searchright_end = guess_end + width

  #print "Search positions are %i--%i %i--%i" % (searchleft_start, searchleft_end, searchright_start, searchright_end)

  startindex = bisect.bisect_left(offsetvec, startbp, searchleft_start, searchleft_end)
  endindex = bisect.bisect_right(offsetvec, endbp, searchright_start, searchright_end)

  if (startindex > endindex):
    offset_slice = numpy.concatenate((offsetvec[startindex:],offsetvec[:endindex]))
    data_slice = numpy.concatenate((datvec[startindex:],datvec[:endindex]))
  else:
    offset_slice = offsetvec[startindex:endindex]
    data_slice = datvec[startindex:endindex]

  #offset_slice = circular_range(offsetvec, startindex, endindex)
  #data_slice = circular_range(datvec, startindex, endindex)

  #print (startindex, endindex)

  if (return_inds):
    return (offset_slice, data_slice, startindex, endindex)
  else:
    return (offset_slice, data_slice)
#
#  
#  # handle if we're off the edge on the left
#  if (startbp < 1):
#    beforebp = lastbp + startbp ;# real, wrapped index of position to begin at 
#    print "*"
#    print lastbp
#    print startbp
#    print beforebp
#    print offsetvec[-1]
#    if (guess_start and beforebp < offsetvec[-1]):
#      startindex = guess_start
#    else:
#      startindex = -1
#  # now if we're off the edge on the right
#  elif (startbp > offsetvec[-1]):
#    if (guess_start):
#      startindex = guess_start
#    else:
#      startindex = len(offsetvec)
#
#
#  # now the normal case
#  else:
#    if (guess_start):
#      startindex = guess_start
#    else:
#      startindex = 0
#
#  #now look for the correct place for startindex
#  if ((startindex%arrlen) < (len(offsetvec) / 2)):
#
#    dist_right = math.abs(offsetvec[startindex%arrlen] - startbp)
#    dist_left = math.abs(offsetvec[startindex%arrlen] - (startbp - lastbp))
#  else:
#    dist_right = math.abs(offsetvec[startindex%arrlen] - (startbp+lastbp))
#    dist_left = math.abs(offsetvec[startindex%arrlen] - startbp)
#
#
#
#
#
#  while (offsetvec[startindex] < startbp):
#    startindex += 1
#
#  if (endbp > lastbp):
#    afterbp = endbp - lastbp
#    if (guess_end):
#      if (guess_end >= len(offsetvec)):
#        guess_end -= len(offsetvec)
#      endindex = guess_end
#      while (((offsetvec[endindex]) > afterbp) and ( endindex >= 0)):
#        endindex -= 1
#    else:
#      endindex = 0
#    while (offsetvec[endindex + 1] < afterbp):
#      endindex += 1
#    endindex += arrlen
#  elif (endbp > offsetvec[-1]):
#    if (endbp == lastbp):
#      endindex = len(offsetvec)
#    else:
#      endindex = len(offsetvec) - 1
#  else:
#    if (guess_end):
#      endindex = guess_end
#      while (offsetvec[endindex] > endbp):
#        endindex -= 1
#    else:
#      endindex = startindex
#    while (offsetvec[endindex + 1] < endbp):
#      endindex += 1
#
#  offset_slice = circular_range(offsetvec, startindex, endindex)
#  data_slice = circular_range(datvec, startindex, endindex)
#
#  print (startindex, endindex)
#
#  if (return_inds):
#    return (offset_slice, data_slice, startindex, endindex)
#  else:
#    return (offset_slice, data_slice)

def identify_epods(epod_data, percentile_data, min_epod_length, epod_outfile, lpod=False, delta=25):
  """
  Look for epods over (100 - delta-th) percentile from percentile_data in epod_data

  Each data input should be a file (NOT a vector)
  We then look for all contiguous regions of length at least min_epod_length  
  A file (epod_outfile) containing the locations of epods will also be written

  If lpod is true, then we instead look for similar regions under 25th percentile occupancy

  All units are IN BASE PAIRS
  """

  
  import scipy.stats

  offsets, epod_vec = numpy.loadtxt(epod_data, unpack=True)
  percentile_vec = numpy.loadtxt(percentile_data, usecols=(1,), unpack=True)
  epod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, 100-delta)
  lpod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, delta)
  percentile_vec = []

  # establish a guess for how many bp between each entry
  stride = offsets[1] - offsets[0]

  def epod_cmp(value):
    if (lpod):
      return (value < lpod_cutoff)
    else:
      return (value > epod_cutoff)

  if (lpod):
    print "Searching for regions at least %i bp long with a median z score below %f" % (min_epod_length, lpod_cutoff)
  else:
    print "Searching for regions at least %i bp long with a median z score above %f" % (min_epod_length, epod_cutoff)

  # Now we go through the raw data and find all 1024 bp windows where the
  #   median value is at least percentile_vec
  # Each time we find such a window, it is expanded until we can no longer
  #   keep a valid EPOD by expanding 1 bp at a time
  offset = int(min_epod_length / 2)
  epod_locs = []
  for i in range(int(offsets[-1])):
    if (len(epod_locs) > 0) and i < (epod_locs[-1])[1]:
      continue
#    print "Searching for epod starting at %i" % i
    start = i-offset
    end = i+offset
    curr_median = scipy.median(circular_range_bps(offsets,epod_vec, start, end)[1])
    if (epod_cmp(curr_median)):
      print "Found an epod with median %f between %i and %i" % (curr_median, start, end)
      new_median = curr_median
      newstart = start
      newend = end
      while (new_median > epod_cutoff):
        if len(epod_locs) > 0:
          if start < (epod_locs[-1][1] - 1):
            break

        start=newstart
        newstart -= 1
#        print newstart
        new_median = scipy.median(circular_range_bps(offsets,epod_vec, newstart, end)[1])
#        print new_median

      while (epod_cmp(new_median)):
        end=newend
        newend += 1
#        print newend
        new_median = scipy.median(circular_range_bps(offsets,epod_vec, start, newend)[1])
#        print new_median

      print "After expansion, ipod between %i and %i has median %f" % (start, end, curr_median)
      epod_locs.append( (int(start), int(end)) )
  
  # Take one shot at merging adjacent ipods
  i = 0
  while (i < (len(epod_locs) - 2)):
    j = i+1
    start1, end1 = epod_locs[i]
    start2, end2 = epod_locs[j]
    if (end1 > start2):
      print "Trying to merge epods %s and %s" % (epod_locs[i], epod_locs[j])
      print scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])
      if epod_cmp(scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])):
        print "Merging epods %s and %s" % (epod_locs[i], epod_locs[j])
        epod_locs.pop(j)
        epod_locs[i] = (start1,end2)
        continue

    i += 1

  print epod_locs

  # write a file containing 1s at the positions involved in epods
  epod_loc_vec = scipy.zeros(len(epod_vec))
  for i in range(len(epod_vec)):
    iloc = offsets[i]
    for (start,end) in epod_locs:
      if (start > end):
        if (end < 0):
          end += len(epod_vec)
        else:
          start -= len(epod_vec)
      if (iloc>=start and iloc<= end):
        epod_loc_vec[i] = 1
        continue

  write_grfile(offsets, epod_loc_vec, epod_outfile)


  return epod_locs

def identify_epods_v2(epod_data, percentile_data, min_epod_length, epod_outfile, lpod=False, delta=25):
  """
  Look for epods over (100 - delta-th) percentile from percentile_data in epod_data

  Each data input should be a file (NOT a vector)
  We then look for all contiguous regions of length at least min_epod_length  
  A file (epod_outfile) containing the locations of epods will also be written

  If lpod is true, then we instead look for similar regions under 25th percentile occupancy

  All units are IN BASE PAIRS

  This version of the function tries to expand each epod symmetrically starting from local maxima , which is less likely
   to yield asymmetric epods than is the function above

  """

  
  import scipy.stats

  offsets, epod_vec = numpy.loadtxt(epod_data, unpack=True)
  percentile_vec = numpy.loadtxt(percentile_data, usecols=(1,), unpack=True)
  epod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, 100-delta)
  lpod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, delta)
  percentile_vec = []

  # establish a guess for how many bp between each entry
  stride = offsets[1] - offsets[0]

  def epod_cmp(value):
    if (lpod):
      return (value < lpod_cutoff)
    else:
      return (value > epod_cutoff)

  if (lpod):
    print "Searching for regions at least %i bp long with a median z score below %f" % (min_epod_length, lpod_cutoff)
  else:
    print "Searching for regions at least %i bp long with a median z score above %f" % (min_epod_length, epod_cutoff)

  epod_pot_arr = numpy.zeros_like(epod_vec) - 1000 ;# contains the values of all 1024 bp windows that are potential epods
  # we follow a two-pass approach to find epods
  # first we go through the raw data and find all 1024 bp windows where the
  #   median value is at least percentile_vec
  # Each time we find such a window, we add to epod_pot_arr the score at that location
  #  this way we know the relative heights of various windows
  #  note that we add the plain score, and not the window median, to minimize ties

  offset = int(min_epod_length / 2)
  for i in range(len(offsets)):
    start = offsets[i]-offset
    end = offsets[i]+offset
    curr_median = scipy.median(circular_range_bps(offsets,epod_vec, start, end)[1])
    if (epod_cmp(curr_median)):
      epod_pot_arr[i] = epod_vec[i]

  # now, we find the BEST locations to start potential epods, and try expanding them in either direction

  #numpy.save('test_centers.npy',epod_pot_arr)
  epod_pot_centers = scipy.signal.argrelextrema(epod_pot_arr, numpy.greater_equal,mode='wrap')[0]
  epod_abovezero = numpy.argwhere(epod_pot_arr > -100)

  epod_centers = numpy.intersect1d(epod_pot_centers, epod_abovezero)

  epod_locs = []

  #print "Potential epod centers found at: %s" % epod_centers


  for center_i in epod_centers:
      epod_start = offsets[center_i] - offset
      epod_end = offsets[center_i] + offset

      expand_left = True
      expand_right = True

      while (expand_left or expand_right):
          # we try expanding this epod as far as we can

          if expand_left:
              # try expanding to the left
              trial_median = scipy.median(circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end)[1])
              if trial_median > epod_cutoff:
                  epod_start -= 1
              else:
                  expand_left = False

          if expand_right:
              # try expanding to the right
              trial_median = scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end+1)[1])
              if trial_median > epod_cutoff:
                  epod_end += 1
              else:
                  expand_right = False
      
      print "After expansion, ipod between %i and %i has median %f" % (epod_start, epod_end, scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end)[1]) )
      epod_locs.append( (int(epod_start), int(epod_end)) )
  
  # Take one shot at merging adjacent ipods
  i = 0
  while (i < (len(epod_locs) - 2)):
    j = i+1
    start1, end1 = epod_locs[i]
    start2, end2 = epod_locs[j]
    if (end1 > start2):
      print "Trying to merge epods %s and %s" % (epod_locs[i], epod_locs[j])
      print scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])
      if epod_cmp(scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])):
        print "Merging epods %s and %s" % (epod_locs[i], epod_locs[j])
        epod_locs.pop(j)
        epod_locs[i] = (start1,end2)
        continue

    i += 1

  #print epod_locs

  # write a file containing 1s at the positions involved in epods
  epod_loc_vec = scipy.zeros(len(epod_vec))
  for i in range(len(epod_vec)):
    iloc = offsets[i]
    for (start,end) in epod_locs:
      if (start > end):
        if (end < 0):
          end += len(epod_vec)
        else:
          start -= len(epod_vec)
      if (iloc>=start and iloc<= end):
        epod_loc_vec[i] = 1
        continue

  write_grfile(offsets, epod_loc_vec, epod_outfile)


  return epod_locs

def identify_epods_v3(epod_data, percentile_data, min_epod_length, epod_outfile, lpod=False, delta=25):
  """
  Look for epods over (100 - delta-th) percentile from percentile_data in epod_data

  Each data input should be a file (NOT a vector)
  We then look for all contiguous regions of length at least min_epod_length  
  A file (epod_outfile) containing the locations of epods will also be written

  If lpod is true, then we instead look for similar regions under 25th percentile occupancy

  All units are IN BASE PAIRS

  This version of the function tries to expand each epod symmetrically starting from local maxima , which is less likely
   to yield asymmetric epods than is the original identify_epods function

  In addition, this has been modified relative to v2 to make it less tolerant of drops in the occupancy trace below 0, which
    seem like they ought to break the epod no matter what

  """

  
  import scipy.stats

  offsets, epod_vec = numpy.loadtxt(epod_data, unpack=True)
  percentile_vec = numpy.loadtxt(percentile_data, usecols=(1,), unpack=True)
  epod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, 100-delta)
  lpod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, delta)
  percentile_vec = []

  # establish a guess for how many bp between each entry
  stride = offsets[1] - offsets[0]

  def epod_cmp(value):
    if (lpod):
      return (value < lpod_cutoff)
    else:
      return (value > epod_cutoff)

  if (lpod):
    print "Searching for regions at least %i bp long with a median z score below %f" % (min_epod_length, lpod_cutoff)
  else:
    print "Searching for regions at least %i bp long with a median z score above %f" % (min_epod_length, epod_cutoff)

  epod_pot_arr = numpy.zeros_like(epod_vec) - 1000 ;# contains the values of all 1024 bp windows that are potential epods
  # we follow a two-pass approach to find epods
  # first we go through the raw data and find all 1024 bp windows where the
  #   median value is at least percentile_vec
  # Each time we find such a window, we add to epod_pot_arr the score at that location
  #  this way we know the relative heights of various windows
  #  note that we add the plain score, and not the window median, to minimize ties

  offset = int(min_epod_length / 2)
  for i in range(len(offsets)):
    start = offsets[i]-offset
    end = offsets[i]+offset
    curr_median = scipy.median(circular_range_bps(offsets,epod_vec, start, end)[1])
    if (epod_cmp(curr_median)):
      epod_pot_arr[i] = epod_vec[i]

  # now, we find the BEST locations to start potential epods, and try expanding them in either direction

  #numpy.save('test_centers.npy',epod_pot_arr)
  epod_pot_centers = scipy.signal.argrelextrema(epod_pot_arr, numpy.greater_equal,mode='wrap')[0]
  epod_abovezero = numpy.argwhere(epod_pot_arr > -100)

  epod_centers = numpy.intersect1d(epod_pot_centers, epod_abovezero)

  epod_locs = []

  #print "DEBUG: Potential epod centers found at: %s" % epod_centers


  for center_i in epod_centers:
      # we start at just the centers, which are peaks in the occupancy trace, and make sure that we can expand to
      #  be large enough for an epod without hitting any zeroes
      epod_start = offsets[center_i] - 1
      epod_end = offsets[center_i] + 1

      expand_left = True
      expand_right = True

      while (expand_left or expand_right):
          # we try expanding this epod as far as we can

          if expand_left:
              # try expanding to the left
              # we break if either the window median drops too low, or the value at the position of interest drops below 0
              trial_median = scipy.median(circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end)[1])
              new_value = circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end)[1][0]
              if new_value < 0:
                  expand_left = False
              elif trial_median > epod_cutoff:
                  epod_start -= 1
              else:
                  expand_left = False

          if expand_right:
              # try expanding to the right
              trial_median = scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end+1)[1])
              new_value = circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end)[1][-1]
              if new_value < 0:
                  expand_right = False
              elif trial_median > epod_cutoff:
                  epod_end += 1
              else:
                  expand_right = False
      
      #print "DEBUG: After expansion, ipod between %i and %i has median %f" % (epod_start, epod_end, scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end)[1]) )
      epod_locs.append( (int(epod_start), int(epod_end)) )
  
  # Take one shot at merging adjacent ipods
  i = 0
  while (i < (len(epod_locs) - 2)):
    j = i+1
    start1, end1 = epod_locs[i]
    start2, end2 = epod_locs[j]
    if (end1 > start2):
      #print "DEBUG: Trying to merge epods %s and %s" % (epod_locs[i], epod_locs[j])
      #print scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])
      if epod_cmp(scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])):
        #print "DEBUG: Merging epods %s and %s" % (epod_locs[i], epod_locs[j])
        epod_locs.pop(j)
        epod_locs[i] = (min(start1,start2),max(end1,end2))
        continue

    i += 1

  #print "DEBUG:" 
  #print epod_locs

  # write a file containing 1s at the positions involved in epods
  epod_loc_vec = scipy.zeros(len(epod_vec))
  for i in range(len(epod_vec)):
    iloc = offsets[i]
    for (start,end) in epod_locs:
      if numpy.abs(start-end) < min_epod_length:
          continue

      if (start > end):
        if (end < 0):
          end += len(epod_vec)
        else:
          start -= len(epod_vec)
      if (iloc>=start and iloc<= end):
        epod_loc_vec[i] = 1
        continue

  write_grfile(offsets, epod_loc_vec, epod_outfile)


  return epod_locs

def identify_epods_v4(epod_data, percentile_data, min_epod_length, epod_outfile, lpod=False, delta=25):
  """
  Look for epods over (100 - delta-th) percentile from percentile_data in epod_data

  Each data input should be a file (NOT a vector)
  We then look for all contiguous regions of length at least min_epod_length  
  A file (epod_outfile) containing the locations of epods will also be written

  If lpod is true, then we instead look for similar regions under 25th percentile occupancy

  All units are IN BASE PAIRS

  This version of the function tries to expand each epod symmetrically starting from local maxima , which is less likely
   to yield asymmetric epods than is the original identify_epods function

  this version is similar to v2 calling, but we trim the very ends of EPODs so that
    there are no values less than a certain threshold on the ends
    (in this version, though, they are allowed internally)

  """

  
  import scipy.stats

  offsets, epod_vec = numpy.loadtxt(epod_data, unpack=True)
  percentile_vec = numpy.loadtxt(percentile_data, usecols=(1,), unpack=True)
  epod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, 100-delta)
  epod_cutoff_lo = scipy.stats.scoreatpercentile(percentile_vec, 100-(2*delta))
  lpod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, delta)
  percentile_vec = []

  # establish a guess for how many bp between each entry
  stride = offsets[1] - offsets[0]

  def epod_cmp(value):
    if (lpod):
      return (value < lpod_cutoff)
    else:
      return (value > epod_cutoff)

  if (lpod):
    print "Searching for regions at least %i bp long with a median z score below %f" % (min_epod_length, lpod_cutoff)
  else:
    print "Searching for regions at least %i bp long with a median z score above %f" % (min_epod_length, epod_cutoff)

  epod_pot_arr = numpy.zeros_like(epod_vec) - 1000 ;# contains the values of all 1024 bp windows that are potential epods
  # we follow a two-pass approach to find epods
  # first we go through the raw data and find all 1024 bp windows where the
  #   median value is at least percentile_vec
  # Each time we find such a window, we add to epod_pot_arr the score at that location
  #  this way we know the relative heights of various windows
  #  note that we add the plain score, and not the window median, to minimize ties

  offset = int(min_epod_length / 2)
  for i in range(len(offsets)):
    start = offsets[i]-offset
    end = offsets[i]+offset
    curr_median = scipy.median(circular_range_bps(offsets,epod_vec, start, end)[1])
    if (epod_cmp(curr_median)):
      epod_pot_arr[i] = epod_vec[i]


  #numpy.save('test_centers.npy',epod_pot_arr)
  #epod_pot_centers = scipy.signal.argrelextrema(epod_pot_arr, numpy.greater_equal,mode='wrap')[0]
  epod_centers = numpy.argwhere(epod_pot_arr > -100)

  epod_locs = []

  #print "DEBUG: Potential epod centers found at: %s" % epod_centers


  for center_i in epod_centers:
      # we start at just the centers, which are peaks in the occupancy trace, and make sure that we can expand to
      #  be large enough for an epod without hitting any zeroes
      epod_start = offsets[center_i] - offset
      epod_end = offsets[center_i] + offset

      expand_left = True
      expand_right = True

      while (expand_left or expand_right):
          # we try expanding this epod as far as we can

          if expand_left:
              # try expanding to the left
              trial_median = scipy.median(circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end)[1])
              if trial_median > epod_cutoff:
                  epod_start -= 1
              else:
                  expand_left = False

          if expand_right:
              # try expanding to the right
              trial_median = scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end+1)[1])
              if trial_median > epod_cutoff:
                  epod_end += 1
              else:
                  expand_right = False

      # now, trim back th ends to remove any zero points
      this_vals = circular_range_bps(offsets,epod_vec, epod_start, epod_end)[1]
      start_i = 0
      end_i = -1
      while(this_vals[start_i] < epod_cutoff_lo):
          start_i += 1

      while(this_vals[end_i] < epod_cutoff_lo):
          end_i -= 1

      epod_start +=  start_i
      epod_end += end_i


      
      print "DEBUG: After expansion, ipod between %i and %i has median %f" % (epod_start, epod_end, scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end)[1]) )
      epod_locs.append( (int(epod_start), int(epod_end)) )
  
  # Take one shot at merging adjacent ipods
  i = 0
  while (i < (len(epod_locs) - 2)):
    j = i+1
    start1, end1 = epod_locs[i]
    start2, end2 = epod_locs[j]
    if (end1 > start2):
      print "DEBUG: Trying to merge epods %s and %s" % (epod_locs[i], epod_locs[j])
      print scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])
      if epod_cmp(scipy.median(circular_range_bps(offsets,epod_vec, start1, end2)[1])):
        print "DEBUG: Merging epods %s and %s" % (epod_locs[i], epod_locs[j])
        epod_locs.pop(j)
        epod_locs[i] = (min(start1,start2),max(end1,end2))
        continue

    i += 1

  print "DEBUG:" 
  print epod_locs

  # write a file containing 1s at the positions involved in epods
  epod_loc_vec = scipy.zeros(len(epod_vec))
  for i in range(len(epod_vec)):
    iloc = offsets[i]
    for (start,end) in epod_locs:
      if numpy.abs(start-end) < min_epod_length:
          continue

      if (start > end):
        if (end < 0):
          end += len(epod_vec)
        else:
          start -= len(epod_vec)
      if (iloc>=start and iloc<= end):
        epod_loc_vec[i] = 1
        continue

  write_grfile(offsets, epod_loc_vec, epod_outfile)


  return epod_locs

def identify_epods_v4a(epod_data, percentile_data, min_epod_length, epod_outfile, lpod=False, delta=25,genomelength=4641652):
  """
  Look for epods over (100 - delta-th) percentile from percentile_data in epod_data

  Each data input should be a file (NOT a vector)
  We then look for all contiguous regions of length at least min_epod_length  
  A file (epod_outfile) containing the locations of epods will also be written

  If lpod is true, then we instead look for similar regions under 25th percentile occupancy

  All units are IN BASE PAIRS

  This version of the function tries to expand each epod symmetrically starting from local maxima , which is less likely
   to yield asymmetric epods than is the original identify_epods function

  this version is similar to v2 calling, but we trim the very ends of EPODs so that
    there are no values less than a certain threshold on the ends
    (in this version, though, they are allowed internally)

  """

  
  import scipy.stats

  offsets, epod_vec = numpy.loadtxt(epod_data, unpack=True)
  percentile_vec = numpy.loadtxt(percentile_data, usecols=(1,), unpack=True)
  epod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, 100-delta)
  epod_cutoff_lo = scipy.stats.scoreatpercentile(percentile_vec, 100-(2*delta))
  lpod_cutoff = scipy.stats.scoreatpercentile(percentile_vec, delta)
  percentile_vec = []

  # establish a guess for how many bp between each entry
  stride = offsets[1] - offsets[0]

  def epod_cmp(value):
    if (lpod):
      return (value < lpod_cutoff)
    else:
      return (value > epod_cutoff)

  if (lpod):
    print "Searching for regions at least %i bp long with a median z score below %f" % (min_epod_length, lpod_cutoff)
  else:
    print "Searching for regions at least %i bp long with a median z score above %f" % (min_epod_length, epod_cutoff)

  epod_pot_arr = numpy.zeros_like(epod_vec) - 1000 ;# contains the values of all 1024 bp windows that are potential epods
  # we follow a two-pass approach to find epods
  # first we go through the raw data and find all 1024 bp windows where the
  #   median value is at least percentile_vec
  # Each time we find such a window, we add to epod_pot_arr the score at that location
  #  this way we know the relative heights of various windows
  #  note that we add the plain score, and not the window median, to minimize ties

  offset = int(min_epod_length / 2)
  for i in range(len(offsets)):
    start = offsets[i]-offset
    end = offsets[i]+offset
    curr_median = scipy.median(circular_range_bps(offsets,epod_vec, start, end,genomelength=genomelength)[1])
    if (epod_cmp(curr_median)):
      epod_pot_arr[i] = epod_vec[i]


  #numpy.save('test_centers.npy',epod_pot_arr)
  #epod_pot_centers = scipy.signal.argrelextrema(epod_pot_arr, numpy.greater_equal,mode='wrap')[0]
  epod_centers = numpy.argwhere(epod_pot_arr > -100)

  epod_locs = []

  #print "DEBUG: Potential epod centers found at: %s" % epod_centers


  for center_i in epod_centers:
      # we start at just the centers, which are peaks in the occupancy trace, and make sure that we can expand to
      #  be large enough for an epod without hitting any zeroes

      # skip if the center of this epod is already contained in the previous one
      if (len(epod_locs) > 1) and (offsets[center_i] < epod_locs[-1][1]):
          continue

      epod_start = offsets[center_i] - offset
      epod_end = offsets[center_i] + offset

      expand_left = True
      expand_right = True

      start_index_l = 0
      end_index_l = len(offsets)
      start_index_r = 0
      end_index_r = len(offsets)
      print "starting to do expansion"
      while (expand_left or expand_right):
          # we try expanding this epod as far as we can

          if expand_left:
              # try expanding to the left
              slice_offsets,trial_vals,start_index_l, end_index_l = circular_range_bps(offsets,epod_vec, epod_start - 1, epod_end,genomelength=genomelength,return_inds=True, guess_start = start_index_l, guess_end=end_index_l)
              trial_median = scipy.median(trial_vals)
              if trial_median > epod_cutoff:
                  epod_start -= 1
              else:
                  expand_left = False

          if expand_right:
              # try expanding to the right
              slice_offsets,trial_vals,start_index_r, end_index_r = circular_range_bps(offsets,epod_vec, epod_start, epod_end+1,genomelength=genomelength,return_inds=True, guess_start = start_index_r, guess_end=end_index_r)
              trial_median = scipy.median(trial_vals)
              if trial_median > epod_cutoff:
                  epod_end += 1
              else:
                  expand_right = False

      print "starting to do trimming"
      # now, trim back th ends to remove any zero points
      this_vals = circular_range_bps(offsets,epod_vec, epod_start, epod_end,genomelength=genomelength)[1]
      start_i = 0
      end_i = -1
      while(this_vals[start_i] < epod_cutoff_lo):
          start_i += 1

      while(this_vals[end_i] < epod_cutoff_lo):
          end_i -= 1

      epod_start +=  start_i
      epod_end += end_i


      
      print "DEBUG: After expansion, ipod between %i and %i has median %f" % (epod_start, epod_end, scipy.median(circular_range_bps(offsets,epod_vec, epod_start, epod_end,genomelength=genomelength)[1]) )
      epod_locs.append( (int(epod_start), int(epod_end)) )
  
  # Take one shot at merging adjacent ipods
  i = 0
  while (i < (len(epod_locs) - 2)):
    j = i+1
    start1, end1 = epod_locs[i]
    start2, end2 = epod_locs[j]
    if (end1 > start2):
      print "DEBUG: Trying to merge epods %s and %s" % (epod_locs[i], epod_locs[j])
      print scipy.median(circular_range_bps(offsets,epod_vec, start1, end2,genomelength=genomelength)[1])
      if epod_cmp(scipy.median(circular_range_bps(offsets,epod_vec, start1, end2,genomelength=genomelength)[1])):
        print "DEBUG: Merging epods %s and %s" % (epod_locs[i], epod_locs[j])
        epod_locs.pop(j)
        epod_locs[i] = (min(start1,start2),max(end1,end2))
        continue

    i += 1

  print "DEBUG:" 
  print epod_locs

  # write a file containing 1s at the positions involved in epods
  epod_loc_vec = scipy.zeros(len(epod_vec))
  for i in range(len(epod_vec)):
    iloc = offsets[i]
    for (start,end) in epod_locs:
      if numpy.abs(start-end) < min_epod_length:
          continue

      if (start > end):
        if (end < 0):
          end += len(epod_vec)
        else:
          start -= len(epod_vec)
      if (iloc>=start and iloc<= end):
        epod_loc_vec[i] = 1
        continue

  write_grfile(offsets, epod_loc_vec, epod_outfile)


  return epod_locs


def analyze_tf_occupancy(annfile, ipoddat, tfstr, outfile=None):
  """
  Print or write to a file some notes on ipod occupancy at TF binding sites

  The data that will be written are the name, start, end, average score, and whether the site is occupied by anything else (Y/N)
  """

  from sys import stdout
  import re

  if (outfile):
    outstr = open(outfile, "w")
  else:
    outstr = stdout

  tfregexp = re.compile(tfstr)

  # First get a list of ranges containing binding sites for the specified tf
  my_bindsites, my_anns = get_range_matching_annotations(annfile, tfregexp, ipoddat, True)

  all_locs = numpy.loadtxt(annfile, unpack=True, usecols=(3,4), delimiter="\t")

  # load all of the data
  locs, data = numpy.loadtxt(ipoddat, unpack=True)

  # Now run through the annotations and write data on each of them
  for bs,ann in zip(my_bindsites, my_anns):
    #print bs
    site_start, site_end = bs
    site_scores = data[site_start:(site_end+1)]
    site_meanscore = scipy.mean(site_scores)
    #print "Looking for overlaps"
    overlaps = find_overlaps_bps((locs[bs[0]], locs[bs[1]]), all_locs)
    #print overlaps
    overlapping = "N"
    if (len(overlaps) > 0):
      overlapping = "Y"

    outstr.write("%s\t%i\t%i\t%f\t%s\n" % (ann.replace("\n", ""), locs[site_start], locs[site_end], site_meanscore, overlapping))

  if (outfile):
    outstr.close()


def check_occupancy_at_annotations(annfile, ipoddat, outfile):
  """
  For each annotation in annfile, write the average z score to outfile

  While this will generally be used with z-scores, it can be used with
  any .gr-formatted file
  """

  annstr = open(annfile, 'r')
  locs, data = numpy.loadtxt(ipoddat, unpack=True)
  outstr = open(outfile, 'w')

  print "Done loading"

  for line in annstr:
    linearr = line.split()
    mystart = int(linearr[3])
    myend = int(linearr[4])
    mystrands = linearr[6]
    myname = linearr[8]

    if (mystart > myend):
      tmp=mystart
      mystart=myend
      myend = tmp

    startind = -1
    endind = -1
    for i in range(len(locs)):
      if (locs[i] < mystart):
        startind = i+1
      if (locs[i] < myend):
        endind = i+1

    print "start/end: %i/%i" % (mystart, myend)
    print "Range for %s is %i - %i" % (myname, startind, endind)

    myavg = scipy.mean(data[startind:endind])
    outstr.write("%s (%i - %i): %f\n" % (myname, mystart, myend, myavg))

  outstr.close()
  annstr.close()

def get_range_matching_annotations(annotationfile, annotationregexp, seqfile, return_annotations = False):
  """
  Return a list of ranges from a .gr file matching an annotation regexp

  Returned values are in units of bp

  If return_annotations is true, the annotation at each position is also returned
  """

  import re

  instr = open(annotationfile, 'r')
  offsets = numpy.loadtxt(seqfile, unpack=True, usecols=(0,))
  num_offsets = len(offsets)

  ranges = []
  annotations = []
  numcan = 0

  for line in instr:
    linearr = line.split("\t")
    annotation = " ".join(linearr[8:])
    if (annotationregexp.search(annotation)):
      numcan += 1
      #print "Looking for site (%i candidate)" % numcan
      start = int(linearr[3])
      end = int(linearr[4])
      if (start > end):
        tmp=start
        start=end
        end=tmp
      start_gr = 0
      end_gr = 0

      jumpsize = num_offsets / 2
      curr_ind = jumpsize
      jumpsize /= 2
      
      while (jumpsize > 4):
        #print "Looking for start %s at currval %s with jumpsize %s" % (start, offsets[curr_ind], jumpsize)
        if (start < offsets[curr_ind]):
          curr_ind -= jumpsize
        else:
          curr_ind += jumpsize

        jumpsize /= 2

      while (offsets[start_gr] > start):
        start_gr -= jumpsize

      while (offsets[start_gr] < start):
        start_gr += 1

      #print "Found value %s at %s" % (offsets[start_gr], start_gr)

      if (offsets[start_gr] > end):
        continue

      end_gr = start_gr
      while (offsets[end_gr] <= end and end_gr < len(offsets)):
        #print "Looking for end %s at currval %s" % (end, offsets[end_gr])
        end_gr += 1

      if (offsets[end_gr] > end):
        end_gr -= 1

      if (offsets[end_gr] - offsets[start_gr] > 1000):
        continue

      #print "Using values %s/%s for %s/%s" % (offsets[start_gr], offsets[end_gr], start, end)

      ranges.append( (start_gr, end_gr) )
      annotations.append(annotation)

      #print len(annotations)


  if (return_annotations):
    return (ranges, annotations)
  else:
    return(ranges)

def find_annotations_at_positions(genomefile, rangelist, seqconvfile, genomelength=4641652, enclosed=True, offsets_bp = False):
  """
  Return full information on annotations in each range in rangelist

  Primarily intended as a helper function for write_annotations_at_positions

  Arguments:
    genomefile -- a .gff file containing the annotations of interest
    rangelist -- a list of (start,end) tuples corresponding to the
      index ranges in seqconvfile which should be considered
    seqconvfile -- a .gr file where the first column contains bp positions
    genomelengh -- the length of the genome (should be very large if not circular)
    enclosed -- a boolean flag. If true, only annotations fully enclosed
      by the range will be returned. If false, any overlap will 
      cause a peak to be written

    offsets_bp: true if the input pairs are in bp instead of indices

  Return values:
    A sorted list containing rangestart,rangeend,annotation_start,annotation_end,
      annotation_strand,annotation_name tuples
    rangestart and rangeend are in bp, not the original index units
    The list is sorted by rangestart
  """
  genomefilestr = open(genomefile, 'r')
  seqnums = numpy.loadtxt(seqconvfile, unpack=True, usecols=(0,))

  annotated_epods = []

  for line in genomefilestr:
    linearr = line.split("\t")
    mystart = int(linearr[3])
    myend = int(linearr[4])
    mystrand = linearr[6]
    myname = linearr[8].rstrip()


    if (mystart > myend):
      tmp=mystart
      mystart=myend
      myend=tmp

    if (mystart < 3000):
      print "Working on %s" % line

    for (start_nonconv,end_nonconv) in rangelist:
      if (not offsets_bp):
        start = seqnums[start_nonconv]
        end = seqnums[end_nonconv]
      else:
        start = start_nonconv
        end = end_nonconv

      if (enclosed):
        if (mystart >= start and myend <= end):
          annotated_epods.append((start, end, mystart, myend, mystrand, myname))
          continue
      else:
        #print "Checking against %i/%i" % (start, end)
        if ( (mystart >= start and mystart <= end) or (myend >= start and myend <= end) or (mystart <= start and myend >= end)):
          annotated_epods.append((start, end, mystart, myend, mystrand, myname))
          continue


  genomefilestr.close()

  annotated_epods.sort(lambda x, y : cmp(x[0], y[0]))

  return annotated_epods

def write_annotations_at_positions(genomefile, rangelist, seqconvfile, outfile, genomelength=4641652, enclosed=True, offsets_bp = False):
  """
  For each range of start,end pairs, write the annotations enclosed by the given ranges

  rangelist should be a list of start,end pairs
  genomefile is the file containing the annotations of interest
  genomelength is the total length of the genome in bp
  seqconvfile is a file containing the offsets (in bp) corresponding to
    pairs in rangelist
  enclosed is the enclosed flag that will be passed to 
    find_annotations_at_positions
    offsets_bp: true if the input pairs are in bp instead of indices
  """

  import re

  genere = re.compile('Gene (\w+);')

  genomefilestr = open(genomefile, 'r')
  seqnums = numpy.loadtxt(seqconvfile, unpack=True, usecols=(0,))

  annotated_epods = find_annotations_at_positions(genomefile, rangelist, seqconvfile, genomelength, enclosed=enclosed, offsets_bp = offsets_bp)

  #  for line in genomefilestr:
  #    linearr = line.split()
  #    mystart = int(linearr[3])
  #    myend = int(linearr[4])
  #    mystrand = linearr[6]
  #    myname = linearr[8]
  #
  #    if (mystart > myend):
  #      tmp=mystart
  #      mystart=myend
  #      myend=tmp
  #
  #    for (start_nonconv,end_nonconv) in rangelist:
  #      start = seqnums[start_nonconv]
  #      end = seqnums[end_nonconv]
  #      if (mystart >= start and myend <= end):
  #        annotated_epods.append(start, end, mystart, myend, mystrand, myname)
  #        continue
  #
  #  genomefilestr.close()
  #
  #  annotated_epods.sort(lambda x, y : cmp(x[0], y[0]))

  currstart = -1
  currend = -1

  epodtab = table.Table()
  epodtab.init_keys(["Region start", "Region end", "Annotations"])

  thisrow = None

  for epod in annotated_epods:
    start, end, mystart, myend, mystrand, myname = epod
    if (genere.search(myname)):
      myname = genere.search(myname).groups()[0]

    if ((start != currstart) or (end != currend)):

      if (thisrow is not None):
        thisrow["Annotations"] = ",".join(thisrow["Annotations"])
        epodtab.add_row(thisrow)

      thisrow = {"Region start" : start, "Region end" : end, "Annotations" : [myname]}
      currstart = start
      currend = end

    else:
      thisrow["Annotations"] += [myname]

  table.write_table(outfile, epodtab, sep="\t")


def identify_peaks(infile, medfile, minlength, maxlength, percentile, outfile, troughs = False, genlength =4641652 , outfile_startend = None):
  """
  Identify regions in the specified length range with high occupancy

  The peaks will be taken from the data in infile
  The data will be searched for regions between $minlength and $maxlength bp
    long that have median occupancy at least that of $percentile from $medfile

  The identified peaks will both be returned and written to a file

  If troughs is True, we instead look for regions with *low* occupancy

  All output is in index units, NOT bp units
  """

  import scipy.stats
  import bisect
  import copy

  offsets, datvec = numpy.loadtxt(infile, unpack=True)
  percdat = numpy.loadtxt(medfile, usecols=(1,), unpack=True)
  perc_cutoff = scipy.stats.scoreatpercentile(percdat, percentile)
  percdat = []


  #print "Looking for peaks of minimum length %i and maximum length %i with minimum value %f" % (minlength, maxlength, perc_cutoff)



  def peak_cmp(value):
    if (troughs):
      return (value < perc_cutoff)
    else:
      return (value > perc_cutoff)


  # Now go through the data. At each point, we look at a window 20 bp long
  #  starting at that position, and if it has the required median value,
  #  then we extend probe by probe, until the median drops below the required value
  startind = 0
  peaklocs = []
  while (startind < len(datvec)):
    #print "Looking at startind %i" % startind
    #if (startind < 1000):
    #  print "on startind %i" % startind

    offset = 0
    while offsets[startind + offset] - offsets[startind] < minlength:
      offset += 1
      if (startind + offset) == len(offsets):
        while(genlength + offsets[startind + offset - len(offsets)] - offsets[startind] < minlength):
          offset += 1
        break
    offset -= 1

    firstind_peak = startind
    lastind_peak = startind + offset

    curr_offsets = circular_range(offsets, startind, startind+offset)
    curr_data = circular_range(datvec, startind, startind+offset)
    curr_data_sorted = list(scipy.sort(curr_data))
    curr_median = sorted_list_median(curr_data_sorted)
    
    #print "Considering range %i to %i" % (firstind_peak, lastind_peak)
    #print peaklocs


    if peak_cmp(curr_median):
      # find the full extent of this peak
      # First expand as far as possible upstream

      while (1):
        if (firstind_peak >= len(offsets)):
          break

        trial_ind = firstind_peak - 1
        if trial_ind < 0:
          trial_ind = len(datvec) - 1

        if (lastind_peak >= len(offsets)):
          lastind_peak -= len(offsets)

        if (lastind_peak < firstind_peak):
          corr = genlength
        else:
          corr = 0

        if (offsets[lastind_peak] + corr - offsets[firstind_peak - 1] > maxlength): 
          break

        trial_value = datvec[trial_ind]
        #print curr_data_sorted
        #print "New value is %s" % trial_value
        #trial_data_sorted = copy.deepcopy(curr_data_sorted)
        #bisect.insort_right(trial_data_sorted, trial_value)
        newloc = bisect.bisect_right(curr_data_sorted, trial_value)
        curr_data_sorted.insert(newloc,trial_value)
        new_median = sorted_list_median(curr_data_sorted)
        #print curr_data_sorted


        if peak_cmp(new_median):
          #print "New firstind %i is ok (median %f)" % (trial_ind, new_median)
          firstind_peak = trial_ind
          curr_offsets = scipy.insert(curr_offsets, 0, offsets[firstind_peak])
          curr_data = scipy.insert(curr_data, 0, datvec[firstind_peak])
          #curr_data_sorted = trial_data_sorted
          #print "Compare %f with %f" % (new_median, scipy.median(curr_data))
          #print curr_data
          #print curr_data_sorted
        else:
          curr_data_sorted.pop(newloc)
          break

      # now try to move upwards
      while (1):
        #print "Considering range %i to %i" % (firstind_peak, lastind_peak)
        trial_ind = lastind_peak + 1
        if trial_ind == (len(datvec)-1):
          trial_ind = 0

        if (lastind_peak + 1 >= len(offsets)):
          lastind_peak -= len(offsets)


        if (lastind_peak < firstind_peak):
          corr = genlength
        else:
          corr = 0

        if (offsets[lastind_peak+1]+corr - offsets[firstind_peak] > maxlength): 
          break

        trial_value = datvec[trial_ind]
        #trial_data_sorted = curr_data_sorted
        #bisect.insort_right(trial_data_sorted, trial_value)
        newloc = bisect.bisect_right(curr_data_sorted, trial_value)
        curr_data_sorted.insert(newloc,trial_value)
        new_median = sorted_list_median(curr_data_sorted)

        if peak_cmp(new_median):
          lastind_peak = trial_ind
          curr_offsets = scipy.append(curr_offsets, offsets[lastind_peak])
          curr_data = scipy.append(curr_data, datvec[lastind_peak])
          #curr_data_sorted = trial_data_sorted
          #print "Compare %f with %f" % (new_median, scipy.median(curr_data))
          #print curr_data
          #print curr_data_sorted
        else:
          curr_data_sorted.pop(newloc)
          break


      #print "Appending peak %i %i" % (firstind_peak, lastind_peak)
      peaklocs.append( (firstind_peak, lastind_peak) )
      if (lastind_peak < startind):
        #print "Ending because %i < %i" % (lastind_peak, startind)
        break
      else:
        startind = lastind_peak + 1

    else:
      startind += 1

    #print "New startind will be %i" % startind


  # take one shot at merging adjacent peaks
  i = 0
  while (i < (len(peaklocs) - 2)):
    j = i+1
    start1, end1 = peaklocs[i]
    start2, end2 = peaklocs[j]
    if (end1 >= start2 and ((end2 - start1) > maxlength)):
      #print "Trying to merge peaks %s and %s" % (peaklocs[i], peaklocs[j])
      #print scipy.median(circular_range(datvec, start1, end2))
      if peak_cmp(scipy.median(circular_range(datvec, start1, end2))):
        #print "Merging peaks %s and %s" % (peaklocs[i], peaklocs[j])
        peaklocs.pop(j)
        peaklocs[i] = (start1,end2)
        continue
  
    i += 1


  # write an array with the peaks flagged
  numpos = len(datvec)
  peak_loc_vec = scipy.zeros(numpos)

  for (start,end) in peaklocs:
    if (end > numpos):
      peak_loc_vec[start:] = 1
      peak_loc_vec[:(end - numpos)] = 1
    else:
      peak_loc_vec[start:end] = 1

  write_grfile(offsets, peak_loc_vec, outfile)

  if (outfile_startend):
    outstr = open(outfile_startend, 'w')
    for loc in peaklocs:
      outstr.write("%i %i\n" % loc)

    outstr.close()

  return peaklocs

  
       
def calc_epod_rna_stats(epodlocs, rnafile, genomelength, outfile):
  """
  Write a table containing statistics on epods
  
  For each identified epod, we will write the start and end positions (in bp),
    mean, median, and log2 median RNA occupancy (read from the data in rnagrfile),
    and the 2-means cluster id for this value (0 for tsepod, 1 for hepod)
  
  Arguments:
    epodlocs -- list of start,end pairs containing the epods of interest
      should be in units corresponding to the offsets in rnagrfile 
    rnagrfile -- a .gr file containing the bp offsets and rna scores
    genomelength -- the length of the genome (if circular)
    outfile -- name of the file to which output should be written

  Return values:
    None. A table is written to outfile.
  """

  from scipy.cluster.vq import *

  rnaoffsets, rnavals = read_grfile(rnafile)

  epod_starts = []
  epod_ends = []
  epod_mean_rnas = []
  epod_median_rnas = []
  epod_log_median_rnas = []
  epod_annotations = []


  for start, end in epodlocs:
    #print "considering peak %i %i" % (start,end)
    epod_starts.append(rnaoffsets[start])
    epod_ends.append(rnaoffsets[end])
    epod_rna_scores = rnavals[start:end+1]
    epod_mean_rnas.append(numpy.mean(epod_rna_scores))
    epod_median_rnas.append(numpy.median(epod_rna_scores))
    epod_log_median_rnas.append(numpy.median(scipy.log2(epod_rna_scores)))


  # Do k-means clustering to find hepods vs tsepods based on median occupancy
  means, distortion = kmeans(epod_median_rnas, 2)
  means = scipy.sort(means)
  epod_expression_flag = 1 * (scipy.greater( scipy.absolute(allvars - means[0]), scipy.absolute(allvars-means[1])) )

  outtab = table.Table()
  #outtab.initialize_columns(["Start", "End","Mean", "Median", "Log2 Median","ORFs"])
  outtab.add_column(epod_starts, "Start (bp)")
  outtab.add_column(epod_ends, "End (bp)")
  outtab.add_column(epod_mean_rnas, "Mean")
  outtab.add_column(epod_median_rnas, "Median")
  outtab.add_column(epod_log_median_rnas, "Log2 Median")
  outtab.add_column(epod_expression_flag, "Expression flag")
  #outtab.add_column(epod_annotations, "ORFs")

  table.write_table(outfile, outtab, "\t")



def mask_array(data, masklocs):
  """
  Return an array with the specified regions removed
  """

  maskarr = scipy.ones(len(data))
  for begin,end in masklocs:
    maskarr[begin:end] = 0

  return data[scipy.greater(maskarr, 0)]


def calculate_corrcoefs_all(filelist, namelist):
  """
  Calculate correlation coefficient between all pairs of listed files
  """

  for i in range(len(filelist)):
    for j in range(i+1,len(filelist)):
      print "DEBUG: %s %s" % (filelist[i], filelist[j])
      calculate_corrcoeffs(filelist[i], filelist[j], namelist[i], namelist[j])

def calculate_corrcoeffs(file1, file2, name1, name2):
  """
  Calculate the correlation coefficient for a pair of files and print it
  """

  file1_filtered2 = mktemp()
  file2_filtered1 = mktemp()

  filter_grfile_by_bpranges(file2, file1, file1_filtered2)
  filter_grfile_by_bpranges(file1_filtered2, file2, file2_filtered1)

  mat1 = numpy.loadtxt(file1_filtered2, unpack=True, usecols=(1,))
  mat2 = numpy.loadtxt(file2_filtered1, unpack=True, usecols=(1,))
  os.remove(file1_filtered2)
  os.remove(file2_filtered1)
#    mat2[scipy.less(mat2, minval)] = minval
#    mat2[scipy.greater(mat2, maxval)] = maxval
#    mat2 = ipod_utils.mask_array(mat2, mask_ranges)

  print "%s - %s: %f" % (name1, name2, scipy.corrcoef(mat1, mat2)[1,0])

def hierclust_arrays(allfiles, allnames, outfile):
  """
  Do hierarchical clustering on a set of .gr files
  """

  allmats = []
  for file in allfiles:
    allmats.append(numpy.loadtxt(file, unpack=True, usecols=(1,)))

  allmats = scipy.array(allmats)
  Y = hcluster.pdist(allmats, 'correlation')
  Z=hcluster.linkage(Y,'single')
  hcluster.dendrogram(Z, labels=allnames, leaf_rotation=90, leaf_font_size=8)
  pylab.savefig(outfile)

def normalize_mean_value(infile, outfile, targetval=2500.0):
  """
  Write a new .gr file in which the old file is scaled to a mean value of targetval
  """

  offsets, data = read_grfile(infile)

  curr_mean = scipy.mean(data)
  data *= (targetval / curr_mean)

  write_grfile(offsets, data, outfile)

def normalize_mean_value_withoffset(infile, outfile, targetval=2500.0,offset=0.25):
  """
  Write a new .gr file in which the old file is scaled to a mean value of targetval

  We apply  a small offset to avoid zero entries -- note that we calculate the weighting before the
    offset is applied, so that we're normalizing only the fraction of the values that come from the reads themselves
  """

  offsets, data = read_grfile(infile)

  curr_mean = scipy.mean(data)
  data *= ( (targetval-offset) / curr_mean)
  data += offset

  write_grfile(offsets, data, outfile)

def normalize_by_countfile(infile, countfile, outfile, countscale=10000, offset=0.25):
  """
  Write a new .gr file in which the old file is scaled by a total count
  The count is half of the number of read ends counted in countfile
  We also apply a small offset to each data point (this is useful if one is going to take ratios,
   to prevent zeroes)
  We also multiply by countscale at the end to avoid getting very small/large numbers
  """

  locs, data = read_grfile(infile)
  count_locs, count_data = read_grfile(countfile)
  fullcount = float(numpy.sum(count_data))

  newdat = (data+offset)/fullcount

  write_grfile(locs, newdat*countscale, outfile)

def offset_mean_value(infile, outfile, targetval=2500.0):
  """
  Write a new .gr file in which the old file is offset to a mean value of targetval
  """

  offsets, data = read_grfile(infile)

  curr_mean = scipy.mean(data)
  data = data + targetval - curr_mean

  write_grfile(offsets, data, outfile)

def normalize_mean_value_withgenome(infile, outfile, targetval=2500.0):
  """
  Write a new file in which the old file is scaled to a mean value of targetval

  Here we assume that the the file has genome/offset/value lines
  """

  genomes, offsets, data = numpy.loadtxt(infile, unpack=True, dtype="|S8")

  data = data.astype(float)

  curr_mean = scipy.mean(data)

  data *= (targetval / curr_mean)

  ostr = open(outfile, "w")

  for i in range(len(genomes)):
    ostr.write("%s %s %f\n" % (genomes[i], offsets[i], data[i]))

  ostr.close()


def normalize_median_value(infile, outfile, targetval=2500.0):
  """
  Write a new .gr file in which the old file is scaled to a median value of targetval
  """

  offsets, data = read_grfile(infile)

  curr_median = scipy.median(data)
  data *= (targetval / curr_median)

  write_grfile(offsets, data, outfile)

def offset_median_value(infile, outfile, targetval=2500.0):
  """
  Write a new .gr file in which the old file is offset to a median
    value of targetval
  """

  offsets, data = read_grfile(infile)

  curr_median = scipy.median(data)
  data += (targetval-curr_median)

  write_grfile(offsets, data, outfile)

def normalize_percentile_value(infile, outfile, targetval=2500.0,percentile=50.0,offset=0.0):
  """
  Write a new .gr file in which the old file is scaled so that a given percentile in the
  distribution is equal to targetval

  the value of offset is added to each value, WHILE preserving the target median
  """

  locs, data = read_grfile(infile)

  curr_quantile = scipy.stats.scoreatpercentile(data,percentile)
  newdata = data * ((targetval-offset) / curr_quantile) + offset

  write_grfile(locs, newdata, outfile)

def find_overlaps(targetsite, targetname, allsites, allnames):
  """
  Return a list of annotations overlapping with a given target

  The returned values are a list of (start, end, name) tuples, excluding
    any identical to the target

  We assume that the entries are sorted in ascending order along the chromosome
  """

  site_start, site_end = targetsite

  overlaps = []

  # Find a chunk of the genome which *could* have overlaps
  ind_in_other = 0

  i_under = ind_in_other

  while(allsites[i_under][1] < targetsite[0] - 50):
    i_under += 1

  i_over = i_under

  while (allsites[i_over][0] <= targetsite[1]):
    i_over += 1

  for i in range(i_under, i_over):
    other_start, other_end = allsites[i]

    if (other_start == targetsite[0] and other_end == targetsite[1] and allnames[i] == targetname):
      continue

    if (other_start < targetsite[1] and other_end > targetsite[0]):
      overlaps.append((other_start, other_end, allnames[i]))

  return overlaps


def find_overlaps_bps(targetsite, allsites, guess_start = None):
  """
  Return a list of annotations overlapping with a given target

  Unlike find_overlaps, find_overlaps_bps works with input in units of bps
  targetsite is a 2-element tuple containing the start and end (in bp)
  allsites is an nx2 scipy array containing starts and ends of each sites

  guess_start is a guess at the offset of the site to search for, in index units (not bp!!!)
  """

  site_start, site_end = targetsite

  overlaps = []

  # Find a chunk of the genome which *could* have overlaps
  if (guess_start):
    i = guess_start
  else:
    i = 0

  #print allsites
  #print scipy.shape(allsites)
  while (i < len(allsites[0]) and allsites[0,i] < site_end):
    #print i
    #print allsites[:,i]
    other_start, other_end = allsites[:,i]
    if (other_start < site_end and other_end > site_start):
      overlaps.append((other_start, other_end))

    i += 1

  return overlaps

def filter_grfile_by_bpranges(refgrfile, ingrfile, outgrfile):
  """
  Write to outgrfile only those entries from ingrfile corresponding to
    entries from refgrfile
  """

  inrefstr = open(refgrfile)
  indatstr = open(ingrfile)
  ostr = open(outgrfile, "w")

  datoffset = -1
  datdata = 0
  datline = "A"

  for refline in inrefstr:
    refoffset = int(refline.split()[0])
    while (datoffset < refoffset):
      datline = indatstr.readline()
      if (datline == ""):
        break

      datoffset = int(datline.split()[0])

    if (datline == ""):
      break

    if (datoffset == refoffset):
      ostr.write(datline)

  inrefstr.close()
  indatstr.close()
  ostr.close()




def filter_annotations_by_bpranges(grfile, annfile, newannfile, gaplength=500):
  """
  Write to newannfile only those annotations which are not in a gap in the grfile

  A gap is defined as a break in the numbering of grfile that is at least
    gaplength bp long
  """

  from sys import maxint

  outstr = open(newannfile, "w")
  offsets, data = read_grfile(grfile)
  instr = open(annfile)

  gaps = []

  num_probes = len(offsets)

  # First check for a starting gap
  if (offsets[0] > gaplength):
    gaps.append((0,offsets[0]-gaplength / 2))

  # Now go through the offsets and identify any gaps
  for i in range(1,num_probes):
    if ( (offsets[i] - offsets[i-1]) > gaplength ):
      gaps.append( (offsets[i-1] + gaplength / 2, offsets[i+1] + gaplength / 2) )

  # add in a gap at the end
  gaps.append( (offsets[i-1] + gaplength/2, maxint) )

  # Now that we have the gaps, look at each tf binding site
  #   We exclude all TF binding sites with **any** overlap with a gap

  # we take advantage of the fact that both the input file and the tf file are already ordered by starting index
  gapind = 0

  for line in instr:
    linearr = line.split()
    startind = int(linearr[3])
    endind = int(linearr[4])
    
    is_good = True


    while (int((gaps[gapind])[1]) < startind and gapind < len(gaps)):
      gapind += 1

    my_gapind = gapind
    while (my_gapind < len(gaps) and gaps[my_gapind][0] < endind):
      gapstart,gapend = gaps[my_gapind]
      if ( (startind > gapstart and startind < gapend) or (endind > gapstart and endind < gapend) ):
        is_good = False
        break
      my_gapind += 1

    #for (gapstart, gapend) in gaps:
    #  if ( (startind > gapstart and startind < gapend) or (endind > gapstart and endind < gapend) ):
    #    is_good = False
    #    break

    if (is_good):
      outstr.write(line)

  outstr.close()
    
def normalize_magnitude_to_unity(ingrfile, outgrfile):
  """
  Normalize a file so that < v^2 > is 1
  """

  offsets, data = read_grfile(ingrfile)

  data /= scipy.sqrt(scipy.mean(data*data))

  write_grfile(offsets, data, outgrfile)

def diff_gr_files(grfile1, grfile2, outgrfile):
  """
  Write the difference between two gr files to a third
  """

  offsets, data1 = read_grfile(grfile1)
  offsets, data2 = read_grfile(grfile2)

  data = data1 - data2

  write_grfile(offsets, data, outgrfile)

def divide_gr_files(grfile1, grfile2, outgrfile):
  """
  Write the difference between two gr files to a third
  """

  offsets, data1 = read_grfile(grfile1)
  offsets, data2 = read_grfile(grfile2)

  data = data1 / data2

  write_grfile(offsets, data, outgrfile)

def gen_flagfile_from_grfile(grfile, outfile, threshold=0):
  """
  Write a flag file (gr formatted) containing 1 or 0 depending on whether each value is greater or less than threshold
  """

  offsets, data = read_grfile(grfile)
  newdata = scipy.zeros(scipy.shape(data)) + 1 * scipy.greater(data, threshold)
  write_grfile(offsets, newdata, outfile)
  
def gen_flagfile_from_gff(gff, grfile, outfile):
  """
  Write a gr file containing 1 at positions with annotations in the gff file and 0 elsewhere

  The set of offsets in grfile will be considered

  This will go fastest of the contents of the gff file are sorted, and the gr file *must* be
  """

  offsets, grdata  = read_grfile(grfile)

  all_locs = numpy.loadtxt(gff, unpack=True, usecols=(3,4),dtype="int", delimiter="\t")
  all_names = numpy.loadtxt(gff, unpack=True, usecols=(0,), dtype="|S8", delimiter="\t")

  flagdat = scipy.zeros(scipy.shape(offsets))

  curr_offset = 0
  num_offsets = len(offsets)

  dummy, numentries = scipy.shape(all_locs)

  for i in range(numentries):
    start = all_locs[0,i]
    end = all_locs[1,i]
    while (curr_offset > 0 and offsets[curr_offset] > start):
      curr_offset -= 1

    while (curr_offset < num_offsets and offsets[curr_offset] < start):
      curr_offset += 1

    flagdat[curr_offset] = 1

    while (curr_offset < num_offsets and offsets[curr_offset] < end):
      flagdat[curr_offset] = 1
      curr_offset += 1


  write_grfile(offsets, flagdat, outfile)

def plot_grfile_comparison(grfile1, grfile2, outfig, title1, title2):
  """
  Plot a comparison of pairwise values from 2 gr files
  """

  import scipy.stats

  offsets1, data1 = read_grfile(grfile1)
  offsets2, data2 = read_grfile(grfile2)

  min1 = scipy.stats.scoreatpercentile(data1,2)
  max1 = scipy.stats.scoreatpercentile(data1,98)
  min2 = scipy.stats.scoreatpercentile(data2,2)
  max2 = scipy.stats.scoreatpercentile(data2,98)

  samples, xedges, yedges = pylab.histogram2d(data1, data2, bins=50, range=( (min1, max1), (min2,max2)))
  pylab.figure()
  pylab.pcolor(xedges, yedges, samples)

#  leftcornerx = xedges[0]
#  leftcornery = yedges[0]
#  rightcornerx = xedges[-1]
#  rightcornery = yedges[0] + xedges[-1] - xedges[0]
#  if (rightcornery > yedges[-1]):
#    rightcornery = yedges[-1]
#    rightcornerx = xedges[0] + (yedges[-1] - yedges[0])
#
#  pylab.plot( [leftcornerx, rightcornerx], [leftcornery, rightcornery], 'k--')

  pylab.xlabel(title1)
  pylab.ylabel(title2)

  pylab.colorbar()

  pylab.axis('tight')

  pylab.savefig(outfig)

  
def make_flags_for_strand(celtxtfile, outgrfile, strandname, genomename="E_coli_genome"):
  """
  Write a gr file with the entries from a specified strand and genome

  celtxtfile must be raw output from bpmap_bcel_join;
    its first field should be where strandname is searched for (it will be
    only one character
  """

  instr = open(celtxtfile, 'r')
  outstr = open(outgrfile, 'w')

  for line in instr:
    linearr = line.split("\t")
    if (linearr[0] == strandname and linearr[1] == genomename):
      outstr.write("%s %s\n" % (linearr[2], linearr[3]))

  instr.close()
  outstr.close()

def get_stats_on_top_percentile(grfile, flagfile, perc=99):
  """
  Count how many of the top perc-percentile entries have flags
  """
  
  import scipy.stats
  offsets, vals = read_grfile(grfile)
  flagarr = numpy.loadtxt(flagfile, dtype='int')

  flagstr = open(flagfile, 'r')
  firstline = flagstr.readline()
  linearr = firstline.split()
  flagnames = linearr[1:]
  flagstr.close()

  flagsbyprobe = scipy.zeros( (len(offsets), len(flagnames) ) )

  flagind = 0

  for j in range(len(offsets)):
    myoffset = offsets[j]

    while (flagarr[flagind, 0] < myoffset):
      flagind += 1

    if (flagarr[flagind,0] > myoffset):
      print "Warning: couldn't find offset %i" % myoffset
      continue
      
    flagsbyprobe[j,:] = flagarr[flagind,:]

  goodinds = scipy.greater(vals, scipy.stats.scoreatpercentile(vals,perc))

  print "Flag numfalse numtrue numfalse(enriched) numtrue(enriched)"
  for i in range(1,len(flagnames)):
    myflag = flagnames[i]
    numentries = len(offsets)
    numgood = len(vals[scipy.greater(flagsbyprobe[:,i], 0)])
    probeflags = flagsbyprobe[:,i]
    goodprobeflags = probeflags[goodinds]
    numenriched = len(goodprobeflags)
    numenrichedgood = len(goodprobeflags[scipy.greater(goodprobeflags, 0)])
    print "%s %i %i %i %i %6f %6f %6f %6f" % (myflag, numentries-numgood, numgood, numenriched-numenrichedgood, numenrichedgood, (numentries-numgood) / float(numentries), numgood / float(numentries), (numenriched-numenrichedgood) / float(numenriched), (numenrichedgood) / float(numenriched))
   



def print_zscore_stats_manyflags(grfile, flagfile, outprefix):
  """
  Print statistics for each flag in a given flag file

  Here we do all combinations of columns besides the first
  """

  import scipy.stats
  import numpy
  import pylab

  collist = ['0.0', '0.0', '0.55', '0.55']
  linelist = ["dashed", "solid", "dashed", "solid"]


  offsets, vals = read_grfile(grfile)
  flagarr = numpy.loadtxt(flagfile, dtype='int')

  flagstr = open(flagfile, 'r')
  firstline = flagstr.readline()
  linearr = firstline.split()
  flagnames = linearr[1:]
  flagstr.close()

  print "Names: %s" % flagnames

  # First do stats for each individual flag

  print "  %4s | %8s | %8s | %8s | %8s | %8s | %8s" % ("Flag", "N", "Mean", "Std", "Median", "Skew", "Kurtosis")
  print "=" * 80
  print " "
  print "Statistics on individual flags:"

  for i in range(1,len(flagnames)):
    print "  %s" % flagnames[i]
    truevals = []
    falsevals = []

    flagind = 0

    for j in range(len(offsets)):
      myoffset = offsets[j]

      while (flagarr[flagind, 0] < myoffset):
        flagind += 1

      if (flagarr[flagind,0] > myoffset):
        print "Warning: couldn't find offset %i" % myoffset
        continue
        
      #flagind = scipy.searchsorted(flagarr[:,0], offsets[j])
      flagval = flagarr[flagind,i]

      if (flagval > 0):
        truevals.append(vals[j])
      else:
        falsevals.append(vals[j])


    truevals = scipy.sort(truevals)
    falsevals = scipy.sort(falsevals)

    istrue = "Yes"
    for myvals in (truevals, falsevals):
      numvals = len(myvals)
      mean = scipy.mean(myvals[int(0.05*numvals):int(0.95*numvals)])
      median = numpy.median(myvals)
      std = scipy.std(myvals)
      skew = scipy.stats.skew(myvals)
      kurt = scipy.stats.kurtosis(myvals)
      print "    %4s: %8i | % 7.3f | % 7.3f | % 7.3f | % 7.3f | % 7.3f" % (istrue, numvals, mean, std, median, skew, kurt)
      istrue = "No"

    print "U-test: %s %s" % scipy.stats.mannwhitneyu(truevals, falsevals)

    myfig = pylab.figure()
    pylab.hist(truevals, bins=5000, normed=True, cumulative=True, histtype='step', label="Yes", range=(-5,5))
    pylab.hist(falsevals, bins=5000, normed=True, cumulative=True, histtype='step', label="No", range=(-5,5))
    pylab.title("Histogram for %s" % flagnames[i])
    pylab.legend(loc=2)
    pylab.xlim( (-4,4) )
    pylab.savefig("%s_%s_histogram.pdf" % (outprefix, flagnames[i]))
    print "SAVING TO %s in print_zscore_stats_manyflags A" % ("%s_%s_histogram.pdf" % (outprefix, flagnames[i]))
    pylab.xlim( (0,4) )
    pylab.ylim( (0.65,1.0) )
    pylab.savefig("%s_%s_histogram_zoom.pdf" % (outprefix, flagnames[i]))
    pylab.close(myfig)



  # Now do pairwise combinations
  for i in range(1,len(flagnames)):
    for k in range(i+1,len(flagnames)):
      print "  %s-%s" % (flagnames[i], flagnames[k])
      flagvals1 = []
      flagvals2 = []

      flagind = 0

      for j in range(len(offsets)):
        myoffset = offsets[j]

        while (flagarr[flagind, 0] < myoffset):
          flagind += 1

        if (flagarr[flagind,0] > myoffset):
          print "Warning: couldn't find offset %i" % myoffset
          continue
          
        #flagind = scipy.searchsorted(flagarr[:,0], offsets[j])
        flagvals1.append(flagarr[flagind,i])
        flagvals2.append(flagarr[flagind,k])

      set00 = []
      set01 = []
      set10 = []
      set11 = []

      for ind in range(len(flagvals1)):
        flag1 = flagvals1[ind]
        flag2 = flagvals2[ind]
        value = vals[ind]

        if (flag1 > 0):
          if (flag2 > 0):
            set11.append(value)
          else:
            set10.append(value)
        else:
          if (flag2 > 0):
            set01.append(value)
          else:
            set00.append(value)

      myfig = pylab.figure()

      stringind = 0
      #descrips = ["No,No", "No,Yes", "Yes,No", "Yes,Yes"]
      descrips = ["Noncoding, No TF", "Noncoding, TF site", "Coding, No TF", "Coding, TF site"]
      for myvals in (set00,set01,set10,set11):
        numvals = len(myvals)
        mean = scipy.mean(myvals[int(0.05*numvals):int(0.95*numvals)])
        median = numpy.median(myvals)
        std = scipy.std(myvals)
        skew = scipy.stats.skew(myvals)
        kurt = scipy.stats.kurtosis(myvals)
        print "    % 4s: %8i | % 7.3f | % 7.3f | % 7.3f | % 7.3f | % 7.3f" % (descrips[stringind], numvals, mean, std, median, skew, kurt)
        pylab.hist(myvals, bins=5000, normed=True, cumulative=True, histtype='step', range=(-5,5), lw=2, ec=collist[stringind], fc=collist[stringind], ls=linelist[stringind])
        pylab.plot( [-11,-10], [-11,-10], ls=linelist[stringind], c=collist[stringind], lw=2, label=descrips[stringind])
        stringind += 1

      pylab.title("Four way histograms for (%s,%s)" % (flagnames[i],flagnames[k]))
      pylab.legend(loc=4)
      pylab.xlim((-5,5))
      pylab.ylim( (0,1) )
      print "SAVING TO %s in print_zscore_stats_manyflags B" % ("%s_%s_histogram.pdf" % (outprefix, flagnames[i]))
      pylab.savefig("%s_%s_vs_%s_histogram.pdf" % (outprefix, flagnames[i], flagnames[k]))
      pylab.xlim((0,5.0))
      pylab.ylim( (0.70,1) )
      pylab.savefig("%s_%s_vs_%s_histogram_zoom.pdf" % (outprefix, flagnames[i], flagnames[k]))
      pylab.close(myfig)

      print "T-tests:"
      print "     %8s | %8s | %8s | %8s" % ("00", "01", "10", "11")
      print " 00  %8e | %8e | %8e | %8e" % (scipy.stats.mannwhitneyu(set00,set00)[1],scipy.stats.mannwhitneyu(set00,set01)[1], scipy.stats.mannwhitneyu(set00,set10)[1], scipy.stats.mannwhitneyu(set00,set11)[1])
      print " 01  %8e | %8e | %8e | %8e" % (scipy.stats.mannwhitneyu(set01,set00)[1],scipy.stats.mannwhitneyu(set01,set01)[1], scipy.stats.mannwhitneyu(set01,set10)[1], scipy.stats.mannwhitneyu(set01,set11)[1])
      print " 10  %8e | %8e | %8e | %8e" % (scipy.stats.mannwhitneyu(set10,set00)[1],scipy.stats.mannwhitneyu(set10,set01)[1], scipy.stats.mannwhitneyu(set10,set10)[1], scipy.stats.mannwhitneyu(set10,set11)[1])
      print " 11  %8e | %8e | %8e | %8e" % (scipy.stats.mannwhitneyu(set11,set00)[1],scipy.stats.mannwhitneyu(set11,set01)[1], scipy.stats.mannwhitneyu(set11,set10)[1], scipy.stats.mannwhitneyu(set11,set11)[1])



def print_zscore_stats(grfile, flagfile, valnamedict, retvals=True):
  """
  Print statistics for each flag in a given flag file

  We separate the values in grfile into sets based on the levels in flagfile,
    and then calculate statistics for each, and do wilcox tests for means 
    between them

  valnamedict should have the names that should be printed for each set, keyed
    by the corresponding value in flagfile

  If retvals is true, we also return a dictionary containing the statistics
    for each file
  """

  import scipy.stats

  offsets, vals = read_grfile(grfile)
  flagoffsets, flagvals = read_grfile(flagfile)

  print flagvals.dtype
  flagvals = scipy.around(flagvals)
  flagvals = scipy.array(flagvals, dtype='int')

  write_grfile(flagoffsets, flagvals, "tmpflags.gr")

  flaglevels = list(set(flagvals))

  flagnames = [valnamedict[name] for name in flaglevels]

  print "Statistics for file %s with flags %s" % (grfile, flagfile)
  print "  %10s  %8s %8s %8s %8s %8s" % ("Name", "Mean", "Std", "Median", "Skew", "Kurtosis")

  retdict = {}

  for flag,name in zip(flaglevels, flagnames):
    print "Working on flag %s, name %s" % (flag, name)
    myvals = []

    for i in range(len(offsets)):
      off = offsets[i]
      val = vals[i]
      flagoff = flagoffsets[i]
      flagval = flagvals[i]

      if (flagval == flag):
        myvals.append(val)


    #myvals = scipy.sort(vals[scipy.equal(flagvals, flag)])
    myvals = scipy.sort(scipy.array(myvals))
    #myvals = scipy.sort(vals[scipy.equal(flagvals, flag)])

    #return(myvals)

    #myvals_write = vals[scipy.equal(flagvals, flag)]
    #myinds_write = flagoffsets[scipy.equal(flagvals, flag)]

    #write_grfile(myinds_write, myvals_write, "%s_test_flag%s.gr" % (grfile, flag))

    numvals = len(myvals)

    mean = scipy.mean(myvals[int(0.05*numvals):int(0.95*numvals)])
    median = numpy.median(myvals)
    std = scipy.std(myvals)
    skew = scipy.stats.skew(myvals)
    kurt = scipy.stats.kurtosis(myvals)
    print "  %8s: % 10.6f % 7.3f % 10.6f % 7.3f % 7.3f" % (name,mean, std, median, skew, kurt)
    retdict[(name,'mean')] = mean
    retdict[(name,'median')] = median
    retdict[(name,'std')] = std
    retdict[(name,'skew')] = skew
    retdict[(name,'kurt')] = kurt

  #print "  Wilcox tests for equality of mean:"
  #for i in range(len(flaglevels)):
  #  valsi = vals[scipy.equal(flagvals, flaglevels[i])]
  #  for j in range(i+1, len(flaglevels)):
  #    valsj = vals[scipy.equal(flagvals, flaglevels[j])]
      #print "    %s - %s : %e" % (flagnames[i], flagnames[j], rpy.r.wilcox_test(valsi, valsj)['p.value'])

  if (retvals):
    return retdict
  
def shuffle_vals(invalues):
  
  import random
  random.seed()

  newvals = []
  oldvals = list(invalues)
  for i in range(len(oldvals)):
    j = random.randrange(len(oldvals))
    val = oldvals.pop(j)
    newvals.append(val)

  return scipy.array(newvals)


def calc_shuffled_peaks(infile, minlength, maxlength, percentile, smooth1length, smooth2length, outdir, outprefix, genlength=4641652, numsamples=100):
  """
  Calculate numsamples random shufflings of the input data and do peak
    calling for each of them
  """
  
  import os

  offsets, vals = read_grfile(infile)
  os.mkdir(outdir)

  for i in range(numsamples):
    shufdat = shuffle_vals(vals)
    thisprefix = "%s/%s_shuffle%i" % (outdir, outprefix, i)
    write_grfile(offsets, shufdat, "%s.gr" % thisprefix)

    for avglen in [smooth1length, smooth2length]:
      do_runningavg_opt("%s.gr" % thisprefix, "%s_runningavg_%ibp.gr" % (prefix,avglen), avglen, genlength)

    identify_peaks("%s_runningavg_21bp.gr" % thisprefix, "%s_runningavg_33bp.gr" % thisprefix, minlength, maxlength, percentile, "%s_peaks_voradef_long.gr" % thisprefix, outfile_startend = "%s_peaks_voradef_long_startend.txt")


def quantile_normalize(origmat):
  """
  Apply quantile normalization to a matrix

  Here the rows are probes and columns are samples
  XXX FIXME: NOT YET TESTED
  """

  inmat = origmat.astype(float)
  sortedmat = numpy.array(inmat, copy=True)

  nrow, ncol = inmat.shape

  sortedmat.sort(axis=0)
  meanvals = numpy.mean(sortedmat, axis=1)

  sortorder = inmat.argsort(axis=0)

  for col in range(ncol):
    sortedmat[sortorder[:,col], col] = meanvals

  return sortedmat

def normalize_by_runavg(infile, avgfile, outfile, scaledval=2500.0):
  """
  Scale a file such that the running average is uniformly scaledval

  Thus each probe is set to infile * scaledval / avgfile
  """

  locs, vals = read_grfile(infile)
  avglocs, avgvals = read_grfile(avgfile)

  newvals = vals * scaledval / avgvals
  write_grfile(locs, newvals, outfile)

def sorted_list_median(mylist):
  # Helper function -- gives the median of a sorted list
    if (len(mylist) % 2 == 0):
      middle_ind = len(mylist) / 2
      return (mylist[middle_ind] + mylist[middle_ind - 1]) / 2.0
    else:
      middle_ind = len(mylist) / 2 
      return mylist[middle_ind]

def score_normed_vals(scorefile, meanstdref, outfile):
  """
  Apply 09032010 scoring to scorefile

  Each score is log(val) - mean(log(refs)) / sigma(log(refs))
  meanstdref must contain columns with index, mean(log(refs)), and sigma(log(refs))
  Note that the input data to be scored is **not** expected to be log transformed--
    we do that inside of this function
  """

  inds, vals = read_grfile(scorefile)
  refinds, means, sigmas = numpy.loadtxt(meanstdref, unpack=True)
  # Make sure we don't get nans -- replace zeroes with the lowest nonzero value
  #smallestval = numpy.min(vals[vals > 0.001])
  #vals[vals < smallestval] = smallestval
  zscores = (numpy.log2(vals) - means) / sigmas
  write_grfile(inds, zscores, outfile)

def score_normed_vals_nolog(scorefile, meanstdref, outfile):
  """
  Apply 09032010 scoring to scorefile

  Each score is val - mean(refs) / sigma(refs)
  meanstdref must contain columns with index, mean(refs), std(refs)
  Note that the input data to be scored is **not** expected to be log transformed--
    we do that inside of this function
  """

  inds, vals = read_grfile(scorefile)
  refinds, means, sigmas = numpy.loadtxt(meanstdref, unpack=True)
  zscores = (vals - means) / sigmas
  write_grfile(inds, zscores, outfile)


def apply_index_offset(infile, outfile, offset):
  """
  Apply a uniform offset to all indices in infile and write the new data to outfile
  """

  inds, vals = read_grfile(infile)
  write_grfile(inds+offset, vals, outfile)

def spline_correct_genome(infile,outfile,plot=False,genome_length = 4641652, oriCloc=3923883,perspline=True):
  # do a spline-based correction of any periodicity in the input file

  temp_smoothed = tempfile.NamedTemporaryFile(delete=False)
  temp_smoothed.close()

  spline_smooth_genome( infile, temp_smoothed.name, plot=plot,genome_length = genome_length, oriCloc=oriCloc,perspline=perspline)
  divide_gr_files( infile,temp_smoothed.name, outfile)
  os.remove(temp_smoothed.name)

def spline_correct_genome_fromref(infile,reffile,outfile,plot=False,genome_length = 4641652, oriCloc=3923883):
  # do a spline-based correction of any periodicity in the input file
  # we fit the spline based on reffile, but then smooth the data from outfile
  # typically the "input" sample will be used as reffile

  temp_smoothed = tempfile.NamedTemporaryFile(delete=False)
  temp_smoothed.close()

  spline_smooth_genome( reffile, temp_smoothed.name, plot=plot,genome_length = genome_length, oriCloc=oriCloc)
  divide_gr_files( infile,temp_smoothed.name, outfile)
  os.remove(temp_smoothed.name)


def spline_smooth_genome(infile, outfile, plot=True, plotfile="spline.pdf", genome_length = 4641652, oriCloc=3923883, outfile_all = None, perspline=True):
  """
  Write a spline-smoothed version of the data from infile to outfile

  We do very heavy smoothing; a periodic spline is fitted with four knots, at oriC and quadrants associated with it

  Other than that we use splrep/splenv with standard settings

  If outfile_all is not None, we write a file containing the value
    of the interpolation at all integer positions

  If perspline is false, we do not make the spline periodic
  """

  offs,dat = read_grfile(infile)

  # augment offs/dat with a data point ensuring proper periodicity
  offs_recentered = (offs - oriCloc) % genome_length
  sortorder = numpy.argsort(offs_recentered)

  offs_mod = offs_recentered[sortorder]
  vals_mod = dat[sortorder]
  offs_orig_mod = offs[sortorder]
  resort_order = numpy.argsort(offs_orig_mod)

  offs_aug = numpy.append(offs_mod, genome_length+offs_mod[0])
  dat_aug = numpy.append(vals_mod, vals_mod[0])
  #offs_aug = numpy.append(offs, genome_length+offs[0])
  #dat_aug = numpy.append(dat, dat[0])

  knots=[genome_length / 4, genome_length/2, 3 * genome_length / 4]

  print knots
  print len(offs_aug)
  print len(dat_aug)
  allspl = si.splrep(offs_aug, dat_aug, k = 3, t=knots, task=-1, w=scipy.ones(len(offs_aug)), per=perspline, full_output=1)
  myspl = allspl[0]
  print myspl
  print allspl[1:]
  if allspl[2] > 0:
    raise(ValueError("ERROR FITTING SPLINE"))
  print "----"
  splvals = (si.splev(offs_mod, myspl))[resort_order]

  #print si.splev(offs[0], myspl)
  #print si.splev(offs[-2], myspl)
  #print si.splev(offs[0], myspl, 1)
  #print si.splev(offs[-2], myspl, 1)

  write_grfile(offs, splvals, outfile)

  if outfile_all is not None:
    raise("Fix this function to account for altered periodicity")
    alloffs = scipy.arange(genome_length)
    write_grfile(alloffs, si.splev(alloffs, myspl), outfile_all)

  if plot:
    pylab.figure()
    pylab.plot(offs, dat, 'b-')
    pylab.plot(offs, splvals, 'r--')
    knots = myspl[0]
    knots = knots[(knots >= 0) & (knots <= genome_length)]
    knotvals = si.splev(knots, myspl)
    print knots
    pylab.plot((knots+oriCloc)%genome_length, knotvals, 'yo')
    pylab.plot(oriCloc, 0, 'go')
    pylab.savefig(plotfile)

## The following functions make and manipulate copy index files
# The copy index is defined as follows:
# For each probe, approximate the matching DNA : probe ratio for that probe
# This is matchseqs \cdot normamounts / numprobes
# Matchseqs is a vector containing 1 at all positions where the sequence matches the probe, 0 elsewhere
# normamounts is the normalized amount of DNA (based on spline fitting from ipod_utils) present at each position
# Numprobes is the number of probes on the array that match the sequence (but this time **not** including its reverse complement)
# This program expects indexing of the probes and normamounts to agree--they should all either be numbered based
# on the start or the middle of the probe


def make_ci_file(allsplinefile, outfile):
  """
  Make a gr file containing the copy index as described above

  We must input a spline file containing *all* genome locations, not just
    those of probes
  """

  probeseqfile="/home/petefred/ST_research/ipod/tvora_orig_data/E_coli_matchseqs_file.gr"
  probehashfile = "/home/petefred/ST_research/ipod/tvora_orig_data/E_coli_probe_hash_file.txt"
  affyfile = "/home/petefred/ST_research/ipod/tvora_orig_data/probe_design_affy_E_coli.txt"

  probestr = open(affyfile)

  indices, normamounts = read_grfile(allsplinefile)

  probestarts = []
  probeseqs = []
  for line in probestr:
    if line[0] == "#":
      continue

    linearr = line.split()
    probestarts.append(int(linearr[3]))
    probeseqs.append(linearr[0])

  probe_count_hash = pickle.load(open(probehashfile))
  probe_match_seqs_hash = pickle.load(open(probeseqfile))

  probe_cis = []

  #probe_count_hash = pickle.load(match_seqs_file)
  #all_matchseqs = numpy.loadtxt(match_seqs_file)


  for i in range(len(probestarts)):
    if (i % 100 == 0):
      print "Working on probe %i" % i

    matchseq_inds = probe_match_seqs_hash[i]
    num_probe_matches = probe_count_hash[probeseqs[i]]

    dotprod = 0.0
    for ind in matchseq_inds:
      dotprod += normamounts[ind]

    my_probe_ci = dotprod / num_probe_matches

    probe_cis.append(my_probe_ci)

  write_grfile(probestarts, probe_cis / numpy.min(probe_cis), outfile)

def make_matchseqs_list(myseq, genomeseq):
  """
  Return a list showing all indices of positions whose sequence matches this one

  For these purposes we ignore which strand matches
  We do include the trivial identity match
  n.b. we look over the entire genome, not just at probe positions
  """

  myseq_len = len(myseq)
  myseq_rc = reverse_complement(myseq)
  retlist = []

  for seq in [myseq,myseq_rc]:
    startind = 0
    while (startind >= 0):
      startind = genomeseq.seq.find(seq, start=startind)
      if (startind >= 0):
        retlist.append(startind)
        startind += 1

  return retlist


def make_numprobes_hash(probeseqs):
  # Make a dictionary containing the number of probes with identical sequences, indexed by sequence

  probehash = {}

  for ps in probeseqs:
    if probehash.has_key(ps):
      probehash[ps] = probehash[ps] + 1
    else:
      probehash[ps] = 1


  return probehash

def prepare_ref_files(affy_file, genome_fasta, probe_hash_file, matchseqsfile):
  """
  Write the constant files needed to produce copy indices

  This function reads the affy probes and genome sequence, and writes the
    probe count hash file and the matchseq vectors for all probes

  Neither of these things is hyb dependent; only the spline amounts vary
  """

  genome = SeqIO.read(open(genome_fasta), "fasta")

  probestr = open(affy_file)

  probestarts = []
  probeseqs = []
  for line in probestr:
    if line[0] == "#":
      continue

    linearr = line.split()
    probestarts.append(int(linearr[3]))
    probeseqs.append(linearr[0])

  probe_count_hash = make_numprobes_hash(probeseqs)

  outstr = open(probe_hash_file, 'wb')
  pickle.dump(probe_count_hash, outstr)
  outstr.close()

  all_matchseqs = {}

  for i in range(len(probestarts)):
    print "Working on probe %i" % i
    matchseqs = make_matchseqs_list(probeseqs[i], genome)
    all_matchseqs[i] = matchseqs

  outstr = open(matchseqsfile, 'wb')
  pickle.dump(all_matchseqs, outstr)
  outstr.close()

def do_standard_analysis_temp_permute(prefix, pathdict):
  """
  Sandbox function for permuting between tjv and 20100903 analysis

  -first: switch from pm-mm to pm only
  """

  GENOME_LENGTH = 4641652
  REFS_MEANSTD = pathdict['REFS_MEANSTD']
  REFS_MEANSTD_20100903 = pathdict['REFS_MEANSTD_20100903']
  REFS_80KB_20100903 = pathdict['REFS_80KB_20100903']
  REFS_MEANSTD_LOG_20100903 = pathdict['REFS_MEANSTD_LOG_20100903']
  REFS_80KB = pathdict['REFS_80KB']
  BPMAPFILE = pathdict['BPMAPFILE']
  CODINGFLAGFILE = pathdict['CODINGFLAGFILE']

  celfile = "%s.CEL" % prefix
  txtfile = "%s_sandbox.txt" % prefix
  cleantxtfile = "%s_cleaned_sandbox.txt" % prefix
  grfile = "%s_ecolionly_sandbox.gr" % prefix
  grfile_sc2500 = "%s_ecolionly_sc2500_sandbox.gr" % prefix
  gr80file = "%s_ecolionly_sc2500_avg80kb_sandbox.gr" % prefix
  gr80splinefile = "%s_ecolionly_sc2500_avg80kb_spline_sandbox.gr" % prefix
  grfile_localnorm = "%s_localnorm_sandbox.gr" % prefix
  normdat_raw = "%s_ecolionly_normalized_sandbox.gr" % prefix
  normdat_zscore = "%s_standard-analysis_ecolionly_zscores_sandbox.gr" % prefix
  normdat_zscore_log = "%s_standard-analysis_ecolionly_logspace_zscores_sandbox.gr" % prefix

  convert_cel_to_text(celfile, BPMAPFILE, txtfile)
  clean_text(txtfile, cleantxtfile)
  extract_genome_flag(cleantxtfile, grfile)
  offset_median_value(grfile, grfile_sc2500, 2500.0)
  do_runningavg_opt(grfile, gr80file, width=80001, genomelength=GENOME_LENGTH)
  spline_smooth_genome(gr80file, gr80splinefile, plot=False)
  #make_ci_file(splinefile_all, ci_file)
  #normalize_by_runavg(grfile, gr80splinefile, grfile_localnorm, 2500.0)
  normalize_vs_ref(grfile_sc2500, gr80file, REFS_MEANSTD_20100903, REFS_80KB, normdat_raw, normdat_zscore)
  #normalize_by_runavg(grfile, gr80file, grfile_localnorm, 2500.0)
  #normalize_median_value(grfile_localnorm, normdat_raw, 2500.0)
  #score_normed_vals_nolog(normdat_raw, REFS_MEANSTD_20100903,  normdat_zscore)
  #apply_index_offset(normdat_zscore, normdat_zscore_offset, CENTER_OFFSET)
  #score_normed_vals(normdat_raw, REFS_MEANSTD_LOG_20100903,  normdat_zscore_log)
  #apply_index_offset(normdat_zscore_log, normdat_zscore_offset_log, CENTER_OFFSET)

  #convert_cel_to_text(celfile, BPMAPFILE, txtfile)
  #clean_and_corretext(txtfile, cleantxtfile)
  #extract_genome_flag(cleantxtfile, grfile)
  #offset_median_value(grfile, grfile_sc2500, 2500.0)
  #do_runningavg_opt(grfile_sc2500, gr80file, width=80001, genomelength=GENOME_LENGTH)
  #normalize_vs_ref(grfile_sc2500, gr80file, REFS_MEANSTD, REFS_80KB, normdat_raw, normdat_zscore)
  #do_runningavg_opt(normdat_zscore, "%s_runavg128.gr" % smoothed_prefix, width=129, genomelength=GENOME_LENGTH)
  #plot_zscore_hist("%s_runavg128.gr" % smoothed_prefix, hist_fig_file, CODINGFLAGFILE, zerolabel="Coding", onelabel="Noncoding")
  #plot_zscore_cum_hist("%s_runavg128.gr" % smoothed_prefix, chist_fig_file, CODINGFLAGFILE, zerolabel="Coding", onelabel="Noncoding")
  #convolve_file(normdat_zscore, gaussed_zscore_file, genomelength=GENOME_LENGTH)

def variance_normalize_file(infile, outfile):
  """
  Variance normalize the scores in a gr file
  """

  offs, dat = read_grfile(infile)
  centerdat = dat - scipy.mean(dat)
  vardat = centerdat / scipy.std(dat)

  write_grfile(offs, vardat, outfile)

def gaussian_normalize_file(infile,outfile,target_median, target_iqr):
  """
  Normalize a file based on its median and iqr

  The median is set to target_median, and then the data scaled to have
    an iqr of target_iqr
  """

  offs, dat = read_grfile(infile)
  centerdat = dat - scipy.median(dat) + target_median
  iqr = scipy.stats.scoreatpercentile(centerdat, 75) - scipy.stats.scoreatpercentile(centerdat, 25)
  scaledat = centerdat * target_iqr / iqr

  write_grfile(offs, scaledat, outfile)


def find_min_max_ind(values):
  """
  Find and return the indices of the min and max values in an array
  """

  minind = 0
  minval = values[0]
  maxind = 0
  maxval = values[0]

  for i in range(1, len(values)):
    if (values[i]) > maxval:
      maxval = values[i]
      maxind = i

    if (values[i] < minval):
      minval = values[i]
      minind = i

  return (minind, maxind)

def find_index_subset(infile, targetfile, outfile, skiprows=0):
  # Given 3 gr files, write only lines from infile to outfile matching indices from targetfile
  locs, vals = read_grfile(infile)
  targetlocs, targetvals = read_grfile(targetfile, skiprows=skiprows)
  targetloc_set = set(targetlocs)
  newlocs = []
  newvals = []

  for i in range(len(locs)):
    if locs[i] in targetloc_set:
      newlocs.append(locs[i])
      newvals.append(vals[i])

  write_grfile(newlocs, newvals, outfile)

def pc_subtract_counts(file1,file2,outfile,offset=0.0):
  # given two gr files, log2 scale them (after adding an offset to the values), and then subtract values from
  #  the second file from those in the first file, after rescaling according to the ratios of their contributions to
  #  the first principle component
  #  see my notes from 10/16/17 for the justification for this

  locs1,vals1=read_grfile(file1)
  locs2,vals2=read_grfile(file2)

  log_val1 = numpy.log2(vals1 + offset)
  log_val2 = numpy.log2(vals2 + offset)

  X = numpy.vstack( (log_val1,log_val2) )
  pca=sklearn.decomposition.PCA(n_components=2,whiten=False) 
  pca.fit(X.transpose())
  
  c_mat = pca.components_
  scalefac = c_mat[0,1] / c_mat[0,0]

  rescaled_1 = 2**(log_val1 * scalefac)
  rescaled_2 = 2**log_val2

  newvals = rescaled_1 - rescaled_2

  write_grfile(locs1,newvals,outfile)


def loess_subtract_counts(file1,file2,outfile,outplot=None,lscalefac=1.0):
  # given two gr files, do a loess fit to find the best mapping of file2 to file1, and then subtract from file1
  #  the fitted value at each point
  # write the difference to outfile
  # we multiple the loess-normalized counts by lscalefac prior to subtraction to account for some spread

  locs1,vals1=read_grfile(file1)
  locs2,vals2=read_grfile(file2)

  allinds = range(len(locs1))
  randinds = range(0,len(locs1),len(locs1)/500) + [len(locs1)-1]

  randinds.sort()
  ordered_inds = numpy.array(randinds)


  sortord = numpy.argsort(vals2)

  v1_forloess = vals1[sortord][ordered_inds]
  v2_forloess = vals2[sortord][ordered_inds]

  growthcurves.do_loess_fit(v2_forloess,v1_forloess,'myfit')
  xmin = numpy.min(vals2)
  xmax = numpy.max(vals2)
  xvals_test_plot = numpy.linspace(xmin,xmax)
  print "done with fit"
  predvals = numpy.array(growthcurves.do_loess_interpolation_array(list(vals2),'myfit') )
  predvals_forplot = growthcurves.do_loess_interpolation_array(list(xvals_test_plot),'myfit')
  print predvals[:10]
  pylab.figure()
  pylab.bone()
  #pylab.hexbin(vals2,vals1,gridsize=50,xscale='log',yscale='log')
  pylab.plot(vals2,vals1,'bo')
  pylab.plot(xvals_test_plot,predvals_forplot,'g-')
  pylab.plot(xvals_test_plot,2*numpy.array(predvals_forplot),'b-')
  pylab.plot(xvals_test_plot,0.5*numpy.array(predvals_forplot),'r-')
  pylab.plot(xvals_test_plot,lscalefac*numpy.array(predvals_forplot),'y-')
  pylab.xlabel(file2)
  pylab.ylabel(file1)
  if (outplot is not None):
    pylab.savefig(outplot + "_2dhist.png")

  # FIX THIS
  #pylab.figure()
  #pylab.hist(numpy.log2( v1_forloess / v2_forloess), bins=25)
  #if (outplot is not None):
  #  pylab.savefig(outplot + "_1dhist.png")

  #print "A"
  #pylab.figure()
  #tmpvals = vals1 / predvals
  #print "B"
  #pylab.hist(numpy.log2(tmpvals[numpy.isfinite(tmpvals)]),bins=25)
  #if (outplot is not None):
  #  pylab.savefig(outplot + "_1dhist_selfsub.png")
#
#  print "C"

  

  print type(vals1)
  print type(predvals)
  print type(vals1[0])
  print type(predvals[0])
  print vals1.shape
  print predvals.shape
  predvals.shape = vals1.shape

  newvals = vals1 - lscalefac*numpy.maximum(numpy.zeros_like(predvals),predvals)
  write_grfile(locs1,newvals,outfile)

def slope_subtract_counts_log(file1, file2, outfile,offset=0):
  # find the linear fit that minimizes the mean distance (rather than msd) of the line to a scatter plot of file2 vs file1, in log2 space
  # then transform the data from file1 according to that fit, and subtract the file2 data from the rescaled file1 data IN THE ORIGINAL LINEAR SPACE
  # we write that result to outfile
  # offset is added to all values before log scaling

  locs1,vals1=read_grfile(file1)
  locs2,vals2=read_grfile(file2)

  log_file_1 = numpy.log2(vals1 + offset)
  log_file_2 = numpy.log2(vals2 + offset)

  # now I search a bunch of slope/intercept combinations to find the best one
  bestints = []
  fitscores = []
  slopes = numpy.arange(0,5,0.01)
  #slopes = numpy.arange(0,5,0.1)


  for slope in slopes:
    print "working on slope %f" % slope
    bestscore = 0.0
    bestint = 0.0

    xterm = slope * log_file_1 - log_file_2
    denom = numpy.sqrt(slope*slope + 1)

    #print xterm
    #print denom

    for yint in numpy.arange(-50,50,0.01):
    #for yint in numpy.arange(-50,50,0.1):
      #print "working on int %f" % yint
      numerators = numpy.abs( xterm + yint)
      distmed = numpy.mean( numerators / denom )
      #distsum = numpy.sum(log_chip_vals - slope * log_ipod_vals + yint)
      #distsum = 0
      #for i in range(len(log_ipod_vals)):
      #    distsum += d_point_line( log_ipod_vals[i], log_chip_vals[i], slope, yint)

      score = 1.0/distmed
      #print score

      if score > bestscore:
         bestscore = score
         bestint = yint

    fitscores.append(bestscore)
    bestints.append(bestint)

  # now renormalize and process the data
  bestind = numpy.argmax(fitscores)

  bestslope = slopes[bestind]
  bestint = bestints[bestind]
  
  rescaled_f1 = numpy.power(2, log_file_1 * bestslope + bestint)
  new_data = rescaled_f1 - vals2
  write_grfile(locs1,new_data,outfile)

 

def lin_subtract_counts_log(file1,file2,outfile,outplot=None,soffset=0.0):
  # given two gr files, do a linear fit to find the best mapping of file2 to file1, and then subtract from file1 (note that we only fit things with an rna polymerase value over 0)
  #  the fitted value at each point
  # write the difference to outfile
  # soffset is added to all of the fitted values prior to subtraction
  # This is basically designed for finding differences between log scaled data sets (e.g., ipod vs inp log ratio files)

  print "entering lin_subtract_counts_log with %s %s %s %s %s" % (file1,file2,outfile,outplot,soffset)

  locs1,vals1=read_grfile(file1)
  locs2,vals2=read_grfile(file2)

  sortord = numpy.argsort(vals2)

  v1_forfit = vals1[sortord]
  v2_forfit = vals2[sortord]

  # apply a correction to prevent negative rnapol values to correspond to positive expected ipod occupancies
  goodflags = (v2_forfit > soffset)
  v1_forfit = v1_forfit[goodflags]
  v2_forfit = v2_forfit[goodflags]


  #myspl = growthcurves.fit_bspline(v2_forfit,v1_forfit,numknots=5)
  slope,intercept,rval,pval,stderr = scipy.stats.linregress( v2_forfit,v1_forfit)
  zeroval = -1 * intercept / slope
  #growthcurves.do_loess_fit(v2_forfit,v1_forfit,'myfit')

  #def get_predval(xval):
  #  if xval < zeroval:
  #    return 0.0
  #  else:
  #    return slope * xval + intercept

  xmin = numpy.min(vals2)
  xmax = numpy.max(vals2)
  xvals_test_plot = numpy.linspace(xmin,xmax)
  print "done with fit"
  #predvals = numpy.array( [ get_predval(x) for x in vals2 ] )
  #predvals_forplot = numpy.array( [ get_predval(x) for x in xvals_test_plot] )
  predvals = numpy.maximum(slope * vals2 + intercept, numpy.zeros_like(vals2))
  predvals_forplot = numpy.maximum(slope * xvals_test_plot + intercept, numpy.zeros_like(xvals_test_plot))
  #predvals = growthcurves.eval_bspline(vals2,myspl)
  #predvals_forplot = growthcurves.eval_bspline(xvals_test_plot,myspl)

  print xvals_test_plot
  pylab.figure()
  pylab.bone()
  pylab.hexbin(vals2,vals1,gridsize=50,xscale='linear',yscale='linear',bins='log')
  #pylab.plot(vals2,vals1,'bo')
  pylab.plot(xvals_test_plot,predvals_forplot,'g-',linewidth=2)
  #pylab.plot(xvals_test_plot,1 + predvals_forplot,'b-')
  #pylab.plot(xvals_test_plot,predvals_forplot - 1,'r-')
  pylab.plot(xvals_test_plot,soffset+predvals_forplot,'y-',linewidth=2)
  pylab.xlabel(file2)
  pylab.ylabel(file1)
  if (outplot is not None):
    pylab.savefig(outplot + "_2dhist.png")

  #pylab.figure()
  #pylab.hist(numpy.log2( v1_forloess / v2_forloess), bins=25)
  #if (outplot is not None):
  #  pylab.savefig(outplot + "_1dhist.png")

  print "A"
  pylab.figure()
  tmpvals = vals1 - (predvals+soffset)
  print "B"
  pylab.hist(tmpvals, bins=25)
  if (outplot is not None):
    pylab.savefig(outplot + "_1dhist_subvals.png")

  print "C"

  print vals1.shape
  print predvals.shape
  predvals.shape = vals1.shape

  newvals = vals1 - (predvals + soffset)
  write_grfile(locs1,newvals,outfile)

def log_subtract_counts_linfit(file1,file2,outfile,soffset=0.0):
  # given two gr files, do a linear fit to find the best mapping of file2 to file1, and then put 2 to the power of both of them and take the log2 of the difference
  # write the log2 difference to outfile


  locs1,vals1=read_grfile(file1)
  locs2,vals2=read_grfile(file2)

  #myspl = growthcurves.fit_bspline(v2_forfit,v1_forfit,numknots=5)
  slope,intercept,rval,pval,stderr = scipy.stats.linregress( vals2, vals1)
  predvals_2 = vals2 * slope + intercept + soffset

  subvals = (2**vals1) - (2**predvals_2)
  logsubvals = numpy.log2(subvals)

  write_grfile(locs1,logsubvals,outfile)

def subtract_counts_linfit(file1,file2,outfile,soffset=0.0):
  # given two gr files, do a linear fit to find the best mapping of file2 to file1, and then write the difference


  locs1,vals1=read_grfile(file1)
  locs2,vals2=read_grfile(file2)

  #myspl = growthcurves.fit_bspline(v2_forfit,v1_forfit,numknots=5)
  slope,intercept,rval,pval,stderr = scipy.stats.linregress( vals2, vals1)
  predvals_2 = vals2 * slope + intercept + soffset

  subvals = vals1 - predvals_2

  write_grfile(locs1,subvals,outfile)

def spline_subtract_counts(file1,file2,outfile,outplot=None,soffset=0.0):
  # given two gr files, do a spline fit to find the best mapping of file2 to file1, and then subtract from file1
  #  the fitted value at each point
  # write the difference to outfile
  # soffset is added to all of the fitted values prior to subtraction
  # This is basically designed for finding differences between log scaled data sets (e.g., ipod vs inp log ratio files)

  print "entering spline_subtract_counts with %s %s %s %s %s" % (file1,file2,outfile,outplot,soffset)

  locs1,vals1=read_grfile(file1)
  locs2,vals2=read_grfile(file2)

  sortord = numpy.argsort(vals2)

  v1_forfit = vals1[sortord]
  v2_forfit = vals2[sortord]

  # apply a correction to prevent negative rnapol values to correspond to positive expected ipod occupancies
  v1_forfit[v2_forfit < 0.0] = 0.0


  myspl = growthcurves.fit_bspline(v2_forfit,v1_forfit,numknots=5)

  xmin = numpy.min(vals2)
  xmax = numpy.max(vals2)
  xvals_test_plot = numpy.linspace(xmin,xmax)
  print "done with fit"
  predvals = growthcurves.eval_bspline(vals2,myspl)
  predvals_forplot = growthcurves.eval_bspline(xvals_test_plot,myspl)

  print xvals_test_plot
  pylab.figure()
  pylab.bone()
  pylab.hexbin(vals2,vals1,gridsize=50,xscale='linear',yscale='linear',bins='log')
  #pylab.plot(vals2,vals1,'bo')
  pylab.plot(xvals_test_plot,predvals_forplot,'g-')
  pylab.plot(xvals_test_plot,1 + predvals_forplot,'b-')
  pylab.plot(xvals_test_plot,predvals_forplot - 1,'r-')
  pylab.plot(xvals_test_plot,soffset+predvals_forplot,'y-')
  pylab.xlabel(file2)
  pylab.ylabel(file1)
  if (outplot is not None):
    pylab.savefig(outplot + "_2dhist.png")

  #pylab.figure()
  #pylab.hist(numpy.log2( v1_forloess / v2_forloess), bins=25)
  #if (outplot is not None):
  #  pylab.savefig(outplot + "_1dhist.png")

  print "A"
  pylab.figure()
  tmpvals = vals1 - (predvals+soffset)
  print "B"
  pylab.hist(tmpvals, bins=25)
  if (outplot is not None):
    pylab.savefig(outplot + "_1dhist_subvals.png")

  print "C"

  print vals1.shape
  print predvals.shape
  predvals.shape = vals1.shape

  newvals = vals1 - (predvals + soffset)
  write_grfile(locs1,newvals,outfile)

def clamp_grfile(infile,outfile,minval=0.0):
  # set to minval all values below it
  
  locs,vals=read_grfile(infile)
  vals[vals<minval] = minval
  write_grfile(locs,vals,outfile)
  
def calc_rzscores(infile,outfile):
    # calculate robust z-scores for a given input file
    # here we calculate the rzscore as the number of mads away from the median each data point
    #  we use scaling of the mad by 1.4826 to use it as an estimator for the standard deviation


    locs,vals=read_grfile(infile)
    mad = 1.4826 * numpy.median( numpy.abs( vals - numpy.median(vals) ) )
    zscores = (vals - numpy.median(vals)) / mad
    write_grfile(locs, zscores, outfile)
