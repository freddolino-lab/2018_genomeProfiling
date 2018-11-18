#!/usr/bin/python

import numpy
import scipy
import pylab
import scipy.interpolate as si
import growthcurves

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

def divide_gr_files(grfile1, grfile2, outgrfile):
  """
  Write the difference between two gr files to a third
  """

  offsets, data1 = read_grfile(grfile1)
  offsets, data2 = read_grfile(grfile2)

  data = data1 / data2

  write_grfile(offsets, data, outgrfile)

def read_grfile(filename, skiprows=0):
  """
  Read data from a .gr file into a numpy array

  In the process we transpose the array, so that it has a[0] the indices and a[1] the values
  """

  locs = numpy.loadtxt(filename, dtype='int', usecols=(0,), skiprows=skiprows)
  vals = numpy.loadtxt(filename, usecols=(1,), skiprows=skiprows)
  return (locs,vals)

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


