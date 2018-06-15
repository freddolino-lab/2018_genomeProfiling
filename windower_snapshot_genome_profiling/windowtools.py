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
import bisect
import numpy as np
import gfftools
import random


def gff_median_center(gff_entry):
    """
    Given a gff entry, find the median location between start and end 
    coordinates

    Returns a tuple of (median, median)
    """


    median = int(np.median([gff_entry.start, gff_entry.end]))
    return (median, median)

def gff_threeprime_center(gff_entry):
    """
    Given a gff entry, find the most 3' locations. 

    Returns a tuple of (3'end, 3'end)
    i.e.
    If the strand is - then
    returns a tuple (gff.start, gff.start)
    else
    returns a tuple (gff.end, gff.end)


    """
    if gff_entry.direction == "-":
        return(gff_entry.start, gff_entry.start)
    else:
        return(gff_entry.end, gff_entry.end)

def gff_fiveprime_center(gff_entry):
    """
    Given a gff entry, find the most 5' locations. 

    Returns a tuple of (5'end, 5'end)
    i.e.
    If the strand is - then
    returns the (gff.end, gff.end)
    else
    returns the (gff.start, gff.start)
    """
    if gff_entry.direction == "-":
        return(gff_entry.end, gff_entry.end)
    else:
        return(gff_entry.start, gff_entry.start)

def gff_identity_center(gff_entry):
    """
    Given a gff entry, return the window the gff spans

    Returns (gff.start, gff.end)
    """
    return(gff_entry.start, gff_entry.end)


def get_gff_window(gff_entry, fiveprime, threeprime, center_func):
    """
    Given a gff entry, return a padded window around the "center". Takes
    into account the strand of the gff entry and will add padding to the 5'
    or 3' direction relative to the strand.

    Inputs:
        gff_entry - a GffEntry object from gfftools
        fiveprime - int number of bp to pad in the fiveprime direction
        threeprime - int number of bp to pad in the threeprime direction
        center_func - function that takes a gff_entry and returns a window tuple
                      of (start, end) Where start < end.
    Returns:
        a tuple (start, end) containing the new window.
    """

    start, end = center_func(gff_entry)
    if not gff_entry.direction == "-":
        return (start-fiveprime, end+threeprime)
    else:
        return(start-threeprime, end+fiveprime)

def gff_median_center(gff_entry):
    """
    Given a gff entry, find the median location between start and end 
    coordinates

    Returns a tuple of (median, median)
    """


    median = int(np.median([gff_entry.start, gff_entry.end]))
    return (median, median)

def gff_threeprime_center(gff_entry):
    """
    Given a gff entry, find the most 3' locations. 

    Returns a tuple of (3'end, 3'end)
    i.e.
    If the strand is - then
    returns a tuple (gff.start, gff.start)
    else
    returns a tuple (gff.end, gff.end)


    """
    if gff_entry.direction == "-":
        return(gff_entry.start, gff_entry.start)
    else:
        return(gff_entry.end, gff_entry.end)

def gff_fiveprime_center(gff_entry):
    """
    Given a gff entry, find the most 5' locations. 

    Returns a tuple of (5'end, 5'end)
    i.e.
    If the strand is - then
    returns the (gff.end, gff.end)
    else
    returns the (gff.start, gff.start)
    """
    if gff_entry.direction == "-":
        return(gff_entry.end, gff_entry.end)
    else:
        return(gff_entry.start, gff_entry.start)

def gff_identity_center(gff_entry):
    """
    Given a gff entry, return the window the gff spans

    Returns (gff.start, gff.end)
    """
    return(gff_entry.start, gff_entry.end)


def get_gff_window(gff_entry, fiveprime, threeprime, center_func):
    """
    Given a gff entry, return a padded window around the "center". Takes
    into account the strand of the gff entry and will add padding to the 5'
    or 3' direction relative to the strand.

    Inputs:
        gff_entry - a GffEntry object from gfftools
        fiveprime - int number of bp to pad in the fiveprime direction
        threeprime - int number of bp to pad in the threeprime direction
        center_func - function that takes a gff_entry and returns a window tuple
                      of (start, end) Where start < end.
    Returns:
        a tuple (start, end) containing the new window.
    """

    start, end = center_func(gff_entry)
    if not gff_entry.direction == "-":
        return (start-fiveprime, end+threeprime)
    else:
        return(start-threeprime, end+fiveprime)

def discretize_window(window, strand, bins):
    """ Given a [start, end] window (assuming 1 based coordinates with inclusive
    ends) break the window into bins and return them in the 5' to 3' direction.

    Inputs:
        window - tuple(start, end)
        strand - strand associated with the window
        bins - number of bins to break the window into
    Returns:
        a list of bins # of windows in the 5' to 3' direction according to
        strand

    >>> discretize_window((1, 10), "+", 2)
    [(1, 5), (6, 10)]
    >>> discretize_window((1, 10), "-", 2)
    [(6, 10), (1, 5)]
    >>> discretize_window((1, 9), "+", 2)
    [(1, 4), (5, 9)]
    >>> discretize_window((1, 3), "+", 3)
    [(1, 1), (2, 2), (3, 3)]
    >>> discretize_window((1, 11), "+", 2)
    [(1, 5), (6, 11)]
    >>> discretize_window((1, 3), "+", 4)
    Traceback (most recent call last):
        ...
    ValueError: bins are larger than size of window. Window (1, 3), bins 4, windowlength 3
    """
    outwindows = []
    # since we are in 1 based coordinates we need to add 1 to the length
    window_length = window[1]+1 - window[0]
    if bins > window_length:
        raise ValueError("bins are larger than size of window. Window %s, bins %s, windowlength %s"%(window, bins, window_length))

    # using integer math here which floors the division
    binsize = window_length/bins
    start = window[0]
    # loop through more bins than needed
    for i in xrange(0, bins+1):
        # since theses are inclusive indicies, don't include the last bp
        outwindows.append([start, start+binsize-1])
        # increment the start
        start = start+binsize
    outwindows = outwindows[0:-1]
    outwindows[-1][-1] = window[1]
    # turn them into tuples
    outwindows = [tuple(this_window) for this_window in outwindows]
    # warn if we don't make enough bins, this should never happen with the
    # exception I am raising before
    if len(outwindows) != bins:
        print "warning not enough windows made: %s %s %s %s"%(window, strand,bins, len(outwindows))
        print "%s\n%s"%(window, outwindows)
    # reverse the windows if on the minus strand to keep them in the 3' to 5'
    # direction
    if not strand == "-":
        return outwindows
    else:
        return outwindows[::-1]

def add_window(window, size, strand, fiveprime=True):
    """ Creates a window adjacent to a window on either the fiveprime or
    three prime end of an existing window

    Inputs:
        window - tuple(start, end) window to put window adjacent
        size - int size in bp of new window to create
        strand - strand associated with the input window
        fiveprime - bool add new window to the five prime side of the input
        window

    Returns:
        (start, stop) of newly created window

    >>> add_window((50,70), 20, "+")
    (30, 49)
    >>> add_window((50,70), 20, "-")
    (71, 90)
    >>> add_window((50,70), 20, "+", fiveprime=False)
    (71, 90)
    >>> add_window((50,70), 20, "-", fiveprime=False)
    (30, 49)
    """
        
        
    if not strand == "-":
        if fiveprime:
            return (window[0]-size, window[0]-1)
        else:
            return (window[1]+1, window[1] + size)
    else:
        if fiveprime:
            return (window[1]+1, window[1] + size)
        else:
            return (window[0]-size, window[0]-1)

def get_random_window(window_length, genome_length):
    """
    Find a random window of length in the genome. Assumes a circular
    genome.

    Inputs:
        window_length - int length in bp of window size desired
        genome_length - int length of genome to sample from
    Returns:
        tuple(start, end) of random window in the genome

    >>> x = get_random_window(30, 400)
    >>> x[1] - x[0] +1 == 30
    True
    """
    bin_start = random.randint(0, genome_length)
    # still in 1 inclusive windows so need to subtract 1
    window = (bin_start, bin_start+window_length-1)
    return window

def simple_window_signal_circular(signal, window, genome_size):
    """ Given a signal numpy array from 0-genome size and a window in
    1-based inclusive coordinates, return a slice of the signal array within
    the window. 
    Assumes window is in [1-end] coordinates and signal is in
    [0-end) coordinates and the genome is circular

    Inputs:
        signal - 1d np.array of values for every bp in genome_size
        window - tuple(start, end) window to get slice of array from
        genome_size - size of the genome

    Output:
        slice of signal that window overlaps
    >>> signal = np.array([1, 2, 3, 4, 5, 6])
    >>> simple_window_signal_circular(signal, (1, 2), 6)
    array([1, 2])
    >>> simple_window_signal_circular(signal, (-1, 2), 6)
    array([5, 6, 1, 2])
    >>> simple_window_signal_circular(signal, (6, 8), 6)
    array([6, 1, 2])
    >>> simple_window_signal_circular(signal, (1, 1), 6)
    array([1])
    """
    # convert gff coordinate to 0 based coordinate
    window = list(window)
    window[0] = window[0] -1
    if window[0] < 0:
        window[0] = genome_size + window[0]
        window[1] = genome_size + window[1]
    if window[1] >= genome_size:
        window_sig = np.append(signal[window[0]:genome_size], signal[0:window[1]-genome_size])
    else:
        window_sig = signal[window[0]:window[1]]

    return window_sig

def complex_window_signal_circular(data, window, genome_size, return_locs=False):
    """ Given a signal numpy array from 0-genome size and a window in
    1-based inclusive coordinates, return a slice of the signal array within
    the window. 
    Assumes window is in [1-end] coordinates and signal is in
    [0-end) coordinates and the genome is circular

    Inputs:
        signal - 2d np.array of values for every bp in genome_size
        window - tuple(start, end) window to get slice of array from
        genome_size - size of the genome

    Output:
        slice of signal that window overlaps


    >>> locs = np.array([0, 1, 5, 35, 55, 99])
    >>> signal = np.array([1, 2, 3, 4, 5, 6])
    >>> data = np.column_stack([locs, signal])
    >>> complex_window_signal_circular(data, (1, 2), 100)
    array([1, 2])
    >>> complex_window_signal_circular(data, (37, 55), 100)
    array([], dtype=int64)
    >>> complex_window_signal_circular(data, (37, 56), 100)
    array([5])
    >>> complex_window_signal_circular(data, (36, 55), 100)
    array([4])
    >>> complex_window_signal_circular(data, (36, 56), 100)
    array([4, 5])
    >>> complex_window_signal_circular(data, (98, 102), 100)
    array([6, 1, 2])
    >>> complex_window_signal_circular(data, (-2, 2), 100)
    array([6, 1, 2])
    >>> complex_window_signal_circular(data, (1, 1), 100)
    array([1])
    >>> complex_window_signal_circular(data, (6, 6), 100)
    array([3])
    """
    signal = data[:,1]
    locs = data[:,0]
    # first location in window
    # change window to 0 based index
    startbp = window[0] - 1
    # last location in window
    # change window to 0 based index
    endbp = window[1]
    
    # search where in the locs array?
    searchleft_start = 0
    searchleft_end = len(locs)
    searchright_start = searchleft_start
    searchright_end = searchleft_end

    # make sure we have a sorted window:
    if startbp > endbp:
        startbp, endbp = endbp, startbp

    # deal with being outside the genome in a circular manner
    if startbp < 1:
        startbp += genome_size
    if endbp > locs[-1]:
        endbp -= genome_size
    # find what locations to pull
    startindex = bisect.bisect_left(locs, startbp, searchleft_start, searchleft_end)
    endindex = bisect.bisect_left(locs, endbp, searchright_start, searchleft_end)

    if (startindex > endindex):
        locs_slice = np.concatenate((locs[startindex:], locs[:endindex]))
        data_slice = np.concatenate((signal[startindex:], signal[:endindex]))
    else: 
        locs_slice = locs[startindex:endindex]
        data_slice = signal[startindex:endindex]

    if return_locs:
        return(np.column_stack((locs_slice, data_slice)))
    else:
        return(data_slice)


def sliding_bounds(size, length, slide_by=1):
    """
    Generator for a sliding window across a length.

    Input:
        size - int bp for the size of the window
        length - int length of the genome to slide over
        slide_by - int num of bp to slide over
    Output:
        windows in 1 based coordinates
    >>> gen = sliding_bounds(2, 4, 1)
    >>> gen.next()
    (1, 2)
    >>> gen.next()
    (2, 3)
    >>> gen.next()
    (3, 4)
    """
    # convert to 1 based coordinates
    for val in range(1, length+1, slide_by):
       yield (val, val+size-1)
