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
"""
fasta

Module to read, write, and manipulate fasta files.

"""

def complement(sequence):
    """Complement a nucleotide sequence
    >>> complement("AGTC")
    'TCAG'
    >>> complement("AGNT")
    'TCNA'
    >>> complement("AG-T")
    'TC-A'
    """
    # create a dictionary to act as a mapper
    comp_dict = {'A': 'T', 'G':'C', 'C':'G', 'T': 'A', 'N':'N', '-':'-'}
    # turn the sequence into a list
    sequence = list(sequence)
    # remap it to the compelmentary sequence using the mapping dict
    sequence = [comp_dict[base] for base in sequence]
    # join the new complemented sequence list into a string
    sequence = ''.join(sequence)
    return sequence

class FastaEntry(object):
    """ 
    Stores all the information for a single fasta entry. 

    An example of a fasta entry is below:

    >somechromosomename
    AGAGATACACACATATA...ATACAT #typically 50 bases per line

    Args:
        header (str): The complete string for the header ">somechromosomename" 
                      in the example above. Defaults to ">"
        seq (str): The complete string for the entire sequence of the entry

    Attributes:
        header (str): The complete string for the header ">somechromosomename" 
                      in the example above.
        seq (str): The complete string for the entire sequence of the entry

    """
    def __init__(self, header = ">", seq = ""):
        self.header = header
        self.seq = seq
        self.length = None

    def __str__(self):
        return "<FastaEntry>" + self.chrm_name() + ":" + str(len(self))

    def write(self,fhandle):
        fhandle.write(self.header+"\n")
        for i in range(0,len(self), 70):
            try:
                fhandle.write(self.seq[i:i+70]+"\n")
            except IndexError:
                fhandle.write(self.seq[i:-1] + "\n")


    def set_header(self, header):
        self.header = header

    def set_seq(self, seq):
        self.seq = seq

    def __len__(self):
        if self.length:
            return self.length
        else:
            return(len(self.seq))

    def pull_seq(self, start, end, circ=False, rc=False):
        """ 
        Obtain a subsequence from the fasta entry sequence

        Args:
            start (int)    : A start value for the beginning of the slice. Start
                             coordinates should always be within the fasta entry
                             sequence length.
            end   (int)    : An end value for the end of the slice. If circ is
                             True then end coordinates can go beyond the fasta
                             entry length.
            circ  (boolean): Flag to allow the end value to be specified beyond
                             the length of the sequence. Allows one to pull
                             sequences in a circular manner.
        Returns:
            A subsequence of the fasta entry as specified by the start and end

        Raises:
            ValueError: If start < 0 or >= fasta entry sequence length.
            ValueError: If circ is False and end > fasta entry sequence length.
        """

        seq_len = len(self)
        if start < 0 or start >= seq_len:
            if circ:
                if start < 0:
                    start = seq_len + start
                    end = seq_len + end
                elif start >= seq_len:
                    start = start - seq_len
                    end = end - seq_len
            else:
                raise ValueError("Start %s is outside the length of the sequence %s"%(start,seq_len))
        if end > seq_len:
            if circ:
                seq = self.seq[start:seq_len] + self.seq[0:(end-seq_len)]
            else: 
                raise ValueError("End %s is outside length of sequence %s"%(end,seq_len))
        else:
            seq = self.seq[start:end]
        if rc:
            return complement(seq)[::-1]
        else:
            return seq.upper()

    def chrm_name(self):
        """
        Pulls the chromosome name from the header attribute.

        Assumes the header is of the type ">chromosomename" and nothing else
        is in the header.

        Returns:
            chromosome name
        """
        return self.header[1:]
                

class FastaFile(object):
    """ 
    Stores all the information for a single fasta file.

    An example of a fasta file is below

    >somechromosomename
    AGAGATACACACATATA...ATACAT 
    GGGAGAGAGATCTATAC...AGATAG
    >anotherchromosomename
    AGAGATACACACATATA...ATACAT #typically 50 bases per line

    Attributes:
        data (dict): where the keys are the chromosome names and the entries are
                     FastaEntry objects for each key
    """

    def __init__(self):
        self.data = {}
        self.names = []

    def __iter__(self):
        for name in self.names:
            yield self.pull_entry(name)


    def read_whole_file(self, fhandle):
        """ 
        Read an entire fasta file into memory and store it in the data attribute
        of FastaFile

        Args:
            fhandle (File)    : A python file handle set with mode set to read
        Returns:
            None

        Raises:
            ValueError: If fasta file does not start with a header ">"
        """

        line = fhandle.readline().strip()
        if line[0] != ">":
            raise ValueError("File is missing initial header!")
        else:
            curr_entry = FastaEntry(header = line.rstrip().split()[0])
        line = fhandle.readline().strip()
        curr_seq = []
        while line != '':
            if line[0] == ">":
                curr_entry.set_seq(''.join(curr_seq))
                self.data[curr_entry.chrm_name()] = curr_entry
                self.names.append(curr_entry.chrm_name())
                curr_seq = []
                curr_entry = FastaEntry(line)
            else:
                curr_seq.append(line)

            line = fhandle.readline().strip()

        curr_entry.set_seq(''.join(curr_seq))
        self.data[curr_entry.chrm_name()] = curr_entry

    def pull_entry(self, chrm):
        """
        Pull a FastaEntry out of the FastaFile

        Args:
            chrm (str): Name of the chromosome that needs pulled
        Returns:
            FastaEntry object
        """
        try:
            return self.data[chrm]
        except KeyError:
            raise KeyError("Entry for chromosome %s does not exist in fasta file"%chrm)

    def add_entry(self, entry):
        """
        Add a FastaEntry to the object

        Args:
            entry (FastaEntry): FastaEntry to add
        Returns:
            None
        """
        self.data[entry.chrm_name()]= entry
        self.names.append(entry.chrm_name())

    def chrm_names(self):
        return self.data.keys()

    def write(self, fhandle):
        """ 
        Write the contents of self.data into a fasta format

        Args:
            fhandle (File)    : A python file handle set with mode set to write
        Returns:
            None

        """
        for chrm in self.chrm_names():
            entry = self.pull_entry(chrm)
            entry.write(fhandle)
