#!usr/bin/python

# Tools for managing gff files
# These are really just files that contain a set of 9 tab delimited fields:
#  Genome_name Data_origin Site_type start end . +/-(direction) . comments
#
# Copyright (c) 2018 Peter Freddolino Michael Wolfe University of Michigan. All
# rights reserved.
#
#
#Developed by: Peter Freddolino, Michael Wolfe
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
#the distribution.  Neither the names of Peter Freddolino, Michael Wolfe,
#University of Michigan, nor the names of its contributors may be used to
#endorse or promote products derived from this Software without specific prior
#written permission.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS
#OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.

class GffEntry:
  """
  A simple container class for the equivalent of a gff file line
  """

  FORMAT_STRING = "%s\t%s\t%s\t%i\t%i\t.\t%s\t.\t%s"

  def __init__(self, line=None):

    if (line):
      self.parse_gff_line(line)
    else:
      self.genome_name = ""
      self.data_origin = ""
      self.site_type = ""
      self.start = 0
      self.end = 0
      self.direction = "+"
      self.comments = ""


  def parse_gff_line(self, line):
    """
    Set this entry's values to those of a line from a gff file
    """

    datarray = line.split("\t")
    self.genome_name = datarray[0]
    self.data_origin = datarray[1]
    self.site_type = datarray[2]
    self.start = int(datarray[3])
    self.end = int(datarray[4])
    self.direction = datarray[6]
    self.comments = " ".join(datarray[8:])
    self.comments = self.comments.replace("\t", " ")
    self.comments = self.comments.replace("\n", "")

  def __repr__(self):
    """
    Return a formatted gff line, which can be used to reconstitute the object or be written directly to a gff file
    """

    return GffEntry.FORMAT_STRING % (self.genome_name, self.data_origin, self.site_type, self.start, self.end, self.direction, self.comments)



class GffData:
  """
  Class for storing and manipulating gff data
  """

  def __init__(self):
    self.data = []
    self.index = 0

  def __iter__(self):
    return self

  def next(self):
    if self.index == len(self.data):
      self.index = 0
      raise StopIteration
    self.index = self.index + 1
    return self.data[self.index-1]



  def clear_db(self):
    self.data = []

  def parse_gff_file(self,filename, clear=True):
    """
    Parse a gff file and store the lines in self.data

    The current contents of this object are overwritten iff clear 
    """

    if (clear):
      self.clear_db()

    instr = open(filename, "r")
    for line in instr:
      newline = GffEntry(line)
      self.data.append(newline)

  def cleanup(self):
    """
    Remove all duplicate entries and sort based on starting position
    """

    self.data = list(set(self.data))
    self.data.sort(cmp=lambda a,b: cmp(a.start, b.start))

  def write_gff_file(self, filename):
    """
    Write the current contents of my data to a file
    """

    ostr = open(filename, "w")

    for line in self:
      ostr.write("%s\n" % line)

    ostr.close()

  def addline(self, genome_name, data_origin, site_type, start, end, direction, comments):
    """
    Add a line with the given data
    """

    newobj = GffEntry()
    newobj.genome_name = genome_name
    newobj.data_origin = data_origin
    newobj.site_type = site_type
    newobj.start = start
    newobj.end = end
    newobj.direction = direction
    newobj.comments = comments

    self.data.append(newobj)

  def find_entry(self, findfunc, findall=False):
    """
    Return the line or lines for which findfunc is true when given a gff item

    If findall is false, only the first such entry is returned

    Findfunc should take a gffline object as its only argument
    """

    matches = filter(findfunc, self.data)

    if len(matches) == 0:
      return []

    if (findall):
      return matches
    else:
      return matches[0]

  def add_entry(self, new_entry):
    # Add an externally constructed GffEntry object

    self.data.append(new_entry)


