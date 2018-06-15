# python 2.7
import imp

ipod_utils = imp.load_source('ipod_utils', 'ipod_utils.py')

# fitting HupA
print('fitting HupA')
print(ipod_utils.spline_correct_genome.__doc__)
# do a spline-based correction of any periodicity in the input file
# turn off plotting as unneccessary and avoiding Runtime Error
ipod_utils.spline_correct_genome('65_coverage/HupA/_count.gr', '67_coverageNormalized/HupA/HupA.coverage.spline.gr')

# fitting HupB
print('fitting HupB')
ipod_utils.spline_correct_genome('65_coverage/HupB/_count.gr', '67_coverageNormalized/HupB/HupB.coverage.spline.gr')
