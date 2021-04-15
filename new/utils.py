import os
import shelve

import universe

def update(tool, message):
    print("{:22s} : {:s}".format(tool, message))

def generate_index():
    os.system("gmx make_ndx -f {0} >> builder.log 2>&1 << EOF\nq\nEOF".format(universe.get('d_nameList')[-1]))

def add_to_nameList(name):
    if universe.has_key('d_nameList'):
        temp = universe.get('d_nameList')
        temp.append(name)
        universe.add('d_nameList', temp)
    else:
        universe.add('d_nameList', [name])
