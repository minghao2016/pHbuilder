import os
import shelve

import universe

def update(tool, message):
    print("{:22s} : {:s}".format(tool, message))

def generate_index():
    os.system("gmx make_ndx -f {0} >> builder.log 2>&1 << EOF\nq\nEOF".format(universe.get('d_nameList')[-1]))

def add_to_nameList(name):
    with shelve.open('universe') as shelf:
        # If d_nameList does not exist, add it as an empty list.
        try:
            shelf['d_nameList']
        except KeyError:
            shelf['d_nameList'] = []
    
        # Append name to d_nameList.
        temp = shelf['d_nameList']
        temp.append(name)
        shelf['d_nameList'] = temp
