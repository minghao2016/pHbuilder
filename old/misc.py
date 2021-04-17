import os
import fnmatch

def inferFullName():
    filtered = fnmatch.filter(os.listdir('.'), "*_NPT.pdb")
    return filtered[0]

def inferName():
    filtered = fnmatch.filter(os.listdir('.'), "*_NPT.pdb")
    return filtered[0][0:len(filtered[0])-7]

def backupFile(fname):
    if os.path.isfile(fname):
        num = 1
        while (True):
            if os.path.isfile("#{0}.{1}#".format(fname, num)):
                num += 1
                continue

            os.system("mv {0} '#{0}.{1}#'".format(fname, num))
            break
