import os

def backupFile(fname):
    if os.path.isfile(fname):
        num = 1
        while (True):
            if os.path.isfile("#{0}.{1}#".format(fname, num)):
                num += 1
                continue

            os.system("mv {0} '#{0}.{1}#'".format(fname, num))
            break
