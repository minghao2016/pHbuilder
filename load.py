def Int(fileName, line, col):
    for x, y in enumerate(open(fileName)):
        if x == line - 1:
            return int(y.split()[col-1])

def Float(fileName, line, col):
    for x, y in enumerate(open(fileName)):
        if x == line - 1:
            return float(y.split()[col-1])

def Str(fileName, line, col):
    for x, y in enumerate(open(fileName)):
        if x == line - 1:
            return y.split()[col-1]

def IntList(fileName):
    return [int(val) for val in open(fileName).read().splitlines()]

def FloatList(fileName):
    return [float(val) for val in open(fileName).read().splitlines()]

def StrList(fileName):
    return open(fileName).read().splitlines()

def Col(fileName, col, start=0, stop=0):
    data = []
    
    try:
        for x, y in enumerate(open(fileName)):
            if start == 0 and stop == 0 and (y.split()[0][0] not in ['@','#','&']):
                data.append(float(y.split()[col-1]))

            elif (x >= start-1) and (x <= stop-1):
                data.append(float(y.split()[col-1]))
    
    except IndexError:
        pass

    return data
