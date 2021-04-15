import shelve

# Set/update variable to universe.
def add(varName, value):
    with shelve.open('universe') as shelf:
        shelf[varName] = value

# Check whether universe contains a certain varName
def has(varName):
    with shelve.open('universe') as shelf:
        return varName in shelf

# Retrieve variable from universe.
def get(varName):
    if has(varName):
        return shelve.open('universe')[varName]

    data = eval(input("couldn't retrieve var \"{0}\" from universe. Enter manually: ".format(varName)))
    print("add {0} = {1} {2}".format(varName, data, type(data)))
    add(varName, data)
    return data

# Display all variables (name, data, type) stored in the universe.
def inspect():
    with shelve.open('universe') as shelf:
        # Determine longest valueName for formatting:
        longest = 0
        for item in shelf:
            if (len(item) > longest):
                longest = len(item)
        
        for item in sorted(shelf):
            # If item is a long list, only print first, last element (to save screen space)
            if (type(shelf[item]) == type([]) and len(shelf[item]) > 2):
                print("{0:{arg}s} = [{1}, ..., {2}] ({3}) {4}".format(item, shelf[item][0], shelf[item][-1], len(shelf[item]), type([]), arg=longest).replace('\n', ''))
            else:
                print("{0:{arg}s} = {1} {2}".format(item, shelf[item], type(shelf[item]), arg=longest))
