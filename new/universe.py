import shelve

# Set/update variable to universe.
def add(varName, value):
    with shelve.open('universe') as shelf:
        shelf[varName] = value

# Retrieve variable from universe.
def get(varName):
    with shelve.open('universe') as shelf:
        try:
            return shelf[varName]
        except KeyError:
            data = eval(input("GET couldn't retrieve var \"{0}\" from {1}. Enter manually: ".format(varName, 'universe')))
            shelf[varName] = data
            print("SET {0} = {1}  {2}".format(varName, shelf[varName], type(shelf[varName])))
            return shelf[varName]

# Display all variables (name, data, type) stored in the universe.
def inspect():
    with shelve.open('universe') as shelf:
        # Determine longest valueName for formatting:
        longest = 0
        for item in shelf:
            if (len(item) > longest):
                longest = len(item)
        
        for item in shelf:
            # If item is a long list, only print first, last element (to save screen space)
            if (type(shelf[item]) == type([]) and len(shelf[item]) > 2):
                print("{0:{arg}s} = [{1}, ..., {2}] ({3}) {4}".format(item, shelf[item][0], shelf[item][-1], len(shelf[item]), type([]), arg=longest).replace('\n', ''))
            else:
                print("{0:{arg}s} = {1}  {2}".format(item, shelf[item], type(shelf[item]), arg=longest))
