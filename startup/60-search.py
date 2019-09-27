# this file 60-search.py
# contains functions for searching in databases (of various types, including python dictionary)
# written by David Rais, 2017, david.rais@atlas.cz

def get_recursively(search_dict, field):
    """
    Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    """
    fields_found = []

    for key, value in search_dict.items():

        if key == field:
            fields_found.append(value)

        elif isinstance(value, dict):
            results = get_recursively(value, field)
            for result in results:
                fields_found.append(result)

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    more_results = get_recursively(item, field)
                    for another_result in more_results:
                        fields_found.append(another_result)

    return fields_found
    
def getScanNum(item):
    """ finds integer strings in text, and return the first one as integer """
    re_integer_machine = re.compile(r'\d+') # regular expression to find integer strings 
    integers = re_integer_machine.findall(item)
    if len(integers) != 0:
        num = int(re_integer_machine.findall(item)[0])
    else:
        num = 0
    return num

def get_float(s):
    floats = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", s)
    floats = [float(x) for x in floats]
    return floats

def find_analyses (whichDatasetName, GTArootDir, fname_pattern='overview.xml'):
    """Recursively traverses directory subtree, and searches the overview.xml files for dataset filename information, and extracts fit information from the summary text file.
    inputs: 
        whichDatasetName = part of original dataset file name
        GTArootDir = directory, where the glotaran results are located
    outputs dictionary with fileds:
        'datasetFname', (filename of the original dataset, without the extension)
        'dirName', (directory path to the results file, relative to GTArootDir)
        'fitName', (fit results filename, excluding the .timpres extension)
        'iterNum', (number of performed iterations)
        'residual_stderr', ...
        'kVect','kVect_stderr' (vectors of the best fit parameters and estimated standard errors...)
        'IRF', 'IRF_stderr' -ditto-
        'Parmu', 'Parmu_stderr' -ditto-
    """
    rootDir = GTArootDir + 'results'
    resultsDirs = []
    for dirName, subdirList, fileList in os.walk(rootDir):
        #print('Found directory: %s' % dirName)
        for fname in fileList:
            if fname_pattern in fname:
                datasetFname,resultFname = getDatasetResultFnamePart(dirName+'/'+fname)
                if whichDatasetName in datasetFname:
                    summary = {'datasetFname' : datasetFname,
                              'dirName': dirName[len(GTArootDir):],
                              'fitName':resultFname}
                    summary.update(getResultsFromSummary(dirName))
                    resultsDirs.append(summary)
    return(resultsDirs)    
