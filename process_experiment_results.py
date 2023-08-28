# -*- coding: utf-8 -*-

def process_input(line, results):
    # split line on whitespaces
    line = line.split()
    
    # # ignore lines with less than 29 successful runs. These don't have avgEvals, so line should be shorter.
    # if(len(line) <= 4):
    #     return

    # F1_F4_dict = {} #{"GAMBIT_R": {}, "GAMBIT_K": {}}
    # F5_dict = {} #{"GAMBIT_R": {}, "GAMBIT_K": {}}
    # first element is path to file, split on '/'
    relevant_path = line[0][line[0].find('GAMBIT'):]
    path = relevant_path.split('/')

    # for F1-F4, path is 9 elements long. But to make sure, check the file name
    # if(path[3] in ["F1", "F2", "F3", "F4"]):
        # # process F1-F4. First make sure the appropriate (nested) dictionary already exists
        # if(path[3] not in F1_F4_dict[path[0]]):
        #     F1_F4_dict[path[0]][path[3]] = {}
        # if(path[4] not in F1_F4_dict[path[0]][path[3]]):
        #     F1_F4_dict[path[0]][path[3]][path[4]] = {}
        # if(path[5] not in F1_F4_dict[path[0]][path[3]][path[4]]):
        #     F1_F4_dict[path[0]][path[3]][path[4]][path[5]] = {}

        # now add the data to the dictionary
    currDataSettings = f"{path[0]}__{path[3]}__{path[4]}__{path[5]}"

    # lines that didn't complete at least 29 runs just add empty tuple (if not already in dict)
    if(len(line) <= 4):
        if(currDataSettings not in results):
            results[currDataSettings] = (None, None, None)
        return

    currPopsize = int(line[3])
    currAvgEvals = float(line[5])
    popsize_mult = float(path[1].split("_")[1])
    if((currDataSettings not in results) or (results[currDataSettings][0] is None) or (currPopsize < results[currDataSettings][0]) ):
        results[currDataSettings] = (currPopsize, currAvgEvals, popsize_mult)
        

    # elif(path[5] == "F5"):
    #     # process F5
    #     pass


def main():
    # read in text file line by line
    results = {}
    with open('experiment_results.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            process_input(line, results)
    results_list = []
    for key in results:
        results_list.append(f"{key}: ({results[key][0]}, {results[key][1]}, {results[key][2]})")
    print(*results_list, sep="\n")
    

if __name__ == "__main__":
    main()