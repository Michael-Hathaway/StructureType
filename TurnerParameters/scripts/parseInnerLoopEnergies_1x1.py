'''
Filename: parseInnerLoopEnergies_1x1.py
Author: Michael Hathaway

Description: script to parse the energy parameters for 1x1 internal loops from the turner parameter text file

Usage: python parseInnerLoopEnergies_1x1.py <input filename> <output file name without .py ending>
'''

import sys
import argparse

'''
Function: parseInnerLoopEnergies(filename)
Description: Function parse the energy parameters for 1x1 internal loops from the turner parameter text file
Parameters: (filename) -- string -- name of the file to be parsed
Return Type: dictionary
'''
def parseInnerLoopEnergies(filename):
    try:
        f = open(filename, 'r')
    except:
        print(f'An error occurred when trying to acess the file: {filename}')
        sys.exit()


    #move past first 23 lines in the file
    for i in range(21):
        next(f)

    contents = f.read() #read blocks
    contents = contents.split('\n') #split blocks by new line

    blocks = []
    for i in range(len(contents)): #break file into list of seperate parameter blocks
        if contents[i].count('Y') == 6 and contents[i+1].count('-') > 20:
            blocks.append(contents[i:i+14])

    bases = ['A', 'C', 'G', 'U'] #base label

    energyDict = {} #dictionary to store energy values
    # values will be stored in dictioanry as (5' closing pair, 3' closing pair, x, y)

    for block in blocks: #iterate through each block
        labelBlock1 = list(filter(lambda arg: arg.isalpha(), block[6]))
        labelBlock2 = list(filter(lambda arg: arg.isalpha(), block[7]))

        labels = list(zip(labelBlock1, labelBlock2)) #labels for inner loop closing pairs
        values = block[10:14] #rows with parameter value

        for i in range(len(values)): #iterate through each line of values(4 lines)
            lineContents = values[i].split(' ') #split line of values into a list of values
            lineContents = list(filter(lambda arg: (arg == '') == False, lineContents)) #filter out empty strings
            counter = 0

            for j in range(len(lineContents)): #iterate through the list of values
                # bases[i] will give the value of x
                # j % 4 will gove value of y

                energyDict[labels[counter], labels[counter+1], bases[i], bases[j%4]] = float(lineContents[j])

    return energyDict

'''
Function: writeToFile(dictionary, outputName)
Description: Function to write the parsed file contents to python file
Parameters: (dictionary) - the dictionary to be written to the python file
            (outputName) - Name of the output file without the .py ending
Return Type: None
'''
def writeToFile(dictionary, outputName):
    with open(f'{outputName}.py', 'w') as f:
        f.write(f'{outputName} = {dictionary}')



'''
Function: parseArgs()
Description: Function to handle command line arguments
parameters: None
Return Type: tuple containing the input file name and the output file name
'''
def parseArgs():
    parser = argparse.ArgumentParser(description="Turner Parameters Parser for 1x1 Inner Loop Energies.")
    parser.add_argument('Input_File', help="Specify file to be parsed.", type=str)
    parser.add_argument('Dictionary_Name', help="specify name of the dictionary and file to write dictionary to.", type=str)

    args = parser.parse_args()
    return (args.Input_File, args.Dictionary_Name)



## Main function ##
if __name__ == '__main__':
    args = parseArgs()
    d = parseInnerLoopEnergies(args[0])
    writeToFile(d, args[1])
