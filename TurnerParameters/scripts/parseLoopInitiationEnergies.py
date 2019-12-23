'''
Filename: parseLoopInitiationEnergies.py
Author: Michael Hathaway

Description: python script to parse loop initition turner parameters.
initiation paramters for inner loops, bulges, and hairpins are parsed into
a .py file that cannot be imported and used in other pythin scripts

Usage:
python3 parseLoopInitiationEnergies.py <data text file> <output filename>
'''

import argparse
import sys

'''
Function: parseLoopInitiation(Filename)
Description: Function to parse loop initiation parameters 
Parameters: (filename) -- string -- name of the file to be parsed
Return Type: Tuple of dictionaries
'''
def parseLoopInitiation(filename):
	try:
		f = open(filename, 'r')
	except:
		print(f'An error occurred when trying to acess the file: {filename}')
		sys.exit()

	for i in range(4):
		next(f)

	text = f.read()
	rows = []
	for line in text.split('\n'): #break text into lines
		rows.append(line.rstrip('\r\n').lstrip()) #add lines after stripping extra characters and spaces

	data = [line.split() for line in rows[:-1]] #convert each line into a row of characters

	internalLoopInit = {}
	bulgeInit = {}
	hairpinInit = {}

	for line in data:
		loopLen = int(line[0])

		if line[1] == '.':
			internalLoopInit[loopLen] = None
		else:
			internalLoopInit[loopLen] = float(line[1])

		if line[2] == '.':
			bulgeInit[loopLen] = None
		else:
			bulgeInit[loopLen] = float(line[2])

		if line[3] == '.':
			hairpinInit[loopLen] = None
		else:
			hairpinInit[loopLen] = float(line[3])

	return (internalLoopInit, bulgeInit, hairpinInit)


'''
Function: writeToFile(dictionary, outputName)
Description: Function to write the parsed file contents to python file
Parameters: (dictTuple) - tuple of dictionaries to be written to the python file
			(dictName) - Name of the output file without the .py ending
Return Type: None
'''
def writeToPythonFile(dictTuple, dictName):
	with open(f'{dictName}.py', 'w') as f:
		f.write(f'InternalLoopInit = {dictTuple[0]}')
		f.write('\n')
		f.write(f'BulgeInit = {dictTuple[1]}')
		f.write('\n')
		f.write(f'HairpinInit = {dictTuple[2]}')



'''
Function: parseArgs()
Description: Function to handle command line arguments
parameters: None
Return Type: tuple containing the input file name and the output file name
'''
def parseArgs():
	parser = argparse.ArgumentParser(description="Turner Parameters Parser for loop initiation values.")
	parser.add_argument('Input_File', help="Specify file to be parsed.", type=str)
	parser.add_argument('Dictionary_Name', help="specify name of the dictionary and file to write dictionary to.", type=str)

	args = parser.parse_args()
	return (args.Input_File, args.Dictionary_Name)


## Main function ##
if __name__ == '__main__':
	args = parseArgs()
	dictTuple = parseLoopInitiation(args[0])
	writeToPythonFile(dictTuple, args[1])
