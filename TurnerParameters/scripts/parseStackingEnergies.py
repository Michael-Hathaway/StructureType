'''
Filename: parseWatsonCrickTurnerParameters.py
Author: Michael Hathaway
Date: 11/11/2019

Description: file to parse .txt file with data on watson crick stacking energies.
The script parses the data into a dictionary of dictionaries where the key to the first dictionary
is a two string tuple containing the preceding base pairs, the key in the next dictionary is a
two string tuple containg the following base pairs, and the stored value is the stacking energy associated
with that arrangement of base pairs. This dictionary is written to another .py file that can be imported and used in other python scripts

Usage:
python3 parseWatsonCrickTurnerParameters.py <data text file> <output file name without .py ending>
'''

import argparse
import sys

'''
Function: parseTurnerParametersStacking(filename)
Description: Function splits file contents into 4 blocks that contain the stacking energy data
parameters: (filename) -- string -- name of the file to be parsed
Return Type: Dictionary representing the file information
'''
def parseTurnerParametersStacking(filename):
	try: #try to open the provided file
		f = open(filename, 'r')
	except: #error handling
		print(f'An error occurred when trying to acess the file: {filename}')
		sys.exit()

	#skip header and start at line 18
	for i in range(18):
		next(f)

	#read the rest of the file contents
	text = f.read()
	rows = []
	for line in text.split('\n'): #break text into lines
		rows.append(line.rstrip('\r\n').lstrip()) #add lines after stripping extra characters and spaces

	blocks = []
	for i in range(len(rows)): #break file into list of seperate parameter blocks
		if rows[i].count('Y') == 4 and rows[i+1].count('-') > 20:
			blocks.append([line.split() for  line in (rows[i+5:i+7] + rows[i+8:i+12])]) #only get labels and values

	return _parseBlocks(blocks)


'''
Function: _parseBlocks(blocks)
Description: Function sorts blocks produced in the parseTurnerParametersStacking() function into a dictionary.
parameters: (blocks) -- list -- list of stacking energy data produced by parseTurnerParametersStacking() function
Return Type: Dictionary
'''
def _parseBlocks(blocks):
	stackDict = {}

	#iterate through blocks
	for block in blocks:
		#get labels for the preceding base pairs
		label1 = [pair[0] for pair in block[0]]
		label2 = [pair[0] for pair in block[1]]
		basePairLabels = list(zip(label1, label2))

		for label in basePairLabels:
			if label not in stackDict.keys():
				stackDict[label] = dict()

		nucleotideLabels = ['A', 'C', 'G', 'U'] #labels for rows and columns in file
		for i in range(2, len(block)): #iterate through 4 value lines in block
			for index, value in enumerate(block[i]):
				# i-2 in nucleiotideLabels will be column label
				# index%4 in nucleotide labels will be row label
				# index//4 will be preceding pair from basePairLabels
				precedingLabel = basePairLabels[index//4]
				followingLabel = (nucleotideLabels[i-2], nucleotideLabels[index%4])
				if value != '.':
					stackDict[precedingLabel][followingLabel] = float(value)

	return stackDict

'''
Function: writeToPythonFile(stackDict, dictName)
Description: Function creates a new python file and writes the stacking energy dictionary to it.
parameters: (stackDict) - dictionary contents produced by the parseTurnerParametersStacking() function
			(dictName) - name of the dictionary and the file to be produced
Return Type: None
'''
def writeToPythonFile(stackDict, dictName):
	with open(f'{dictName}.py', 'w') as f:
		f.write(f'{dictName} = {stackDict}')

'''
Function: parseArgs()
Description: Function to handle command line arguments
parameters: None
Return Type: tuple containing the input file name and the output file name
'''
def parseArgs():
	parser = argparse.ArgumentParser(description="Turner Parameters Parser for Watson Crick Stacking values.")
	parser.add_argument('Input_File', help="Specify file to be parsed.", type=str)
	parser.add_argument('Dictionary_Name', help="specify name of the dictionary and file to write dictionary to.", type=str)

	args = parser.parse_args()
	return (args.Input_File, args.Dictionary_Name)


## Main Function ##
if __name__ == '__main__':
	args = parseArgs()
	stackDict = parseTurnerParametersStacking(args[0])
	writeToPythonFile(stackDict, args[1])
