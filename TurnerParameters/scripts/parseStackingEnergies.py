'''
Filename: parseWatsonCrickTurnerParameters.py
Author: Michael Hathaway
Date: 11/11/2019

Description: file to parse .txt file with data on watson crick stacking energies.

file parses the data into a dictionary of dictionaries where the key to the first dictionary
is a two string tuple containing the preceding base pairs, the key in the next dictionary is a 
two string tuple containg the following base pairs, and the stored value is the stacking energy associated
with that arrangement of base pairs.

This dictionary is written to another .py file that can be imported and used in other python scripts

Usage:
python3 parseWatsonCrickTurnerParameters.py <data text file> <output file name>
'''

import argparse

'''
Function: parseTurnerParametersStacking(filename)
Description:
parameters:
Return Type:
'''
def parseTurnerParametersStacking(filename):
	with open(filename, 'r') as f:
		#skip header and start at line 23
		for i in range(23): 
			next(f)

		text = f.read()
		rows = []
		for line in text.split('\n'): #break text into lines
			rows.append(line.rstrip('\r\n').lstrip()) #add lines after stripping extra characters and spaces

		block1 = [line.split() for  line in (rows[0:2] + rows[3:7])]
		block2 = [line.split() for  line in (rows[15:17] + rows[18:22])]
		block3 = [line.split() for  line in (rows[30:32] + rows[33:37])]
		block4 = [line.split() for  line in (rows[44:46] + rows[47:51])]
	
		blocks = [block1, block2, block3, block4]
		
		return _parseBlocks(blocks)
		

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
		

def writeToPythonFile(stackDict, dictName):
	with open(f'{dictName}.py', 'w') as f:
		f.write(f'{dictName} = {stackDict}')


def parseArgs():
	parser = argparse.ArgumentParser(description="Turner Parameters Parser for Watson Crick Stacking values.")	
	parser.add_argument('Input_File', help="Specify file to be parsed.", type=str)
	parser.add_argument('Dictionary_Name', help="specify name of the dictionary and file to write dictionary to.", type=str)

	args = parser.parse_args()
	return (args.Input_File, args.Dictionary_Name)



if __name__ == '__main__':
	args = parseArgs()
	stackDict = parseTurnerParametersStacking(args[0])
	writeToPythonFile(stackDict, args[1])
	


