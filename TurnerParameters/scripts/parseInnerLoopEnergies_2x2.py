'''
Filename: parseInnerLoopEnergies_2x2.py
Description: script to parse the energy parameters for 2x2 internal loops from the turner parameter text file
Author: Michael Hathaway
Date: 11/29/2019
'''

import sys
import argparse

'''
Function:
Description:
Parameters:
Return Type:
'''
def parseInnerLoopEnergies(filename):
	try:
		f = open(filename, 'r')
	except:
		print(f'An error occurred when trying to acess the file: {filename}')
		sys.exit()


	#move past first 23 lines in the file
	for i in range(32):
		next(f)

	contents = f.read() #read blocks
	contents = contents.split('\n') #split blocks by new line

	blocks = []
	for i in range(len(contents)): #break file into list of seperate parameter blocks
		if contents[i].count('Y') == 1 and contents[i+1].count('-') > 20:
			blocks.append(contents[i:i+25])
	
	bases = [('A', 'A'), ('A', 'C'), ('A', 'G'), ('A', 'U'), ('C', 'A'), ('C', 'C'), ('C', 'G'), ('C', 'U'), ('G', 'A'), 
			('G', 'C'), ('G', 'G'), ('G', 'U'), ('U', 'A'), ('U', 'C'), ('U', 'G'), ('U', 'U')] #base label

	energyDict = {} #dictionary to store energy values
	# values will be stored in dictioanry as (5' closing pair, 3' closing pair, x pair, y pair)

	for b in range(len(blocks)):
		labelBlock1 = list(filter(lambda arg: arg.isalpha(), blocks[b][6])) #filter out spaces 
		labelBlock2 = list(filter(lambda arg: arg.isalpha(), blocks[b][7])) #filter out spaces

		labels = list(zip(labelBlock1, labelBlock2)) #get labels for closing pairs
		values = blocks[b][9:25] #rows with parameter value
		
		for i in range(len(values)): #iterate through each line of values(16 lines)
			lineContents = values[i].split(' ') #split line of values into a list of values
			lineContents = list(filter(lambda arg: (arg == '') == False, lineContents)) #filter out empty strings
			counter = 0
			
			for j in range(len(lineContents)): #iterate through the list of values	
				energyDict[labels[counter], labels[counter+1], bases[i], bases[j]] = float(lineContents[j])
		

	return energyDict

'''
Function:
Description:
Parameters:
Return Type:
'''
def writeToFile(dictionary, outputName):
	with open(f'{outputName}.py', 'w') as f:
		f.write(f'{outputName} = {dictionary}')

'''
Function:
Description:
Parameters:
Return Type:
'''
def parseArgs():
	parser = argparse.ArgumentParser(description="Turner Parameters Parser for 2x2 Inner Loop Energies.")	
	parser.add_argument('Input_File', help="Specify file to be parsed.", type=str)
	parser.add_argument('Dictionary_Name', help="specify name of the dictionary and file to write dictionary to.", type=str)

	args = parser.parse_args()
	return (args.Input_File, args.Dictionary_Name)


if __name__ == '__main__':
	args = parseArgs()
	d = parseInnerLoopEnergies(args[0])
	writeToFile(d, args[1])


