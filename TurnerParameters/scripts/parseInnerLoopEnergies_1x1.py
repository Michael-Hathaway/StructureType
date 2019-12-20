'''
Filename: parseInnerLoopEnergies_1x1.py
Description: script to parse the energy parameters for 1x1 internal loops from the turner parameter text file
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
	for i in range(21):
		next(f)

	contents = f.read() #read blocks
	contents = contents.split('\n') #split blocks by new line

	blocks = []
	blocks.append(contents[0:14])
	blocks.append(contents[15:29])
	blocks.append(contents[30:44])
	blocks.append(contents[45:59])
	blocks.append(contents[60:74])
	blocks.append(contents[75:89])

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
	parser = argparse.ArgumentParser(description="Turner Parameters Parser for 1x1 Inner Loop Energies.")	
	parser.add_argument('Input_File', help="Specify file to be parsed.", type=str)
	parser.add_argument('Dictionary_Name', help="specify name of the dictionary and file to write dictionary to.", type=str)

	args = parser.parse_args()
	return (args.Input_File, args.Dictionary_Name)


if __name__ == '__main__':
	args = parseArgs()
	d = parseInnerLoopEnergies(args[0])
	writeToFile(d, args[1])


