'''
Filename: StructureType.py
Author: Michael Hathaway

Description: python module that defines the structureType Object.
The structureType Object provides a user friendly mechanism for working with
RNA structure type files in the python programming language.
'''

## Module Imports ##
import numpy as np
import re

## Free Energy Parameter Imports ##
from TurnerParameters.parameters.LoopInitiationEnergy import InternalLoopInit, BulgeInit, HairpinInit
from TurnerParameters.parameters.StackingEnergies import StackingEnergies
from TurnerParameters.parameters.InnerLoop_1x1_Energies import InnerLoop_1x1_Energies
from TurnerParameters.parameters.InnerLoop_1x2_Energies import InnerLoop_1x2_Energies
from TurnerParameters.parameters.InnerLoop_2x2_Energies import InnerLoop_2x2_Energies
from TurnerParameters.parameters.InnerLoopMismatches import InnerLoopMismatches_2x3

## Other Free Energy Parameter Constants ##
INTERMOLECULAR_INIT = 4.09 #intermolecular initiation value
STEM_SYMMETRY_PENALTY = 0.43
STEM_AU_END_PENALTY = 0.45
INNER_LOOP_ASYMMETRY_PENALTY = 0.6
SPECIAL_C_BULGE = -0.9
BULGE_AU_END_PENALTY = 0.45

'''
## About the structureTypeObject ##

'''
class StructureType:

	def __init__(self, filename=None):
        	#RNA Molecule basic info
		#all values are stored as strings
		self._name = None
		self._length = None
		self._pageNum = None

		#structural representations of RNA molecule
		#all values are stored as strings
		self._sequence = None
		self._DBNotation = None
		self._structureArray = None
		self._varna = None

		#secondary structure information
		'''
		secondary structure information for each RNA molecule is stored as a dictionary where
		the key is the label for the secondary structure and the value is an object that contains
		the data for the secondary structure with the given label.

		Note: multiloops, innerloops, and pseudoknots are stored as nested dictionaries where the first
		key accesses the given structure, the second key accesses the particular subunit of the
		structure, and the value being stored is an object that contains the data for the structure subunit

		For information on how each secondary structure class was implemented and how to access
		secondary structure data, see the StructureTypeComponents section below.
		'''
		self._stems = {}
		self._hairpins = {}
		self._bulges = {}
		self._innerLoops = {}
		self._multiLoops = {}
		self._externalLoops = {}
		self._pk = {}
		self._ncbp = {}
		self._ends = {}

		'''
		component Array
		The component array is a numpy array of the same length as the molecule where each index
		contains the label for the secondary structure tha that index is a part of.
		the component array is initialized as None. When the length of the molecule is
		parsed from the .st file, a numpy array of that length is generated
		'''
		self._componentArray = None

		#load data from file if file is specified by user
		if filename != None:
			self._loadFile(filename)


	#if you use the structureTypeObject in a print() statement it will return the RNA molecule name.
	def __str__(self):
		return f'RNA: {self._name}'


###########################
###### Load File ##########
###########################

	'''
	Function Name: loadFile(filename)
	Description: user accessible function that can be used to load data from a structure type file into
	the StructureType object if no file is provided at object instantiation.
	Parameters: (filename) - str - name of the structure type file to be loaded into the object
	Return Type: None
	'''
	def loadFile(self, filename):
		self._loadFile(filename)


	'''
	Function Name: _loadFile(filename)
	Description: Internal method to parse the data in an RNA structure tyoe file into a StructureType object
	Parameters: (filename) - str, name of the file to be parsed
	Return Type: StructureType

	#Note: At this point, the data on segments and pseudoknots is not parsed from the .st file
	'''
	def _loadFile(self, filename):

		# check that file is valid structure type
		if filename[-3::] != '.st':
			print('Must provide a valid structure type file.')
			return None

		#try to open the provided file + error handling
		try:
			f = open(filename, 'r')
		except OSError:
			print('An error ocurred when trying to access the file. Check to make sure that the file exists and that the correct filepath was provided.')
			return None
		except:
			print('Something unexpected ocurred when accessing the file')
			return None

		lineCounter = 1 #using counter to identify lines that do not have label but are always in the same position
		for line in f:
			#get name of RNA molecule
			if line[0:6] == '#Name:':
				self._name = line[7:-1:]

			#get length of the RNA sequence
			elif line[0:8] == '#Length:':
				self._length = int(line[10:-1:])
				self._componentArray = np.empty(self._length, dtype=object)

			#get page number for molecule
			elif line[0:12] == '#PageNumber:':
				self._pageNum = int(line[13:-1:])

			#get actual RNA sequence
			elif lineCounter == 5:
				self._sequence = line[:-1] #drop the newline character

			#get Dot Bracket Notation for the molecule
			elif lineCounter == 6:
				self._DBNotation = line[:-1] #drop the newline character

			#get Annotated symbol form of the molecule
			elif lineCounter == 7:
				self._structureArray = line[:-1] #drop the newline character

			#not really sure what this notation is ???
			elif lineCounter == 8:
				self._varna = line[:-1] #drop the newline characters

				features = f.read() #read the rest of the file into features variable
				break

			lineCounter += 1


		features = features.split('\n') #split rest of file contents into a list of strings
		features = features[:-1] #remove the newline string at the end of the list

		for i in range(len(features)): #iterate through the individual string

			##stems##
			if features[i][0] == 'S' and features[i][1].isdigit():
				self._parseStemData(features[i].split(' '))

			##Hairpins##
			elif features[i][0] == 'H':
				self._parseHairpinData(features[i].split(' '))

			##Bulges##
			elif features[i][0] == 'B':
				self._parseBulgeData(features[i].split(' '))

			##Inner Loops##
			elif features[i][0] == 'I' and re.search('I\d{1,3}.1', features[i]):
				self._parseInnerLoopData(features[i].split(' '), features[i+1].split(' ')) #pass both inner loop components

			##Multiloops##
			elif features[i][0] == 'M':
				self._parseMultiloopData(features[i].split(' '))

			##external loops##
			elif features[i][0] == 'X':
				self._parseExternalLoopData(features[i].split(' '))

			##NCBP##
			elif features[i][0:4] == 'NCBP':
				self._parseNCBPData(features[i].split(' '))

			##Ends##
			elif features[i][0] == 'E':
				self._parseEndData(features[i].split(' '))


		f.close() #close the file



	'''
	Function Name: _parseStemData()
	Description: Internal method used by _loadFile() to parse all the stem information in the
	structure type file into the StructureType object
	Parameters: stemData - list - list of data from a line in the structure type file describing a stem
	Return Type: None
	'''
	def _parseStemData(self, stemData):
		#get stem Label
		stemLabel = stemData[0]

		#get start of first segment
		part1_start = ''
		for char in stemData[1]:
			if char.isnumeric():
				part1_start += char
			else:
				break
		part1_start = int(part1_start)

		#get stop of first segment
		part1_stop = ''
		for char in reversed(stemData[1]):
			if char.isnumeric():
				part1_stop += char
			else:
				break
		part1_stop = int(part1_stop[::-1])

		#get sequence for first segment
		part1_seq = ''
		for char in stemData[2]:
			if char.isalpha():
				part1_seq += char

		#get start of second segment
		part2_start = ''
		for char in stemData[3]:
			if char.isnumeric():
				part2_start += char
			else:
				break
		part2_start = int(part2_start)

		#get stop of second segment
		part2_stop = ''
		for char in reversed(stemData[3]):
			if char.isnumeric():
				part2_stop += char
			else:
				break
		part2_stop = int(part2_stop[::-1])

		#get sequence for second segment
		part2_seq = ''
		for char in stemData[4]:
			if char.isalpha():
				part2_seq += char

		#add data to the stems dictionary
		newStem = Stem(stemLabel, part1_seq, part2_seq, (part1_start, part1_stop), (part2_start, part2_stop))
		self._addStemToComponentArray(newStem)
		self.addStem(stemLabel, newStem)



	'''
	Function Name: _parseHairpinData()
	Description: Internal method used by _loadFile() to parse all the hairpin information in the
	structure type file into the StructureType object
	Parameters: hairpinData - list - list of data from a line in the structure type file describing a hairpin
	Return Type: None
	'''
	def _parseHairpinData(self, hairpinData):
		#get hairpin label
		hairpinLabel = hairpinData[0]

		#get index of start of hairpin
		hairpin_start = ''
		for char in hairpinData[1]:
			if char.isnumeric():
				hairpin_start += char
			else:
				break
		hairpin_start = int(hairpin_start)

		#get index of end of hairpin
		hairpin_stop = ''
		for char in reversed(hairpinData[1]):
			if char.isnumeric():
				hairpin_stop += char
			else:
				break
		hairpin_stop = int(hairpin_stop[::-1])

		#get sequence of hairpin
		hairpin_seq = ''
		for char in hairpinData[2]:
			if char.isalpha():
				hairpin_seq += char

		#get index of 5' closing base
		close_5_prime_index = ''
		for char in hairpinData[3]:
			if char.isnumeric():
				close_5_prime_index += char
			elif char == ',':
				break
		close_5_prime_index = int(close_5_prime_index)

		#get index of 3' closing base
		close_3_prime_index = ''
		for char in reversed(hairpinData[3]):
			if char.isnumeric():
				close_3_prime_index += char
			elif char == ',':
				break
		close_3_prime_index = int(close_3_prime_index[::-1])

		close_5_prime_base = hairpinData[4][0] #get 5' closing base
		close_3_prime_base = hairpinData[4][2] #get 3' closing base

		#get pk info
		if hairpinData[5] != '':
			pk = hairpinData[5][3]
		else:
			pk = None


		newHairpin = Hairpin(hairpinLabel, hairpin_seq, (hairpin_start, hairpin_stop),
			(close_5_prime_base, close_3_prime_base), (close_5_prime_index, close_3_prime_index), pk)
		self._addHairpinToComponentArray(newHairpin)
		self.addHairpin(hairpinLabel, newHairpin)



	'''
	Function Name: _parseBulgeData()
	Description: Internal method used by _loadFile() to parse all the bulge information in the
	structure type file into the StructureType object
	Parameters: bulgeData - list - list of data from a line in the structure type file describing a bulge
	Return Type: None
	'''
	def _parseBulgeData(self, bulgeData):
		#get bulge label
		bulgeLabel = bulgeData[0]

		#get start index of bulge
		bulge_start = ''
		for char in bulgeData[1]:
			if char.isnumeric():
				bulge_start += char
			else:
				break
		bulge_start = int(bulge_start)

		#get stop index of bulge
		bulge_stop = ''
		for char in reversed(bulgeData[1]):
			if char.isnumeric():
				bulge_stop += char
			else:
				break
		bulge_stop = int(bulge_stop[::-1])

		#get bulge sequence
		bulge_seq = ''
		for char in bulgeData[2]:
			if char.isalpha():
				bulge_seq += char

		#get index of 5' base in preceding pair
		precedingPair5pIndex = ''
		for char in bulgeData[3]:
			if char.isnumeric():
				precedingPair5pIndex += char
			elif char == ':':
				break
		precedingPair5pIndex = int(precedingPair5pIndex)

		#get index of 3' base in preceding pair
		precedingPair3pIndex = ''
		for char in reversed(bulgeData[3]):
			if char.isnumeric():
				precedingPair3pIndex += char
			elif char == ':':
				break
		precedingPair3pIndex = int(precedingPair3pIndex[::-1])

		#get 5' base of preceding pair
		precedingPair5pBase = bulgeData[4][0]

		#get 3' base in preceding pair
		precedingPair3pBase = bulgeData[4][2]

		#get index of 5' base in trailing pair
		trailingPair5pIndex = ''
		for char in bulgeData[5]:
			if char.isnumeric():
				trailingPair5pIndex += char
			elif char == ':':
				break
		trailingPair5pIndex = int(trailingPair5pIndex)

		#get index of 3' base in trailing pair
		trailingPair3pIndex = ''
		for char in reversed(bulgeData[5]):
			if char.isnumeric():
				trailingPair3pIndex += char
			elif char == ':':
				break
		trailingPair3pIndex = int(trailingPair3pIndex[::-1])

		#get 5' base of trailing pair
		trailingPair5pBase = bulgeData[6][0]

		#get 3' base in trailing pair
		trailingPair3pBase = bulgeData[6][2]

		#get pk info
		if bulgeData[7] != '':
			pk = bulgeData[7][3]
		else:
			pk = None


		newBulge = Bulge(bulgeLabel, bulge_seq, (bulge_start, bulge_stop), (precedingPair5pBase, precedingPair3pBase),
						(precedingPair5pIndex, precedingPair3pIndex), (trailingPair5pBase, trailingPair3pBase),
						(trailingPair5pIndex, trailingPair3pIndex), pk)
		self._addBulgeToComponentArray(newBulge)
		self.addBulge(bulgeLabel, newBulge)




	'''
	Function Name: _parseInnerLoopData()
	Description: Internal method used by _loadFile() to parse all the inner loop information in the
	structure type file into the StructureType object
	Parameters: innerLoopData - list - list of data from a line in the structure type file describing an inner loop
	Return Type: None
	'''
	def _parseInnerLoopData(self, loop1, loop2):
		#get Inner Loop parent Label
		parentLabel = ''
		for char in loop1[0]:
			if char == '.':
				break
			else:
				parentLabel += char

		#get inner loop subunit label
		loop1SubunitLabel = loop1[0][-1]

		loop2SubunitLabel = loop2[0][-1]

		#get start index of loop subunit
		loop1StartIndex = ''
		for char in loop1[1]:
			if char.isnumeric():
				loop1StartIndex += char
			else:
				break
		loop1StartIndex = int(loop1StartIndex)

		#get stop index of loop subunit
		loop1StopIndex = ''
		for char in reversed(loop1[1]):
			if char.isnumeric():
				loop1StopIndex += char
			else:
				break
		loop1StopIndex = int(loop1StopIndex[::-1])

		loop1Seq = ''
		for char in loop1[2]:
			if char.isalpha():
				loop1Seq += char

		#get start index of loop subunit
		loop2StartIndex = ''
		for char in loop2[1]:
			if char.isnumeric():
				loop2StartIndex += char
			else:
				break
		loop2StartIndex = int(loop2StartIndex)

		#get stop index of loop subunit
		loop2StopIndex = ''
		for char in reversed(loop2[1]):
			if char.isnumeric():
				loop2StopIndex += char
			else:
				break
		loop2StopIndex = int(loop2StopIndex[::-1])

		loop2Seq = ''
		for char in loop2[2]:
			if char.isalpha():
				loop2Seq += char

		#get index of first base in closing pair
		loop1ClosingPairStart = ''
		for char in loop1[3]:
			if char.isnumeric():
				loop1ClosingPairStart += char
			elif char == ',':
				break
		loop1ClosingPairStart = int(loop1ClosingPairStart)

		#get index of second base in closing pair
		loop1ClosingPairEnd = ''
		for char in reversed(loop1[3]):
			if char.isnumeric():
				loop1ClosingPairEnd += char
			elif char == ',':
				break
		loop1ClosingPairEnd = int(loop1ClosingPairEnd[::-1])

		#get index of first base in closing pair
		loop2ClosingPairStart = ''
		for char in loop2[3]:
			if char.isnumeric():
				loop2ClosingPairStart += char
			elif char == ',':
				break
		loop2ClosingPairStart = int(loop2ClosingPairStart)

		#get index of second base in closing pair
		loop2ClosingPairEnd = ''
		for char in reversed(loop2[3]):
			if char.isnumeric():
				loop2ClosingPairEnd += char
			elif char == ',':
				break
		loop2ClosingPairEnd = int(loop2ClosingPairEnd[::-1])

		#store closing pair as a tuple
		closingPairs = ((loop1[4][0], loop1[4][2]), (loop2[4][2], loop2[4][0]))

		newInnerLoop = InnerLoop(parentLabel, loop1SubunitLabel, loop2SubunitLabel, loop1Seq, loop2Seq, (loop1StartIndex, loop1StopIndex),
								(loop2StartIndex, loop2StopIndex), closingPairs,((loop1ClosingPairStart, loop1ClosingPairEnd), (loop2ClosingPairEnd, loop2ClosingPairStart)))
		self._addInnerLoopToComponentArray(newInnerLoop)
		self.addInnerLoop(parentLabel, newInnerLoop)



	'''
	Function Name: _parseExternalLoopData()
	Description: Internal method used by _loadFile() to parse all the external loop information in the
	structure type file into the StructureType object
	Parameters: externalLoopData - list - list of data from a line in the structure type file describing an external loop
	Return Type: None
	'''
	def _parseExternalLoopData(self, externalLoopData):
		#get label for external loop
		externalLoopLabel = externalLoopData[0]

		#get start index of loop subunit
		startIndex = ''
		for char in externalLoopData[1]:
			if char.isnumeric():
				startIndex += char
			else:
				break
		startIndex = int(startIndex)

		#get stop index of loop subunit
		stopIndex = ''
		for char in reversed(externalLoopData[1]):
			if char.isnumeric():
				stopIndex += char
			else:
				break
		stopIndex = int(stopIndex[::-1])

		#get multiloop sequence
		seq = ''
		for char in externalLoopData[2]:
			if char.isalpha():
				seq += char

		#get index of first base in 5' closing pair
		closingPair5pStart = ''
		for char in externalLoopData[3]:
			if char.isnumeric():
				closingPair5pStart += char
			elif char == ',':
				break
		closingPair5pStart = int(closingPair5pStart)

		#get index of second base in 5' closing pair
		closingPair5pEnd = ''
		for char in reversed(externalLoopData[3]):
			if char.isnumeric():
				closingPair5pEnd += char
			elif char == ',':
				break
		closingPair5pEnd = int(closingPair5pEnd[::-1])

		#store closing pair as a tuple
		closingPair5p = (externalLoopData[4][0], externalLoopData[4][2])

		#get index of first base in 3' closing pair
		closingPair3pStart = ''
		for char in externalLoopData[5]:
			if char.isnumeric():
				closingPair3pStart += char
			elif char == ',':
				break
		closingPair3pStart = int(closingPair3pStart)

		#get index of second base in 3' closing pair
		closingPair3pEnd = ''
		for char in reversed(externalLoopData[5]):
			if char.isnumeric():
				closingPair3pEnd += char
			elif char == ',':
				break
		closingPair3pEnd = int(closingPair3pEnd[::-1])

		#store closing pair as a tuple
		closingPair3p = (externalLoopData[6][0], externalLoopData[6][2])

		newExternalLoop = ExternalLoop(externalLoopLabel, seq, (startIndex, stopIndex), closingPair5p,
									 (closingPair5pStart, closingPair5pEnd), closingPair3p, (closingPair3pStart, closingPair3pEnd))
		self._addExternalLoopToComponentArray(newExternalLoop)
		self.addExternalLoop(externalLoopLabel, newExternalLoop)



	'''
	Function Name: _parseMultiloopData()
	Description: Internal method used by _loadFile() to parse all the multiloop information in the
	structure type file into the StructureType object
	Parameters: multiloopData - list - list of data from a line in the structure type file describing an Multiloop
	Return Type: None
	'''
	def _parseMultiloopData(self, multiloopData):
		#get Multiloop parent Label
		parentLabel = ''
		for char in multiloopData[0]:
			if char == '.':
				break
			else:
				parentLabel += char

		#get inner loop subunit label
		subunitLabel = multiloopData[0][-1]

		#get start index of loop subunit
		startIndex = ''
		for char in multiloopData[1]:
			if char.isnumeric():
				startIndex += char
			else:
				break
		startIndex = int(startIndex)

		#get stop index of loop subunit
		stopIndex = ''
		for char in reversed(multiloopData[1]):
			if char.isnumeric():
				stopIndex += char
			else:
				break
		stopIndex = int(stopIndex[::-1])

		#get multiloop sequence
		seq = ''
		for char in multiloopData[2]:
			if char.isalpha():
				seq += char

		#get index of first base in 5' closing pair
		closingPair5pStart = ''
		for char in multiloopData[3]:
			if char.isnumeric():
				closingPair5pStart += char
			elif char == ',':
				break
		closingPair5pStart = int(closingPair5pStart)

		#get index of second base in 5' closing pair
		closingPair5pEnd = ''
		for char in reversed(multiloopData[3]):
			if char.isnumeric():
				closingPair5pEnd += char
			elif char == ',':
				break
		closingPair5pEnd = int(closingPair5pEnd[::-1])

		#store closing pair as a tuple
		closingPair5p = (multiloopData[4][0], multiloopData[4][2])

		#get index of first base in 3' closing pair
		closingPair3pStart = ''
		for char in multiloopData[5]:
			if char.isnumeric():
				closingPair3pStart += char
			elif char == ',':
				break
		closingPair3pStart = int(closingPair3pStart)

		#get index of second base in 3' closing pair
		closingPair3pEnd = ''
		for char in reversed(multiloopData[5]):
			if char.isnumeric():
				closingPair3pEnd += char
			elif char == ',':
				break
		closingPair3pEnd = int(closingPair3pEnd[::-1])

		#store closing pair as a tuple
		closingPair3p = (multiloopData[6][0], multiloopData[6][2])

		newMultiloop = Multiloop(parentLabel, subunitLabel, seq, (startIndex, stopIndex),
								closingPair5p, (closingPair5pStart, closingPair5pEnd),
								closingPair3p, (closingPair3pStart, closingPair3pEnd))
		self._addMultiloopToComponentArray(newMultiloop)
		self.addMultiloop(parentLabel, subunitLabel, newMultiloop)




	'''
	Function Name: _parseNCBPData()
	Description: Internal method used by _loadFile() to parse all the NCBP information in the
	structure type file into the StructureType object
	Parameters: ncbpData - list - list of data from a line in the structure type file describing an NCBP
	Return Type: None
	'''
	def _parseNCBPData(self, ncbpData):
		#get label for ncbp
		ncbpLabel = ncbpData[0]

		#get index of 5' base
		base1ILoc = int(ncbpData[1])

		#get 5' base
		base1 = ncbpData[2]

		#get index 3' base
		base2ILoc = int(ncbpData[3])

		#get 3' base
		base2 = ncbpData[4]

		#get location of NCBP in other secondary structure
		if ncbpData[5] == '':
			loc = None
		else:
			loc = ncbpData[5]


		newNCBP = NCBP(ncbpLabel, (base1, base2), (base1ILoc, base2ILoc), loc)
		self.addNCBP(ncbpLabel, newNCBP)




	'''
	Function Name: _parseEndData()
	Description: Internal method used by _loadFile() to parse all the End information in the
	structure type file into the StructureType object
	Parameters: endData - list - list of data from a line in the structure type file describing an end
	Return Type: None
	'''
	def _parseEndData(self, endData):
		#get label for End
		endLabel = endData[0]

		#get start index
		startIndex = ''
		for char in endData[1]:
			if char.isnumeric():
				startIndex += char
			else:
				break
		startIndex = int(startIndex)

		#get stop index
		stopIndex = ''
		for char in reversed(endData[1]):
			if char.isnumeric():
				stopIndex += char
			else:
				break
		stopIndex = int(stopIndex[::-1])

		#get sequence
		seq = ''
		for char in endData[2]:
			if char.isalpha():
				seq += char


		newEnd = End(endLabel, seq, (startIndex, stopIndex))
		self._addEndToComponentArray(newEnd)
		self.addEnd(endLabel, newEnd)




	'''
	Function Name: _parsePseudoknotData()
	Description:
	Parameters:
	Return Type:
	'''
	def _parsePsuedoknotData(self, pkData):
		pass



	'''
	Function Name: _parseSegmentData
	Description:
	Parameters:
	Return Type:
	'''
	def _parseSegmentData(self, segData):
		pass


#############################
###### COMPONENT ARRAY ######
#############################

	'''
	Function Name: _addStemToComponentArray(stem)
	Description: Internal method used in _loadFile() that adds a given stem to the component array
	Parameters: (stem) - Stem object - stem to be added to the component array
	Return Type: None
	'''
	def _addStemToComponentArray(self, stem):
		for i in range(stem.sequence5pILoc()[0]-1, stem.sequence5pILoc()[1]):
			self._componentArray[i] = stem.label()

		for i in range(stem.sequence3pILoc()[0]-1, stem.sequence3pILoc()[1]):
			self._componentArray[i] = stem.label()

	'''
	Function Name: _addBulgeToComponentArray(bulge)
	Description:  Internal method used in _loadFile() that adds a given bulge to the component array
	Parameters: (bulge) - Bulge object - bulge to be added to the component array
	Return Type: None
	'''
	def _addBulgeToComponentArray(self, bulge):
		for i in range(bulge.sequenceILoc()[0]-1, bulge.sequenceILoc()[1]):
			self._componentArray[i] = bulge.label()

	'''
	Function Name: _addHairpinToComponentArray(hairpin)
	Description: Internal method used in _loadFile() that adds a hairpin to the component array
	Parameters: (hairpin) - Hairpin object - hairpin to be added to the component array
	Return Type: None
	'''
	def _addHairpinToComponentArray(self, hairpin):
		for i in range(hairpin.sequenceILoc()[0]-1, hairpin.sequenceILoc()[1]):
			self._componentArray[i] = hairpin.label()

	'''
	Function Name: _addEndToComponentArray(end)
	Description: Internal method used in _loadFile() that adds an end to the component array
	Parameters: (end) - End object - end to be added to the component array
	Return Type: None
	'''
	def _addEndToComponentArray(self, end):
		for i in range(end.sequenceILoc()[0]-1, end.sequenceILoc()[1],):
			self._componentArray[i] = end.label()

	'''
	Function Name: _addInnerLoopToComponentArray(innerLoop)
	Description: Internal method used in _loadFile() that adds an inner loop to the component array
	Parameters: (innerLoop) - InnerLoop object - inner loop to be added to the component array
	Return Type:
	'''
	def _addInnerLoopToComponentArray(self, innerLoop):
		for pair in innerLoop.loopsILoc():
			for i in range(pair[0]-1, pair[1]):
				self._componentArray[i] = innerLoop.parentLabel()

	'''
	Function Name: _addExternalLoopToComponentArray(el)
	Description: Internal method used in _loadFile() that adds an external loop to the component array
	Parameters: (el) - ExternalLoop object - External loop to be added to the component array
	Return Type: None
	'''
	def _addExternalLoopToComponentArray(self, el):
		for i in range(el.sequenceILoc()[0]-1, el.sequenceILoc()[1]):
			self._componentArray[i] = el.label()

	'''
	Function Name: _addMultiloopToComponentArray(multiloop)
	Description: Internal method used in _loadFile() that adds an multiloop to the component array
	Parameters: (multiloop) - Multiloop object - multiloop to be added to the component array
	Return Type: None
	'''
	def _addMultiloopToComponentArray(self, multiloop):
		for i in range(multiloop.sequenceILoc()[0]-1, multiloop.sequenceILoc()[1]):
			self._componentArray[i] = multiloop.parentLabel() + '.' + multiloop.subunitLabel()

	'''
	Function Name: componentArray()
	Description: function that returns the _componentArray for the StructureType object
	Parameters: None
	Return Type: numpy array
	'''
	def componentArray(self):
		return self._componentArray






###########################
###### SEQUENCE INFO ######
###########################



	'''
	Function Name: Name()
	Description: Function returns that name of the RNA molecule represented in the .st file
	Parameters: None
	Return Type: str
	'''
	def Name(self):
		return self._name

	'''
	Function Name: Length()
	Description: Function returns the length of the RNA sequence represented in the .st file
	Parameters: None
	Return Type: int
	'''
	def Length(self):
		return self._length

	'''
	Function Name: PageNum()
	Description: Function returns the page number for the RNA molecule represented in the .st file
	Parameters: None
	Return Type: int
	'''
	def PageNum(self):
		return self._pageNum


########################################
###### STRUCTURAL REPRESENTATIONS ######
########################################


	'''
	Function Name: Sequence()
	Description: Function returns the RNA sequence(A,U,G,C) for the RNA molecule represented in the .st file
	Parameters: None
	Return Type: str
	'''
	def Sequence(self):
		return self._sequence

	'''
	Function Name: DotBracket()
	Description: function to get the Dot Bracket notation for the StructureType object
	Parameters: None
	Return Type: str
	'''
	def DotBracket(self):
		return self._DBNotation


	'''
	Function Name: StructureArray()
	Description: Function to get the structure Array for the StructureTyoe object
	Parameters: None
	Return Type: str
	'''
	def StructureArray(self):
		return self._structureArray


	'''
	Function Name: VARNA()
	Description:
	Parameters: None
	Return Type: str
	'''
	def VARNA(self):
		return self._varna



###################
###### STEMS ######
###################

	'''
	Function Name: addStem()
	Description: Function to add a new stem to the StructureType object
	Parameters: (stemLabel) - str - key for stem object in self._stems dictionary
				(newStem) - Stem object - Stem object to be stored at given key in the self._stems dictionary
	Return Value: None
	'''
	def addStem(self, stemLabel, newStem):
		self._stems[stemLabel] = newStem


	'''
	Function Name: stemLabels()
	Description: function to access all the stem labels for the StructureType object
	Parameters: None
	Return Type: list
	'''
	def stemLabels(self):
		return list(self._stems.keys())


	'''
	Function Name: stems()
	Description: function to get all the stem objects in the StructureType object
	Parameters: None
	Return Type: list
	'''
	def stems(self):
		return list(self._stems.values())


	'''
	Function Name: numStems()
	Description: function to get the number of stems in a StructureType object
	Parameters: None
	Return Value: int
	'''
	def numStems(self):
		return len(self._stems)


	'''
	Function Name: getStemByLabel(stemLabel)
	Description: Function to get a particular Stem object based on its label
	Parameters: (stemLabel) - str - label for Stem to be accessed
	Return Value: Stem object
	'''
	def getStemByLabel(self, stemLabel):
		try:
			stem = self._stems[stemLabel]
			return stem
		except KeyError:
			print(f'Stem label: {stemLabel} not found.')
			return None



######################
###### HAIRPINS ######
######################



	'''
	Function Name: addHairpin(label, newHairpin)
	Description: Function to add a new hairpin to the StructureType object
	Parameters: (label) - str - key for Hairpin object in self._hairpins dictionary
				(newHairpin) - Hairpin object - Hairpin object to be stored at the given key in the self._haripins dicitonary
	Return Type: None
	'''
	def addHairpin(self, label, newHairpin):
		self._hairpins[label] = newHairpin


	'''
	Function Name: hairpinLabels()
	Description: Function to access all the hairpin labels in the StructureType object
	Parameters: None
	Return Type: list
	'''
	def hairpinLabels(self):
		return list(self._hairpins.keys())


	'''
	Function Name: hairpins()
	Description: Function to get all Hairpin objects in the StructureType object
	Parameters: None
	Return Type: list
	'''
	def hairpins(self):
		return list(self._hairpins.values())


	'''
	Function Name: numHairpins()
	Description: Function to get the number of hairpins in the StructureType object
	Parameters: None
	Return Type: int
	'''
	def numHairpins(self):
		return len(self._hairpins)

	'''
	Function Name: getHairpinByLabel(hairpinLabel)
	Description: Function to access a particular Hairpin Object based on its label
	Parameters: (hairpinLabel) - str - label for hairpin to accessed
	Return Type: Hairpin Object
	'''
	def getHairpinByLabel(self, hairpinLabel):
		try:
			hairpin = self._hairpins[hairpinLabel]
			return hairpin
		except KeyError:
			print(f'Hairpin label: {hairpinLabel} not found')



####################
###### BULGES ######
####################



	'''
	Function Name: addBulge(bulgeLabel, newBulge)
	Description: Function to add a new Bulge object to the StructureType object
	Parameters: (bulgeLabel) - str - the key value to be used for the new Bulge object
				(newBulge) - Bulge Object - Bulge object to be stored at the given key in the self._bulges dictionary
	Return Type: None
	'''
	def addBulge(self, bulgeLabel, newBulge):
		self._bulges[bulgeLabel] = newBulge


	'''
	Function Name: bulgeLabels()
	Description: Function to get all the bulge labels for the StructureType object
	Parameters: None
	Return Type: list
	'''
	def bulgeLabels(self):
		return list(self._bulges.keys())

	'''
	Function Name: bulges()
	Description: Function to get all the Bulge objects for the StructureType object
	Parameters: None
	Return Type: list
	'''
	def bulges(self):
		return list(self._bulges.values())

	'''
	Function Name: numBulges()
	Description: Function to get the number of bulges in a given StructureType object
	Parameters: None
	Return Type: int
	'''
	def numBulges(self):
		return len(self._bulges)

	'''
	Function Name: getBulgeByLabel(bulgeLabel)
	Description: Function to access a particular Bulge object based on its label
	Parameters: (bulgeLabel) - str - label for Bulge object to be accessed
	Return Type: Bulge object
	'''
	def getBulgeByLabel(self, bulgeLabel):
		try:
			bulge = self._bulges[bulgeLabel]
			return bulge
		except KeyError:
			print(f'Bulge label: {bulgeLabel} not found')
			return None


	'''
	Function Name: bulgeNeighbors(bulgeLabel)
	Description: Function returns the labels of the secondary structures that are directly adjacent to
	the bulge with the given label
	Parameters: (bulgeLabel) - str - label for bulge of interest
	Return Type: tupel(str, str)
	'''
	def bulgeNeighbors(self, bulgeLabel):
		bulgeILoc = self.getBulgeByLabel(bulgeLabel).sequenceILoc() #get tuple with bulge start and stop

		# bulgeILoc = (start, stop)
		# in component array bulge will be at (start-1, stop-1) because of 0 indexing
		# so adjacent features will be at (start-2, stop)
		return (self._componentArray[bulgeILoc[0]-2], self._componentArray[bulgeILoc[1]])




###########################
###### Inner Loops ########
###########################



	'''
	Function Name: addInnerLoop(parentLabel, subunitLabel, newInnerLoop)
	Description: function to add a new InnerLoop object to the StructureType object
	Parameters: (parentLabel) - str - parent key value for the Inner loop to be added
				(newInnerLoop) - InnerLoop object - InnerLoop object to be stored at the the given key
				in the self._innerLoops dictionary
	Return Type: None
	'''
	def addInnerLoop(self, parentLabel, newInnerLoop):
		self._innerLoops[parentLabel] = newInnerLoop


	'''
	Function Name: innerLoopLabels()
	Description: Function to return a list of all the inner loop labels in the StructureType object
	Parameters: None
	Return Type: list
	'''
	def innerLoopLabels(self):
		return list(self._innerLoops.keys())


	'''
	Function Name: innerLoops()
	Description: Function to return a list of the InnerLoop objects in the StructureType object
	Parameters: None
	Return Type: list
	'''
	def innerLoops(self):
		return list(self._innerLoops.values())


	'''
	Function Name: numInnerLoops()
	Description: function to get the number of inner loops in a StructureType object
	Parameters: None
	Return Type: int
	'''
	def numInnerLoops(self):
		return len(self._innerLoops)


	'''
	Function Name: getInnerLoopByLabel(label)
	Description: returns InnerLoop object stored at the given key value
	Parameters: (label) - str - the key value for the inner loop to be accessed
	Return Type: tuple(InnerLoop object, InnerLoop object)
	'''
	def getInnerLoopByLabel(self, label):
		try:
			innerLoop = self._innerLoops[label]
			return innerLoop
		except KeyError:
			print(f'Inner Loop: {label} not found.')
			return None


	'''
	Function Name: getInnerLoopSubunitByLabel(parentLabel, subunitLabel)
	Description:
	Parameters: (parentLabel) - str - parent label for inner loop object
				(subunitLabel) - str - subunit label for the inner loop object
	Return Type: dictionary with the following subunit information

				{
					'label' : ,
					'Sequence' : ,
					'SequenceILoc' :
				}
	'''
	def getInnerLoopSubunit(self, parentLabel, subunitLabel):
		innerLoop = None
		try:
			innerLoop = self._innerLoops[parentLabel]
		except KeyError:
			print(f'Inner Loop: {parentLabel} not found.')


		if subunitLabel == '1':
			subunit = {
				'label' : f'{innerLoop._parentLabel}.1',
				'Sequence' : innerLoop._5pSequence,
				'SequenceILoc' : innerLoop._5pSequenceILoc,
			}
			return subunit
		elif subunitLabel == '2':
			subunit = {
				'label' : f'{innerLoop._parentLabel}.2',
				'Sequence' : innerLoop._3pSequence,
				'SequenceILoc' : innerLoop._3pSequenceILoc,
			}
			return subunit
		else:
			return None





###########################
###### MultiLoops #########
###########################


	'''
	Function Name: addMultiloop(parentLabel, subunitLabel, newMultiloop)
	Description: function to add a new Multiloop object to the StructureType object
	Parameters: (parentLabel) - str - parent multiloop label for the multiloop to be added
				(subunitLabel) - str - subunit label for the multiloop to be added
				(newMultiloop) - Multiloop object - Multiloop object to be added at the given key values
	Return Type: None
	'''
	def addMultiloop(self, parentLabel, subunitLabel, newMultiloop):
		if parentLabel not in self._multiLoops.keys():
			self._multiLoops[parentLabel] = dict()

		self._multiLoops[parentLabel][subunitLabel] = newMultiloop

	'''
	Function Name: numMultiloops()
	Description: Function to get the number of multiloops in a StructureTyoe object
	Parameters: None
	Return Type: int
	'''
	def numMultiloops(self):
		return len(self._multiLoops)


	'''
	Function Name: getMultiloopByLabel(label)
	Description: returns a tuple containing the subunits of the Multiloop stored at the given key value
	Parameters: (label) - str - parent label for the multiloop to be accessed
	Return Type: tuple(Multiloop Object, ... , Multiloop Object)
	'''
	def getMultiloopByLabel(self, label):
		try:
			Multiloop = tuple(self._multiLoops[label].values())
			return Multiloop
		except KeyError:
			print(f'Multiloop: {label} not found.')
			return None

	'''
	Function Name: getMultiloopSubunitByLabel(parentLabel, subunitLabel)
	Description: Function to access a particular subunit for a given multiloop
	Parameters: (parentLabel) - str - parent label for the multiloop to be accessed
				(subunitLabel) - str - subunit label for the multiloop to be accessed
	Return Type: Multiloop object
	'''
	def getMultiloopSubunitByLabel(self, parentLabel, subunitLabel):
		try:
			multi = self._multiLoops[parentLabel][subunitLabel]
			return multi
		except KeyError:
			print(f'Multiloop: {parentLabel}.{subunitLabel} not found')
			return None

	'''
	Function Name: multiloopLabels()
	Description: function to return a list of 2-value tuple where the first value is the parent label for a multiloop
	and the second value is the subunit label for a multioop
	Parameters: None
	Return Type: list of tuples
	'''
	def multilooplabels(self):
		labels = []
		for loop in self._multiLoops.keys():
			for subunit in self._multiLoops[loop].values():
				labels.append((subunit.parentLabel(), subunit.subunitLabel()))

		return labels


	'''
	Function Name: multiloops()
	Description: Function to return a list of tuples, where each tuple contains all the multiloop object subunits composing
	each multiloop in the StructureType object
	Parameters: None
	Return Type: list of tuples
	'''
	def multiloops(self):
		loops = []
		for loop in self._multiLoops.values():
			loops.append(tuple(loop.values()))

		return loops



###############################
###### External Loops #########
###############################


	'''
	Function Name: addExternalLoop(elLabel, newEL)
	Description: function to add a new External Loop object to the StructureType object
	Parameters: (elLabel) - str - the key value to be used for the new ExternalLoop object
				(newEL) - ExternalLoop Object - External loop to be stored at given key value
	Return Type: None
	'''
	def addExternalLoop(self, elLabel, newEL):
		self._externalLoops[elLabel] = newEL

	'''
	Function Name: externalLoopLabels()
	Description: Function to return a list of the external loop labels for a StructureType object
	Parameters: none
	Return Type: list
	'''
	def externalLoopLabels(self):
		return list(self._externalLoops.keys())

	'''
	Function Name: externalLoops()
	Description: function to return a list of all the ExternalLoop Objects in a StructureType object
	Parameters: None
	Return Type: list
	'''
	def externalLoops(self):
		return list(self._externalLoops.values())

	'''
	Function Name: numExternalLoops()
	Description: function to return the number of external loops in a StructureType object
	Parameters: None
	Return Type: int
	'''
	def numExternalLoops(self):
		return len(self._externalLoops)

	'''
	Function Name: getExternalLoopByLabel(elLabel)
	Description: Function to access a particular ExternalLoop Object based on its label
	Parameters: (elLabel) - str - label for the ExternalLoop to be accessed
	Return Type: ExternalLoop object
	'''
	def getExternalLoopByLabel(self, elLabel):
		try:
			el = self._externalLoops[elLabel]
			return el
		except KeyError:
			print(f'External Loop: {elLabel} not found.')
			return None



####################
###### NCBP ########
####################



	'''
	Function Name: addNCBP(ncbpLabel, newNCBP)
	Description: Function to add a new NCBP to the StructureType object
	Parameters: (ncbpLabel) - str - label for the new NCBP
				(newNCBP) - NCBP Object - NCBP object to be added
	Return Type: None
	'''
	def addNCBP(self, ncbpLabel, newNCBP):
		self._ncbp[ncbpLabel] = newNCBP

	'''
	Function Name: ncbpLabels()
	Description: function to return a list of all the ncbp labels in a given StructureType object
	Parameters: None
	Return Type: list
	'''
	def ncbpLabels(self):
		return list(self._ncbp.keys())

	'''
	Function Name: NCBPs()
	Description: function to return a list of all the NCBPs in a StructureType object
	Parameters: None
	Return Type: list
	'''
	def NCBPs(self):
		return list(self._ncbp.values())

	'''
	Function Name: numNCBPs()
	Description: Function to get the number of NCBP's in a given StructureType object
	Parameters: None
	Return Type: int
	'''
	def numNCBPs(self):
		return len(self._ncbp)

	'''
	Function Name: getNCBPByLabel(ncbpLabel)
	Description: Function to get a particular NCBP object based on its label
	Parameters: (ncbpLabel) - str - label for NCBP object to be accessed
	Return Type: NCBP object
	'''
	def getNCBPByLabel(self, ncbpLabel):
		try:
			ncbp = self._ncbp[ncbpLabel]
			return ncbp
		except KeyError:
			print(f'NCBP label: {ncbpLabel} not found.')
			return None



####################
###### ENDS ########
####################



	'''
	Function Name: addEnd(endLabel, newEnd)
	Description: Function to add a new End to the StructureType object
	Parameters: (endLabel) - str - label for new End object
				(newEnd) - End Object - new End object to be added
	Return Type: None
	'''
	def addEnd(self, endLabel, newEnd):
		self._ends[endLabel] = newEnd

	'''
	Function Name: endLabels()
	Description: Function to return a list of all the end labels for the StructureType object
	Parameters: None
	Return Type: list
	'''
	def endLabels(self):
		return list(self._ends.keys())

	'''
	Function Name: ends()
	Description: function to return a list of all the End objects for the StructureType object
	Parameters: None list
	Return Type:
	'''
	def ends(self):
		return list(self._ends.values())

	'''
	Function Name: numEnds()
	Description: function to get the number of ends in StructureType object
	Parameters: None
	Return Type: int
	'''
	def numEnds(self):
		return len(self._ends)

	'''
	Function Name: getEndByLabel(endLabel)
	Description: function to access a particular End Object based on its label
	Parameters: (endLabel) - str - the label for the End object to be accessed
	Return Type: End Object
	'''
	def getEndByLabel(self, endLabel):
		try:
			end = self._ends[endLabel]
			return end
		except KeyError:
			print(f'End: {endLabel} not found.')
			return None



################################
######## OTHER FUNCTIONs #######
################################

	'''
	Function Name: getComponentByLabel(label, subLabel=None)
	Description: More general form of get<secondarystructure>ByLabel(). Allows you to access a given
	structure type component based on its label. Useful when using the componentArray because you may not know which
	type of structure type component you will be accessing if you are(for example) looping through the entire array
	Parameters: (label) - str - label of the feature to be accessed
				(subLabel=None) - str - sublabel for components like innerloops and multiloops, default value is None
				and it will not be used unless specified.
	Return Type: returns a structureTypeComponent
	'''
	def getComponentByLabel(self, label, subLabel=None):
		if label[0] == 'S':
			return self.getStemByLabel(label)
		elif label[0] == 'H':
			return self.getHairpinByLabel(label)
		elif label[0] == 'B':
			return self.getBulgeByLabel(label)
		elif label[0] == 'X':
			return self.getExternalLoopByLabel(label)
		elif label[0] == 'E':
			return self.getEndByLabel(label)
		elif label[0] == 'N':
			return self.getNCBPByLabel(label)

		#need special block to handle internal loops
		if label[0] == 'I': #search is for innerloop
			parent = label[:-2]
			subunit = label[-1]
			return self.getInnerLoopSubunitByLabel(parent, subunit)
			'''if subLabel is not None: #subunit ID provided
					return self.getInnerLoopSubunitByLabel(label, subLabel)
			else: #subunit ID not provided
					return self.getInnerLoopByLabel(label)'''

		#need special block for handeling multiloops
		if label[0] == 'M':
			parent = label[:-2]
			subunit = label[-1]
			return self.getMultiloopSubunitByLabel(parent, subunit)
			'''if subLabel is not None: #subunit ID provided
					return self.getMultiloopSubunitByLabel(label, subLabel)
				else: #subunit ID not provided
					return self.getMultiloopByLabel(label)'''

		#if label is not handled by any of these blocks
		print('Label not found in StructureType object.')
		return None






#################################################################################
#################################################################################
#################################################################################
########################### StructureTypeComponents #############################
#################################################################################
#################################################################################
#################################################################################


'''
Description: This Section defines all the classes that represent RNA secoondary
structures in the structure type file. Each object store relevant information
about the structure and will provide functionality for calculating free energies
associated with each structure.

#NOTE# documentaton/commenting is not complete

'''


'''
## STEM OBJECT ##
the Stem object is used to represent RNA secondary structure stems.

Member variable -- data type -- description:

self._label -- String -- the label for the stem as defined in the structure type file.
self._seq_5p -- String -- the 5' portion of the stem sequence.
self._seq_3p -- String -- the 3' portion of the stem sequence.
self._sequenceLen -- Int -- the length of the stem in number of base pairs.
self._seq_5p_index -- (int, int) -- tuple containing the integer value start and stop
	indices for the 5' portion of the stem sequence.
self._seq_3p_index -- (int, int) -- tuple containing the integer value start and stop
	indices for the 3' portion of the stem sequence.


			5' Sequence

			5' ACGUG 3'
			   |||||
			3' UGCAC 5'

			3' Sequence

'''
class Stem:
	# __init__ method for stem object
	def __init__(self, label, seq_5p, seq_3p, seq_5p_index, seq_3p_index):
		self._label = label #sequence label
		self._seq_5p = seq_5p #5' portion of stem
		self._seq_3p = seq_3p #3' portion of stem
		self._sequence = list(zip(list(self._seq_5p), list(self._seq_3p[::-1])))
		self._sequenceLen = len(seq_5p) #sequence length
		self._seq_5p_index = seq_5p_index #tuple containing start and stop indices of 5' prime portion of stem
		self._seq_3p_index = seq_3p_index #tuple containing start and stop indices of 3' prime portion of stem

	#define string representation of object
	def __str__(self):
		return f'Stem: {self._label}'

	#function returns the label for the stem object
	def label(self):
		return self._label

	#function returns the 5' portion of the stem sequence
	def sequence5p(self):
		return self._seq_5p

	#function returns the 3' portion of the stem sequence
	def sequence3p(self):
		return self._seq_3p

	#function returns the stem sequence as a list of tuples containg base pairs. Ex: [('C','G'), ... , ('A', 'U')]
	def sequence(self):
		return self._sequence

	#function returns the length of the stem
	def sequenceLen(self):
		return self._sequenceLen

	#function returns the start and stop indices of the 5' portion of the stem in a tuple. Ex: (start, stop)
	def sequence5pILoc(self):
		return self._seq_5p_index

	#function returns the start and stop indices of the 3' portion of the stem in a tuple. Ex: (start, stop)
	def sequence3pILoc(self):
		return self._seq_3p_index

 	#function calculates the folding free energy change for the stem
	def energy(self):
		energy = 0
		seq = self.sequence() #get stem as list of tuple base pairs

		#intermolecular initiation
		energy += INTERMOLECULAR_INIT

		#check for symmetry
		if self._seq_5p == self._seq_3p:
			energy += STEM_SYMMETRY_PENALTY

		#check for AU end penalty
		if seq[0] == ('A', 'U') or seq[0] == ('U', 'A'):
			energy += STEM_AU_END_PENALTY
		if seq[-1] == ('A', 'U') or seq[-1] == ('U', 'A'):
			energy = STEM_AU_END_PENALTY

		#sum up watson crick stacking interactions
		for i in range(0, self._sequenceLen-1):
			try:
				energy += StackingEnergies[seq[i]][seq[i+1]]
			except KeyError:
				continue

		return energy





'''
## HAIRPIN OBJECT ##
the Hairpin object is used to represent RNA secondary structure hairpins.

Member variable -- data type -- description:

self._label -- string -- the label for the hairpin as defined by the structure type file.
self._sequence -- string -- the RNA sequence for the hairpin.
self._sequenceLen -- Int -- the length of the hairpin as measured in number of nucleotides.
self._sequenceILoc -- (int, int) -- tuple containing the integer start and stop indices
	for the hairpin.
self._closingPair -- (string, string) -- tuple containing two single character strings. The first character
	corresponds to the 5' base in the closing pair. The second character is the 3' base in the closing pair.
self._closing_index -- (int, int) -- tuple containing two integers. The first integer is the index location of
	the 5' base in the closing pair. The second integer is the index location of the 3'base in the closing pair.
self._pk -- Int -- ???


                  C
                A   G
              G       A
               C     G
                A - U <- _Closing Pair = ('A', 'U')
                C - G
                G - C
                5'  3'


'''
class Hairpin:
	# __init__ method for stem object
	def __init__(self, label, seq, seq_index, closing_pair, closing_index, pk):
		self._label = label
		self._sequence = seq
		self._sequenceLen = len(seq)
		self._sequenceILoc = seq_index
		self._closing_pair = closing_pair
		self._closing_index = closing_index
		self._pk = pk

	#define string representation of object
	def __str__(self):
		return f'Hairpin: {self._label}'

	#Function returns the label for the hairpin object
	def label(self):
		return self._label

	#Function returns the sequence that defines the hairpin structure
	def sequence(self):
		return self._sequence

	#function returns the length of the hairpin
	def sequenceLen(self):
		return self._sequenceLen

	#function returns the start and stop indices of the hairpin as a tuple. Ex: (start, stop)
	def sequenceILoc(self):
		return self._sequenceILoc

	#function returns a tuple that contains the closing pair for the hairpin. Ex: (5' closing base, 3' closing base)
	def closingPair(self):
		return self._closing_pair

	#Function returns the index locations of the closing pair bases as a tuple. Ex: (5' closing index, 3' closing index)
	def closingPairILoc(self):
		return self._closing_index

	#function returns the pseadoknot label for the hairpin if it exists
	def hairpinPK(self):
		return self._pk

	#function to calculate folding free energy of hairpin
	def energy(self):
		pass


'''
BULGE OBJECT
The Bulge object is used to represent the bulge RNA secondary structure.

Member Variable -- Data Type -- Description:
self._label -- string -- The label for the bulge as defined by the structure type file.
self._sequence -- string -- The RNA sequence for the bulge.
self._sequenceLen -- Int -- The length of the bulge as measured in nucleotides.
self._sequenceILoc -- (int, int) -- Tuple containing the integer start and stop indices
	for the RNA sequene that defines the bulge.
self._closingPair5p -- (string, string) -- Tuple containing 2 single character strings. The
	first string the the 5' base in 5' closing pair for the bule. The second character is the
	3' base in the 5' closing pair.
self._closingPair5pILoc -- (int, int) -- Tuple containing 2 integers. The first integer
	is the index of the 5' base in 5' closing pair for the bule. The second integer is the
	index of the 3' base in the 5' closing pair
self._closingPair3p -- (string, string) -- Tuple containing 2 single character strings. The
	first string the the 5' base in 3' closing pair for the bule. The second character is the
	3' base in the 3' closing pair.
self._closingPair3pILoc -- (int, int) -- Tuple containing 2 integers. The first integer
	is the index of the 5' base in 3' closing pair for the bule. The second integer is the
	index of the 3' base in the 3' closing pair
self._pk -- int -- ???



                    C
            5'   AGC UAG   3'
                 ||| |||
            3'   UCG-AUC   5'
                     ^ 3' closing pair = ('U', 'A')
                   ^
                    5' closing pair = ('C', 'G')


'''
class Bulge:
	# __init__ method for bulge object
	def __init__(self, label, seq, seq_index, closingPair5, closingPair5ILoc, closingPair3, closingPair3ILoc, pk):
		self._label = label
		self._sequence = seq
		self._sequenceLen = len(seq)
		self._sequenceILoc = seq_index
		self._closingPair5p = closingPair5
		self._closingPair5pILoc = closingPair5ILoc
		self._closingPair3p = closingPair3
		self._closingPair3pILoc = closingPair3ILoc
		self._pk = pk

	#defines string representation for object
	def __str__(self):
		return f'Bulge: {self._label}'

	#Function returns the label for the bulge object
	def label(self):
		return self._label

	#Function returns the sequence that defines the bulge structure
	def sequence(self):
		return self._sequence

	#Function returns the start and stop indices of the bulge as a tuple. Ex: (start, stop)
	def sequenceILoc(self):
		return self._sequenceILoc

	#Function returns the length of the bulge
	def sequenceLen(self):
		return self._sequenceLen

	#Function returns a tuple containg the 5' closing pair for the bulge
	def closingPair5p(self):
		return self._closingPair5p

	#Function returns a tuple containg the indices of the 5' closing pair for the bulge
	def closingPair5pILoc(self):
		return self._closingPair5pILoc

	#Function returns a tuple containg the 3' closing pair for the bulge
	def closingPair3p(self):
		return self._closingPair3p

	#Function returns a tuple containg the indices of the 3' closing pair for the bulge
	def closingPair3pILoc(self):
		return self._closingPair3pILoc

	#function calculates the folding free energy change for the bulge
	def energy(self):
		#get bulge initiation value
		bulgeInitiation = BulgeInit[self._sequenceLen]

		#Add condition for special c bulge case
		cBulge = 0
		if self._sequence == 'C':
			cBulge = SPECIAL_C_BULGE

		#check for AU end penalty for both closing pairs
		endPenalty = 0
		if self._closingPair3p == ('A', 'U') or self._closingPair3p == ('U', 'A'):
			endPenalty += BULGE_AU_END_PENALTY
		if self._closingPair5p == ('A', 'U') or self._closingPair5p == ('U', 'A'):
			endPenalty += BULGE_AU_END_PENALTY

		#sum energy values and return
		energy = bulgeInitiation + endPenalty + INTERMOLECULAR_INIT + cBulge
		return energy



'''
INNER LOOP
The InnerLoop onject is used to represent the inner loop RNA secondary structure.

Member variable -- data type -- description:
self._parentLabel -- string -- parent label for the 2 inner loop subcomponents
self._5pLabel -- string -- label for the 5' inner loop subcomponent
self._3pLabel -- string -- label for the 3' inner loop subcomponent
self._5pLoop -- string -- sequence that defines the 5' inner loop subcomponent
self._3pLoop -- string -- sequence that defines the 3' inner loop subcomponent
self._loopsLen -- tuple(int, int) -- tuple containing the integer value lengths for the 5' and 3' inner loop subcomponents
self._5pLoopILoc -- tuple(int, int) -- tuple containing the integer start and stop locations for the 5' inner loop subcomponent
self._3pLoopILoc -- tuple(int, int) -- tuple containing the integer start and stop locations for the 3' inner loop subcomponent
self._closingPairs -- tuple((string, string), (string, string)) -- tuple with two nested tuples containing the closing pairs for the inner loop
self._closingPairsILoc -- tuple((int, int), (int, int)) -- tuple with two nested tuples containing the index locations of the closing pairs for the inner loop
'''
class InnerLoop:
	# __init__ method for InnerLoop object
	def __init__(self, pLabel, label5p, label3p,  loop5p, loop3p, loop5pILoc, loop3pILoc, closingPairs, closingPairsILoc):
		self._parentLabel = pLabel
		self._5pLabel = label5p
		self._3pLabel = label3p
		self._5pLoop = loop5p
		self._3pLoop = loop3p
		self._loopsLen = (len(loop5p), len(loop3p))
		self._5pLoopILoc = loop5pILoc
		self._3pLoopILoc = loop3pILoc
		self._closingPairs = closingPairs
		self._closingPairsILoc = closingPairsILoc

	#defines the string representation of the object
	def __str__(self):
		return f'Inner Loop: {self._parentLabel}'

	#Function returns the parent label for the inner loop
	def parentLabel(self):
		return self._parentLabel

	#Function returns a tuple containing the labels for the two inner loop subcomponents
	def subunitLabel(self):
		return (self._5pLabel, self._3pLabel)

	#Function returns a tuple containing the sequences for the two inner loop subcomponents
	def loops(self):
		return (self._5pLoop, self._3pLoop)

	#Function returns a tuple containing the the integer value lengths of the two inner loop components
	def loopsLen(self):
		return self._loopsLen

	#Function returns a tuple that contains two tuples containing the integer start and stop positions of the 5' and 3' inner loop components
	def loopsILoc(self):
		return (self._5pLoopILoc, self._3pLoopILoc)

	#Function returns a tuple that contains two tuples containing the closing base pairs of the inner loop components
	def closingPairs(self):
		return self._closingPairs

	#Function returns a tuple that contains two tuples containing the index locations of the closing base pairs of the inner loop components
	def closingPairsILoc(self):
		return self._closingPairsILoc

	#Function to calculate the free energy for the inner loop
	def energy(self):
		#intermolecular Initiation portion
		intermolecularInit = 4.09

		#check for 1x1 - value taken from imported dicitionary
		if len(self._5pLoop) == 1 and len(self._3pLoop) == 1:
			loopEnergy = InnerLoop_1x1_Energies[self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop]
			return intermolecularInit + loopEnergy
		#check for 1x2 - value taken from imported dicitionary
		elif len(self._5pLoop) == 1 and len(self._3pLoop) == 2:
			loopEnergy = InnerLoop_1x2_Energies[self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop[0], self._3pLoop[1]]
			return intermolecularInit + loopEnergy
		#check for 2x2 - value taken from imported dicitionary
		elif len(self._5pLoop) == 2 and len(self._3pLoop) == 2:
			loops = list(zip(list(self._5pLoop), list(self._3pLoop[::-1])))
			loopEnergy = InnerLoop_2x2_Energies[self._closingPairs[0], self._closingPairs[1], loops[0], loops[1]]
			return intermolecularInit + loopEnergy
		#Other cases
		else:
			#Inner loop initiation energy
			loopLength = len(self._5pLoop) + len(self._3pLoop)

			try: #try to get initiation energy from dictionary
				ilInit = InternalLoopInit[loopLength]
			except: #otherwise calculate value
				ilInit = 2.0 + (1.08 * math.log(float(loopLength)))

			#asymmetry penalty
			asym = abs(len(self._5pLoop) - len(self._3pLoop)) * INNER_LOOP_ASYMMETRY_PENALTY

			#AU / GU Closure penalty
			closingPenalty = 0
			closingPair5p, closingPair3p = self.closingPairs() #get the closing pairs for the inner loop
			if closingPair5p in [('A', 'U'), ('G', 'U'), ('U', 'A'), ('U', 'G')]: #check for penalty condition in 5' closing pair
				closingPenalty += 0.7
			if closingPair3p in [('A', 'U'), ('G', 'U'), ('U', 'A'), ('U', 'G')]: #check for penalty in 3' closing pair
				closingPenalty += 0.7

			#Check for mismatches for base pairs on either end of the inner loop
			mismatch = 0
			#2x3 Inner Loops
			if (len(self._5pLoop) == 2 and len(self._3pLoop == 3)):
				pass
			#3x2 Inner Loops
			elif (len(self._5pLoop) == 3 and len(self._3pLoop == 2)):
				pass
			#other inner loops
			else:
				pass

			#sum energy components and return
			return INTERMOLECULAR_INIT + ilInit + asym + closingPenalty + mismatch



'''
External Loops

Member variable -- data type -- description:
self._label -- string -- label for the external loop secondary structure
self._sequence -- string -- base sequence that defines the external loop
self._sequenceILoc --tuple(int, int) -- tuple containing the integer start and stop locations for the external loop sequence
self._closingPair5p -- tuple(string, string) -- tuple that contains the 5' closing base pair for the external loop
self._closingPair5pILoc -- tuple(int, int) -- tuple containg the integer index locations for the 5' closing pair
self._closingPair3p -- tuple(string, string) -- tuple that contains the 3' closing base pair for the external loop
self._closingPair3pILoc -- tuple(int, int) -- tuple containg the integer index locations for the 3' closing pair
'''
class ExternalLoop:
	#__init__() method for the external loop object
	def __init__(self, label, seq, seqILoc, closingPair5p, closingPair5pILoc, closingPair3p, closingPair3pILoc):
		self._label = label
		self._sequence = seq
		self._sequenceILoc = seqILoc
		self._closingPair5p = closingPair5p
		self._closingPair5pILoc = closingPair5pILoc
		self._closingPair3p = closingPair3p
		self._closingPair3pILoc = closingPair3pILoc

	#Defines the string representation of the external loop
	def __str__(self):
		return f'External Loop: {self._label}'

	#Function returns the label for the external loop
	def label(self):
		return self._label

	#Function returns the sequence that defines the external loop
	def sequence(self):
		return self._sequence

	#Function returns a tuple containing the start and stop index locations for the external loop sequence
	def sequenceILoc(self):
		return self._sequenceILoc


'''
ENDS

Member variable -- data type -- description:
self._label -- string -- label for the end objects
self._sequence -- string -- sequence that defines the end objects
self._sequenceILoc -- tuple(int, int) -- tuple containing the integer start and stop locations for the end object
'''
class End:
	#__init__() method for end object
	def __init__(self, label, seq, seqILoc):
		self._label = label
		self._sequence = seq
		self._sequenceILoc = seqILoc

	#define string representation of end object
	def __str__(self):
		return f'End: {self._label}'

	#Function returns the label for the end object
	def label(self):
		return self._label

	#Function returns the sequence that defines the end object
	def sequence(self):
		return self._sequence

	#Function returns a tuple that contains the integer start and stop index locations for the end object
	def sequenceILoc(self):
		return self._sequenceILoc


'''
NON-CANONICAL BASE PAIRINGS

Member variable -- data type -- description:
self._label -- string -- label for the NCBP objects
self._basePair -- tuple(string, string) -- tuple containing the base pairs that define the NCBP object
self._basePairILoc -- tuple(int, int) -- tuple containing the integer locations of the NCBP
self._parentUnit -- string -- label for the secondary structure that the NCBP is located in
'''
class NCBP:
	#__init__() method for the NCBP object
	def __init__(self, label, basePair, basePairILoc, loc):
		self._label = label
		self._basePair = basePair
		self._basePairILoc = basePairILoc
		self._parentUnit = loc

	#Defines the string representation of the NCBP object
	def __str__(self):
		return f'NCBP: {self._label}'

	#Functions returns the label for the NCBP object
	def label(self):
		return self._label

	#Function returns a tuple containing the two base pairs that define the NCBP object
	def basePairs(self):
		return self._basePair

	#Function returns a tuple containing the integer locations of the base pairs that define the NCBP
	def basePairILoc(self):
		return self._basePairILoc

	#Function returns a string the identifies the secondary structure that the NCBP occurs in
	def parentUnit(self):
		return self._parentUnit



'''
Mulitloops -- UNFINISHED
'''
class Multiloop:
	def __init__(self, parentLabel, subunitLabel, sequence, sequenceILoc, closingPair5p, closingPair5pILoc, closingPair3p, closingPair3pILoc):
		self._parentLabel = parentLabel
		self._subunitLabel = subunitLabel
		self._sequence = sequence
		self._sequenceILoc = sequenceILoc
		self._closingPair5p = closingPair5p
		self._closingPair5pILoc = closingPair5pILoc
		self._closingPair3p = closingPair3p
		self._closingPair3pILoc = closingPair3pILoc

	def __str__(self):
		return f'Multiloop: {self._parentLabel}.{self._subunitLabel}'

	def parentLabel(self):
		return self._parentLabel

	def subunitLabel(self):
		return self._subunitLabel

	def sequence(self):
		return self._sequence

	def sequenceILoc(self):
		return self._sequenceILoc



'''
PSEUDOKNOTS -- UNFINISHED
'''
class Pseudoknot:
	pass
