'''
Filename: StructureComponents.py
Author: Michael Hathaway

Description: The Structure Components module defines individual classes for each of the secondary structures defined in the Structure
Type file. These classes are: Stem, Bulge, Hairpin, InnerLoop, ExternalLoop, MultiLoop, PseudoKnot, End, and NCBP.
'''

## Module Imports ##
import numpy as np
import logging

## Free Energy Parameter Imports ##
from LoopInitiationEnergy import InternalLoopInit, BulgeInit, HairpinInit #initiation parameters for internal loops, bulges, and hairpins
from StackingEnergies import StackingEnergies #Watson-Crick stacking interaction parameters
from InnerLoop_1x1_Energies import InnerLoop_1x1_Energies #Stabilities for 1x1 internal loops
from InnerLoop_1x2_Energies import InnerLoop_1x2_Energies #Stabilities for 1x2 internal loops
from InnerLoop_2x2_Energies import InnerLoop_2x2_Energies #Stabilities for 2x2 internal loops
from InnerLoopMismatches import InnerLoopMismatches_2x3, OtherInnerLoopMismtaches #energy values for 2x3 inner loop mismatches
from StackTerminalMismatches import StackTerminalMismatches #stacking terminal mismatches for Hairpin calculations

## Free Energy Parameter Constants ##
INTERMOLECULAR_INIT = 4.09 #intermolecular initiation value
R = 0.001987204258 #source: https://en.wikipedia.org/wiki/Gas_constant
T = 310.15

#Stems(source: https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html)
STEM_SYMMETRY_PENALTY = 0.43
STEM_AU_END_PENALTY = 0.45

#Inner loops
INNER_LOOP_ASYMMETRY_PENALTY = 0.6

#Bulges
SPECIAL_C_BULGE = -0.9
BULGE_AU_END_PENALTY = 0.45

#Hairpins
HAIRPIN_UU_GA_FIRST_MISMATCH_BONUS = -0.9
HAIRPIN_GG_FIRST_MISMATCH_BONUS = -0.8
HAIRPIN_SPECIAL_GU_CLOSURE = -2.2
HAIRPIN_C3 = 1.5
HAIRPIN_C_LOOP_A = 0.3
HAIRPIN_C_LOOP_B = 1.6

#other Constants
CANONICAL_BASE_PAIRS = [('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G')]

#set logging configuration
logging.basicConfig(filename='./StructureComponents.log', level=logging.WARNING, filemode='w', format='%(process)d - %(levelname)s - %(message)s')


'''
## STEM OBJECT ##
the Stem object is used to represent RNA secondary structure stems.

Member variable -- data type -- description:

self._label -- String -- the label for the stem as defined in the structure type file.
self._sequence5p -- String -- the 5' portion of the stem sequence.
self._sequence3p -- String -- the 3' portion of the stem sequence.
self._sequenceLen -- Int -- the length of the stem in number of base pairs.
self._sequence5p_index -- (int, int) -- tuple containing the integer value start and stop
	indices for the 5' portion of the stem sequence.
self._sequence3p_index -- (int, int) -- tuple containing the integer value start and stop
	indices for the 3' portion of the stem sequence.


			5' Sequence

			5' ACGUG 3'
			   |||||
			3' UGCAC 5'

			3' Sequence

'''
class Stem:
	# __init__ method for stem object
	def __init__(self, label="", sequence5p="", sequence3p="", sequence5p_span=(-1, -1), sequence3p_span=(-1, -1)):
		self._label = label #sequence label
		self._sequence5p = sequence5p #5' portion of stem
		self._sequence3p = sequence3p #3' portion of stem
		self._sequence = list(zip(list(self._sequence5p), list(self._sequence3p[::-1])))
		self._sequenceLen = (len(sequence5p) + len(sequence3p)) // 2 #sequence length
		self._sequence5p_span = sequence5p_span #tuple containing start and stop indices of 5' prime portion of stem
		self._sequence3p_span = sequence3p_span #tuple containing start and stop indices of 3' prime portion of stem

	#define string representation of object
	def __str__(self):
		return f'Stem: {self._label}'

	#Internal method to set the object _sequence member variable as a list of tuples based on the 5' and 3' sequence variables
	def _setSequence(self):
		if len(self._sequence5p) == len(self._sequence3p):
			self._sequence = list(zip(list(self._sequence5p), list(self._sequence3p[::-1])))

	#internal method to update the sequenceLen member variable when the sequence is changed by the user
	def _setSequenceLen(self):
		self._sequenceLen = (len(self._sequence5p) + len(self._sequence3p)) // 2

	#function returns the label for the stem object. Also allows for user to change label of stems
	def label(self, newLabel=None):
		if newLabel :
			self._label = newLabel
		else:
			return self._label

	#function returns the 5' portion of the stem sequence
	def sequence5p(self, newSequence=None):
		if (newSequence):  #if new sequence is provided
			if(len(newSequence) == self._sequenceLen): #check that sequence length matchees other 3' sequences
				self._sequence5p = newSequence
				self._setSequence() #reset the self._sequence variable
			else:
				print('Unable to set new 5\' sequence')
		else:
			return self._sequence5p

	#function returns the 3' portion of the stem sequence
	def sequence3p(self, newSequence=None):
		if (newSequence):  #if new sequence is provided
			if(len(newSequence) == self._sequenceLen): #check that sequence length matchees other 5' sequences
				self._sequence3p = newSequence
				self._setSequence() #reset the self._sequence variable
			else:
				print('Unable to set new 3\' sequence')
		else:
			return self._sequence3p

	#function returns the stem sequence as a list of tuples containg base pairs. Ex: [('C','G'), ... , ('A', 'U')]
	def sequence(self, sequence5p=None, sequence3p=None):
		if(sequence5p and sequence3p):
			if(len(sequence5p) == len(sequence3p)):
				self._sequence5p = sequence5p
				self._sequence3p = sequence3p
				self._setSequence()
				self._setSequenceLen()
			else:
				print('Could not set the stem sequence because the 5\' and 3\' sequences are different lengths.')
		else:
			return self._sequence

	#Function to check if all base pairs in a stem are canonical base pairings
	def canonical(self):
		return all(pair in CANONICAL_BASE_PAIRS for pair in self._sequence)

	#function returns the length of the stem
	def sequenceLen(self):
		return self._sequenceLen

	#function returns a tuple containing two tuples that contain start and stop indices for the 5' and 3' sequence of the stem
	def span(self):
		return (self._sequence5p_span, self._sequence3p_span)

	#function returns the start and stop indices of the 5' portion of the stem in a tuple. Ex: (start, stop)
	def sequence5pSpan(self):
		return self._sequence5p_span

	#function returns the start and stop indices of the 3' portion of the stem in a tuple. Ex: (start, stop)
	def sequence3pSpan(self):
		return self._sequence3p_span

 	#function calculates the folding free energy change for the stem
	def energy(self, strict=True, init=False):
		seq = self.sequence() #get stem as list of tuple base pairs

		#check for symmetry
		symmetry = 0
		if self._sequence5p == self._sequence3p:
			symmetry = STEM_SYMMETRY_PENALTY

		#check for AU end penalty
		endPenalty = 0
		if seq[0] == ('A', 'U') or seq[0] == ('U', 'A'):
			endPenalty += STEM_AU_END_PENALTY
		if seq[-1] == ('A', 'U') or seq[-1] == ('U', 'A'):
			endPenalty = STEM_AU_END_PENALTY

		#sum up watson crick stacking interactions
		stack = 0
		for i in range(0, self._sequenceLen-1):
			try:
				stack += StackingEnergies[seq[i]][seq[i+1]]
			except KeyError:
				logging.warning(f'In energy() function for Stem: {self._label}, Stacking energy not found for {seq[i]} and {seq[i+1]}.')
				if strict: #default strict mode - only calculate energy for stems with all valid parameters
					return None
					break
				else:
					continue

		if(init):
			return INTERMOLECULAR_INIT + symmetry + endPenalty + stack
		else:
			return symmetry + endPenalty + stack





'''
## HAIRPIN OBJECT ##
the Hairpin object is used to represent RNA secondary structure hairpins.

Member variable -- data type -- description:

self._label -- string -- the label for the hairpin as defined by the structure type file.
self._sequence -- string -- the RNA sequence for the hairpin.
self._sequenceLen -- Int -- the length of the hairpin as measured in number of nucleotides.
self._span -- (int, int) -- tuple containing the integer start and stop indices
	for the hairpin.
self._closingPair -- (string, string) -- tuple containing two single character strings. The first character
	corresponds to the 5' base in the closing pair. The second character is the 3' base in the closing pair.
self._closing_span -- (int, int) -- tuple containing two integers. The first integer is the index location of
	the 5' base in the closing pair. The second integer is the index location of the 3'base in the closing pair.
self._pk -- Int -- ???


                  C
                A   G
              G       A
               C     G <- first mismatch
                A - U <- _Closing Pair = ('A', 'U')
                C - G
                G - C
                5'  3'


'''
class Hairpin:
	# __init__ method for stem object
	def __init__(self, label="", sequence="", sequenceSpan=(-1, -1), closingPair=('', ''), closingPairSpan=(-1, -1), pk=None):
		self._label = label
		self._sequence = sequence
		self._sequenceLen = len(sequence)
		self._span = sequenceSpan
		self._closingPair = closingPair
		self._closingPairSpan = closingPairSpan
		self._pk = pk

	#define string representation of object
	def __str__(self):
		return f'Hairpin: {self._label}'

	#Function returns the label for the hairpin object. also allows user to define new label
	def label(self, newLabel=None):
		if newLabel:
			self._label = newLabel
		else:
			return self._label

	#Function returns the sequence that defines the hairpin structure. Also allows user to define new sequence
	def sequence(self, newSequence=None):
		if newSequence:
			self._sequence = newSequence #set new sequence
			self._sequenceLen = len(newSequence) #update sequence length
		else:
			return self._sequence

	#function returns the length of the hairpin
	def sequenceLen(self):
		return self._sequenceLen

	#function returns the start and stop indices of the hairpin as a tuple. Ex: (start, stop)
	def span(self):
		return self._span

	#function returns a tuple that contains the closing pair for the hairpin. Ex: (5' closing base, 3' closing base). Also allows user to define new closing pair
	def closingPair(self, newClose=None):
		if newClose:
			self._closingPair = newClose
		else:
			return self._closingPair

	#Function returns the index locations of the closing pair bases as a tuple. Ex: (5' closing index, 3' closing index)
	def closingPairSpan(self):
		return self._closingPairSpan

	#function returns the pseadoknot label for the hairpin if it exists
	def hairpinPK(self):
		return self._pk

	#function to calculate folding free energy of hairpin
	def energy(self, strict=True):
		#get hairpin initiation term
		if self._sequenceLen in HairpinInit: #try to get from dictionary
			init = HairpinInit[self._sequenceLen]
		else: #otherwise calculate
			init = HairpinInit[9] + (1.75 * R * T * np.log(float(self._sequenceLen/9.0)))

		#get terminal mismatch parameter
		firstMismatch = (self._sequence[0], self._sequence[-1])
		try:
			terminalMismatch = StackTerminalMismatches[self._closingPair][firstMismatch]
		except KeyError:
			logging.warning(f'In energy() function for Hairping: {self._label}, terminal mismatch parameters for closing pair: {self._closingPair} and first mismatch: {firstMismatch} not found in Dictionary.')
			if strict:
				return None #strict mode - only calculate energy for hairpins with valid params
			else:
				terminalMismatch = 0

		#UU/GA first mismatch bonus
		uu_ga_bonus = 0
		if firstMismatch == ('U', 'U') or firstMismatch == ('G', 'A'):
			uu_ga_bonus = HAIRPIN_UU_GA_FIRST_MISMATCH_BONUS

		#GG first mismatch bonus
		gg_bonus = 0
		if firstMismatch == ('G', 'G'):
			gg_bonus = HAIRPIN_GG_FIRST_MISMATCH_BONUS

		#Special GU closure
		gu_closure = 0
		if self._closingPair == ('G', 'U') and firstMismatch == ('G', 'G'):
			gu_closure = HAIRPIN_SPECIAL_GU_CLOSURE

		#All C loop penalty
		c_loop_penalty = 0
		if self._sequence.count('C') == self._sequenceLen:
			c_loop_penalty = (self._sequenceLen * HAIRPIN_C_LOOP_A) + HAIRPIN_C_LOOP_B

		return init + terminalMismatch + uu_ga_bonus + gg_bonus + gu_closure + c_loop_penalty


'''
BULGE OBJECT
The Bulge object is used to represent the bulge RNA secondary structure.

Member Variable -- Data Type -- Description:
self._label -- string -- The label for the bulge as defined by the structure type file.
self._sequence -- string -- The RNA sequence for the bulge.
self._sequenceLen -- Int -- The length of the bulge as measured in nucleotides.
self._span -- (int, int) -- Tuple containing the integer start and stop indices
	for the RNA sequene that defines the bulge.
self._closingPair5p -- (string, string) -- Tuple containing 2 single character strings. The
	first string the the 5' base in 5' closing pair for the bule. The second character is the
	3' base in the 5' closing pair.
self._closingPair5pSpan -- (int, int) -- Tuple containing 2 integers. The first integer
	is the index of the 5' base in 5' closing pair for the bule. The second integer is the
	index of the 3' base in the 5' closing pair
self._closingPair3p -- (string, string) -- Tuple containing 2 single character strings. The
	first string the the 5' base in 3' closing pair for the bule. The second character is the
	3' base in the 3' closing pair.
self._closingPair3pSpan -- (int, int) -- Tuple containing 2 integers. The first integer
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
	def __init__(self, label=None, seq='', sequenceSpan=(-1, -1), closingPair5p=('', ''), closingPair5pSpan=(-1, -1), closingPair3p=('', ''), closingPair3pSpan=(-1, -1), pk=None):
		self._label = label
		self._sequence = seq
		self._sequenceLen = len(seq)
		self._span = sequenceSpan
		self._closingPair5p = closingPair5p
		self._closingPair5pSpan = closingPair5pSpan
		self._closingPair3p = closingPair3p
		self._closingPair3pSpan = closingPair3pSpan
		self._pk = pk

	#defines string representation for object
	def __str__(self):
		return f'Bulge: {self._label}'

	#Function returns the label for the bulge object. Also allows user to define new label
	def label(self, newLabel=None):
		if newLabel:
			self._label = newLabel
		else:
			return self._label

	#Function returns the sequence that defines the bulge structure. Also allows user to define new sequence
	def sequence(self, newSequence=None):
		if newSequence:
			self._sequence = newSequence
			self._sequenceLen = len(newSequence)
		else:
			return self._sequence

	#Function returns the start and stop indices of the bulge as a tuple. Ex: (start, stop)
	def span(self):
		return self._span

	#Function returns the length of the bulge
	def sequenceLen(self):
		return self._sequenceLen

	#Function returns a tuple containg the 5' closing pair for the bulge. Also allows user to define new closing pair
	def closingPair5p(self, newClose=None):
		if newClose:
			self._closingPair5p = newClose
		else:
			return self._closingPair5p

	#Function returns a tuple containg the indices of the 5' closing pair for the bulge
	def closingPair5pSpan(self):
		return self._closingPair5pSpan

	#Function returns a tuple containg the 3' closing pair for the bulge. Also allows user to define new closing pair
	def closingPair3p(self, newClose=None):
		if newClose:
			self._closingPair3p = newClose
		return self._closingPair3p

	#Function returns a tuple containg the indices of the 3' closing pair for the bulge
	def closingPair3pSpan(self):
		return self._closingPair3pSpan

	#function calculates the folding free energy change for the bulge
	def energy(self, strict=True):
		if self._sequenceLen == 1: #bulges of length 1
			#check for special C bulge case
			#special C condition = sequence is all 'C' with at least one adjacent 'C'
			specialC = 0
			if self._sequence == 'C' and (self._closingPair5p[0] == 'C' or self._closingPair3p[0] == 'C'):
				specialC = SPECIAL_C_BULGE

			#get base pair stack
			#base pair stack = the stack of the closing base pairs as if the bulge was not present
			try:
				basePairStack = StackingEnergies[self._closingPair5p][self._closingPair3p]
			except KeyError:
				logging.warning(f'In energy() function for Bulge: {self._label}, No base pair stack found for {self._closingPair5p} and {self._closingPair3p}. Energy Value set to float(\'inf\').')

				if strict:
					return None #strict mode - only calculate energy for bulges with valid params
				else:
					basePairStack = 0


			return BulgeInit[1] + specialC + basePairStack

		else: #bulge of length > 1
			if self._sequenceLen in BulgeInit: #try to get value from dictionary
				return BulgeInit[self._sequenceLen]
			else: #otherwise calculate
				return BulgeInit[6] + (1.75 * R * T * np.log(float(self._sequenceLen/6.0)))



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
self._5pLoopSpan -- tuple(int, int) -- tuple containing the integer start and stop locations for the 5' inner loop subcomponent
self._3pLoopSpan -- tuple(int, int) -- tuple containing the integer start and stop locations for the 3' inner loop subcomponent
self._closingPairs -- tuple((string, string), (string, string)) -- tuple with two nested tuples containing the closing pairs for the inner loop
self._closingPairsSpan -- tuple((int, int), (int, int)) -- tuple with two nested tuples containing the index locations of the closing pairs for the inner loop
self._strict -- bool -- boolean used to control whether energy is calculated strictly
'''
class InnerLoop:
	# __init__ method for InnerLoop object
	def __init__(self, pLabel=None, label5p=None, label3p=None,  loop5p='', loop3p='', loop5pSpan=(-1, -1), loop3pSpan=(-1, -1), closingPairs=(('', ''), ('', '')), closingPairsSpan=((-1, -1), (-1, -1))):
		self._parentLabel = pLabel
		self._5pLabel = label5p
		self._3pLabel = label3p
		self._5pLoop = loop5p
		self._3pLoop = loop3p
		self._loopsLen = (len(loop5p), len(loop3p))
		self._span5p = loop5pSpan
		self._span3p = loop3pSpan
		self._closingPairs = closingPairs
		self._closingPairsSpan = closingPairsSpan
		self._strict = True #used for to control energy function

	#defines the string representation of the object
	def __str__(self):
		return f'Inner Loop: {self._parentLabel}'

	#function to update loop lengths upon change
	def _updateLoopLen(self):
		self._loopsLen = (len(self._5pLoop), len(self._3pLoop))

	#Function returns the parent label for the inner loop. Also allows user to set new label
	def label(self, newLabel=None):
		if newLabel:
			self._parentLabel = newLabel
		return self._parentLabel

	#Function returns a tuple containing the labels for the two inner loop subcomponents
	def subunitLabel(self):
		return (self._5pLabel, self._3pLabel)

	#Function returns a tuple containing the sequences for the two inner loop subcomponents. Also allows user to define new loop sequences.
	def loops(self, loop5p=None, loop3p=None):
		if(loop5p and loop3p):
			self._5pLoop = loop5p
			self._3pLoop = loop3p
			self._updateLoopLen()
		else:
			return (self._5pLoop, self._3pLoop)

	#Function that returns the 5' portion of the inner loop. Also allows user to define 5p portion of the loop
	def loop5p(self, loop=None):
		if(loop):
			self._5pLoop = loop
			self._updateLoopLen()
		else:
			return self._5pLoop

	#Function that returns the 3' portion of the inner loop. Also allows user to define 5p portion of the loop
	def loop3p(self, loop=None):
		if(loop):
			self._3pLoop = loop
			self._updateLoopLen()
		else:
			return self._3pLoop

	#Function returns a tuple containing the the integer value lengths of the two inner loop components
	def loopsLen(self):
		return self._loopsLen

	#Function returns a tuple that contains two tuples containing the integer start and stop positions of the 5' and 3' inner loop components
	def span(self):
		return (self._span5p, self._span3p)

	#Function returns a tuple that contains two tuples containing the closing base pairs of the inner loop components
	def closingPairs(self):
		return self._closingPairs

	#Function returns a tuple that contains two tuples containing the index locations of the closing base pairs of the inner loop components
	def closingPairsSpan(self):
		return self._closingPairsSpan


	'''
	Function Name: _getInnerLoopInitEnergy(self)
	Description: Internal method to get the initiation energy parameter for Inner Loop energy function
	Parameters: None
	Return Type: float
	'''
	def _getInnerLoopInitEnergy(self):
		#get total length of inner loop for initiation parameter calculation
		loopLength = len(self._5pLoop) + len(self._3pLoop)
		if loopLength in InternalLoopInit:#try to get initiation energy from dictionary
			return float(InternalLoopInit[loopLength])
		else: #otherwise calculate value
			return InternalLoopInit[6] + (1.08 * np.log(float(loopLength)/6.0))


	'''
	Function Name: _getInnerLoopAsymmetryEnergy(self)
	Description: Internal method to get asymmetry penalty for inner loop energy function
	Parameters: None
	Return Type: float
	'''
	def _getInnerLoopAsymmetryEnergy(self):
		return abs(len(self._5pLoop) - len(self._3pLoop)) * INNER_LOOP_ASYMMETRY_PENALTY


	'''
	Function Name: _getInnerLoopClosingPenalty(self)
	Description: Internal method to get the AU/GU Closing penalty for InnerLoop energy function
	Parameters: None
	Return Type: float
	'''
	def _getInnerLoopClosingPenalty(self):
		closingPenalty = 0
		endPenaltyPairs = [('A', 'U'), ('G', 'U'), ('U', 'A'), ('U', 'G')] #closing pairs that result in end penalty
		closingPair5p, closingPair3p = self.closingPairs() #get the closing pairs for the inner loop
		if closingPair5p in endPenaltyPairs: #check for penalty condition in 5' closing pair
			closingPenalty += 0.7
		if closingPair3p in endPenaltyPairs: #check for penalty in 3' closing pair
			closingPenalty += 0.7

		return float(closingPenalty)


	'''
	Function Name: _getInnerLoopMismatchEnergy_3x2(self)
	Description: Internal method to get the mismatch energy for a 3x2 InnerLoop
	Parameters: None
	Return Type: float
	'''
	def _getInnerLoopMismatchEnergy_3x2(self):
		loop1, loop2 = self.loops()
		mismatch5p = (loop2[0], loop1[-1])
		mismatch3p = (loop1[0], loop2[-1])

		mismatchEnergy_3x2 = 0
		#check for mismatch condition between 5' closing pair and first mismatch
		if ((self._closingPairs[1][1], self._closingPairs[1][0]), mismatch5p) in InnerLoopMismatches_2x3:
			mismatchEnergy_3x2 += InnerLoopMismatches_2x3[(self._closingPairs[0], mismatch5p)]
		else:
			logging.warning(f'In energy() function for 3x2 InnerLoop: {self._parentLabel}, no mismatch parameter for closing pair: {(self._closingPairs[1][1], self._closingPairs[1][0])} and the 5\' mismatch: {mismatch5p}.')
			if (self._strict):
				return None

		#check for mismatch condition between 3'closing pair and mismatch 2
		if ((self._closingPairs[0][1], self._closingPairs[0][0]), mismatch3p) in InnerLoopMismatches_2x3:
			mismatchEnergy_3x2 += InnerLoopMismatches_2x3[(self._closingPairs[1], mismatch3p)]
		else:
			logging.warning(f'In energy() function for 3x2 InnerLoop: {self._parentLabel}, no mismatch parameter for closing pair: {(self._closingPairs[0][1], self._closingPairs[0][0])} and the 3\' mismatch: {mismatch3p}.')
			if (self._strict):
				return None

		return float(mismatchEnergy_3x2)


	'''
	Function Name: _getInnerLoopMismatchEnergy_2x3(self)
	Description: Internal method to get the mismatch energy for a 2x3 InnerLoop
	Parameters: None
	Return Type: float
	'''
	def _getInnerLoopMismatchEnergy_2x3(self):
		loop1, loop2 = self.loops() #get both loops
		mismatch5p = (loop1[0], loop2[-1]) #get 1st mismatch
		mismatch3p = (loop2[0], loop1[-1]) #get 2nd mismatch

		mismatchEnergy_2x3 = 0
		#check for mismatch condition between 5' closing pair and first mismatch
		if (self._closingPairs[0], mismatch5p) in InnerLoopMismatches_2x3:
			mismatchEnergy_2x3 += InnerLoopMismatches_2x3[(self._closingPairs[0], mismatch5p)]
		else:
			logging.warning(f'In energy() function for 2x3 InnerLoop: {self._parentLabel}, no mismatch parameter for closing pair: {self._closingPairs[0]} and the 5\' mismatch: {mismatch5p}.')
			if (self._strict):
				return None

		#check for mismatch condition between 3'closing pair and mismatch 2
		if ((self._closingPairs[1][1], self._closingPairs[1][0]), mismatch3p) in InnerLoopMismatches_2x3:
			mismatchEnergy_2x3 += InnerLoopMismatches_2x3[(self._closingPairs[1], mismatch3p)]
		else:
			logging.warning(f'In energy() function for 2x3 InnerLoop: {self._parentLabel}, no mismatch parameter for closing pair: {(self._closingPairs[1][1], self._closingPairs[1][0])} and the 3\' mismatch: {mismatch3p}.')
			if (self._strict):
				return None

		return float(mismatchEnergy_2x3)


	'''
	Function Name: _getInnerLoopMismatchEnergy_Other(self)
	Description: Internal method to get the inner loop mismtach energy for other inner loops
	Parameters: None
	Return Type: float
	'''
	def _getInnerLoopMismatchEnergy_Other(self):
		loop1, loop2 = self.loops() #get both loops
		mismatch5p = (loop1[0], loop2[-1]) #get 1st mismatch
		mismatch3p = (loop1[-1], loop2[0]) #get 2nd mismatch

		mismatchEnergy_Other = 0
		#check for mismatch 1 for condition
		if mismatch5p in OtherInnerLoopMismtaches:
			mismatchEnergy_Other += OtherInnerLoopMismtaches[mismatch5p]
		elif (self._strict):
			return None

		#check mismatch 2 for condition
		if mismatch3p in OtherInnerLoopMismtaches:
			mismatchEnergy_Other += OtherInnerLoopMismtaches[mismatch3p]
		elif (self._strict):
			return None

		return float(mismatchEnergy_Other)


	'''
	Function Name: _getInnerLoopMismtachEnergy(self)
	Description: Internal method to get the mismatch energy for an inner loop
	Parameters: None
	Return Type: float
	'''
	def _getInnerLoopMismtachEnergy(self):
		#Check for mismatches for base pairs on either end of the inner loop
		mismatchEnergy = 0
		#1 x (n-1) Inner Loops
		loopLength = len(self._5pLoop) + len(self._3pLoop)
		if (len(self._5pLoop) == 1 and len(self._3pLoop) == loopLength-1) or (len(self._5pLoop) == loopLength-1 and len(self._3pLoop) == 1):
			return 0.0 #Mismatch energy is 0, so we dont need to do anything

		#2x3 Inner Loop mismatches
		elif (len(self._5pLoop) == 2 and len(self._3pLoop) == 3):
			return self._getInnerLoopMismatchEnergy_2x3()

		#3x2 inner loop mismatches
		elif (len(self._5pLoop) == 3 and len(self._3pLoop) == 2):
			return self._getInnerLoopMismatchEnergy_3x2()

		#other inner loops
		else:
			return self._getInnerLoopMismatchEnergy_Other()


	'''
	Function Name: _calcEnergy(self)
	Description: Internal method  to calculate the energy for inner loops whose energies are not stored in the imported dictionaries
	Parameters: None
	Return Type: float
	'''
	def _calcEnergy(self):
		#get InnerLoop initiation parameter
		ilInit = self._getInnerLoopInitEnergy()
		if(ilInit is None): #check that parameter is present
			return None

		#asymmetry penalty
		asym = self._getInnerLoopAsymmetryEnergy()
		if(asym is None):#check that parameter is present
			return None

		#AU / GU Closure penalty
		closingPenalty = self._getInnerLoopClosingPenalty()
		if(closingPenalty is None):#check that parameter is present
			return None

		#get mismtach energy
		mismatchEnergy = self._getInnerLoopMismtachEnergy()
		if(mismatchEnergy is None):#check that parameter is present
			return None

		#sum energy components and return
		return ilInit + asym + closingPenalty + mismatchEnergy


	'''
	Function Name: energy(self)
	Description: Function to get the free energy for the inner loop object
	Parameters: None
	Return Type: float
	'''
	def energy(self, strict=True):
		#set mode for energy calculations
		self._strict = strict

		#check for 1x1 - value taken from imported dicitionary
		if len(self._5pLoop) == 1 and len(self._3pLoop) == 1:
			if (self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop) in InnerLoop_1x1_Energies: #check if key in dictionary
				loopEnergy = InnerLoop_1x1_Energies[(self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop)]
				return loopEnergy
			else: #otherwise calculate energy
				logging.warning(f'Inner Loop: {self._parentLabel}, loop is 1x1, but energy parameters is not present in InnerLoop_1x1_Energies dicitonary. Energy value calculated using _calcEnergy() function.')
				if(self._strict):
					return None
				else:
					return self._calcEnergy()

		#check for 1x2 - value taken from imported dicitionary
		elif len(self._5pLoop) == 1 and len(self._3pLoop) == 2:
			if (self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop[1], self._3pLoop[0]) in InnerLoop_1x2_Energies: #check if key in dictionary
				loopEnergy = InnerLoop_1x2_Energies[(self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop[1], self._3pLoop[0])]
				return loopEnergy
			else: #otherwise calculate energy
				logging.warning(f'Inner Loop: {self._parentLabel}, loop is 1x2, but energy parameters is not present in InnerLoop_1x1_Energies dicitonary. Energy value calculated using _calcEnergy() function.')
				if(self._strict):
					return None
				else:
					return self._calcEnergy()

		#check for 2x1 case - value taken from dicitonary
		elif len(self._5pLoop) == 2 and len(self._3pLoop) == 1:
			if ((self._closingPairs[1][1], self._closingPairs[1][0]), (self._closingPairs[0][1], self._closingPairs[0][0]), self._3pLoop, self._5pLoop[1], self._5pLoop[0]) in InnerLoop_1x2_Energies: #check if key in dictionary
				loopEnergy = InnerLoop_1x2_Energies[((self._closingPairs[1][1], self._closingPairs[1][0]), (self._closingPairs[0][1], self._closingPairs[0][0]), self._3pLoop, self._5pLoop[1], self._5pLoop[0])]
				return loopEnergy
			else: #otherwise calculate energy
				logging.warning(f'Inner Loop: {self._parentLabel}, loop is 2x1, but energy parameters is not present in InnerLoop_1x1_Energies dicitonary. Energy value calculated using _calcEnergy() function.')
				if(self._strict):
					return None
				else:
					return self._calcEnergy()

		#check for 2x2 - value taken from imported dicitionary
		elif len(self._5pLoop) == 2 and len(self._3pLoop) == 2:
			loops = list(zip(list(self._5pLoop), list(self._3pLoop[::-1]))) #convert loop sequences to proper format for dictionary
			if (self._closingPairs[0], self._closingPairs[1], loops[0], loops[1]) in InnerLoop_2x2_Energies: #check if key in dictionary
				loopEnergy = InnerLoop_2x2_Energies[(self._closingPairs[0], self._closingPairs[1], loops[0], loops[1])]
				return loopEnergy
			else: #otherwise calculate energy
				logging.warning(f'Inner Loop: {self._parentLabel}, loop is 2x2, but energy parameters is not present in InnerLoop_1x1_Energies dicitonary. Energy value calculated using _calcEnergy() function.')
				if(self._strict):
					return None
				else:
					return self._calcEnergy()

		#Other cases need to be calculated
		else:
			return self._calcEnergy()



'''
External Loops

Member variable -- data type -- description:
self._label -- string -- label for the external loop secondary structure
self._sequence -- string -- base sequence that defines the external loop
self._span --tuple(int, int) -- tuple containing the integer start and stop locations for the external loop sequence
self._closingPair5p -- tuple(string, string) -- tuple that contains the 5' closing base pair for the external loop
self._closingPair5pSpan -- tuple(int, int) -- tuple containg the integer index locations for the 5' closing pair
self._closingPair3p -- tuple(string, string) -- tuple that contains the 3' closing base pair for the external loop
self._closingPair3pSpan -- tuple(int, int) -- tuple containg the integer index locations for the 3' closing pair
'''
class ExternalLoop:
	#__init__() method for the external loop object
	def __init__(self, label, seq, seqSpan, closingPair5p, closingPair5pSpan, closingPair3p, closingPair3pSpan):
		self._label = label
		self._sequence = seq
		self._span = seqSpan
		self._closingPair5p = closingPair5p
		self._closingPair5pSpan = closingPair5pSpan
		self._closingPair3p = closingPair3p
		self._closingPair3pSpan = closingPair3pSpan

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
	def span(self):
		return self._span


'''
ENDS

Member variable -- data type -- description:
self._label -- string -- label for the end objects
self._sequence -- string -- sequence that defines the end objects
self._span -- tuple(int, int) -- tuple containing the integer start and stop locations for the end object
'''
class End:
	#__init__() method for end object
	def __init__(self, label, sequence, span):
		self._label = label
		self._sequence = sequence
		self._span = span

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
	def span(self):
		return self._span


'''
NON-CANONICAL BASE PAIRINGS

Member variable -- data type -- description:
self._label -- string -- label for the NCBP objects
self._basePair -- tuple(string, string) -- tuple containing the base pairs that define the NCBP object
self._basePairSpan -- tuple(int, int) -- tuple containing the integer locations of the NCBP
self._parentUnit -- string -- label for the secondary structure that the NCBP is located in
'''
class NCBP:
	#__init__() method for the NCBP object
	def __init__(self, label, basePair, basePairSpan, loc):
		self._label = label
		self._basePair = basePair
		self._basePairSpan = basePairSpan
		self._parentUnit = loc

	#Defines the string representation of the NCBP object
	def __str__(self):
		return f'NCBP: {self._label}'

	#Functions returns the label for the NCBP object
	def label(self):
		return self._label

	#Function returns a tuple containing the two base pairs that define the NCBP object
	def sequence(self):
		return self._basePair

	#Function returns a tuple containing the integer locations of the base pairs that define the NCBP
	def Span(self):
		return self._basePairSpan

	#Function returns a string the identifies the secondary structure that the NCBP occurs in
	def parentUnit(self):
		return self._parentUnit



'''
Mulitloops -- UNFINISHED
'''
class MultiLoop:
	def __init__(self, parentLabel, subunitLabel, sequence, span, closingPair5p, closingPair5pSpan, closingPair3p, closingPair3pSpan):
		self._parentLabel = parentLabel
		self._subunitLabel = subunitLabel
		self._sequence = sequence
		self._span = span
		self._closingPair5p = closingPair5p
		self._closingPair5pSpan = closingPair5pSpan
		self._closingPair3p = closingPair3p
		self._closingPair3pSpan = closingPair3pSpan

	def __str__(self):
		return f'Multiloop: {self._parentLabel}.{self._subunitLabel}'

	def label(self):
		return self._parentLabel

	def subunitLabel(self):
		return self._subunitLabel

	def sequence(self):
		return self._sequence

	def span(self):
		return self._span



'''
PSEUDOKNOTS -- UNFINISHED
'''
class PseudoKnot:
	pass
