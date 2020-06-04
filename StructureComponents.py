'''
Filename: StructureComponents.py
Author: Michael Hathaway

Description: The Structure Components module defines individual classes for each of the secondary structures defined in the Structure
Type file. These classes are: Stem, Bulge, Hairpin, InternalLoop, ExternalLoop, MultiLoop, PseudoKnot, End, and NCBP.
'''

## Module Imports ##
import numpy as np
import logging

## Free Energy Parameter Imports ##
from parameters.LoopInitiationEnergy import InternalLoopInit, BulgeInit, HairpinInit #initiation parameters for internal loops, bulges, and hairpins
from parameters.StackingEnergies import StackingEnergies #Watson-Crick stacking interaction parameters
from parameters.InnerLoop_1x1_Energies import InnerLoop_1x1_Energies #Stabilities for 1x1 internal loops
from parameters.InnerLoop_1x2_Energies import InnerLoop_1x2_Energies #Stabilities for 1x2 internal loops
from parameters.InnerLoop_2x2_Energies import InnerLoop_2x2_Energies #Stabilities for 2x2 internal loops
from parameters.InnerLoopMismatches import InnerLoopMismatches_2x3, OtherInnerLoopMismtaches #energy values for 2x3 inner loop mismatches
from parameters.StackTerminalMismatches import StackTerminalMismatches #stacking terminal mismatches for Hairpin calculations
from parameters.SpecialHairpins import SpecialHairpins #special case hairpins with precalculated energies

## Free Energy Parameter Constants ##
R = 0.001987204258 #source: https://en.wikipedia.org/wiki/Gas_constant
T = 310.15

#Stems(source: https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html)
INTERMOLECULAR_INIT = 4.09 #intermolecular initiation value
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

'''
-- set logging configuration --
Logging file will be used to record errors associated with the StructureComponent energy functions. These errors
are usually related to missing parameters for the energy calculation.
'''
logging.basicConfig(filename='./StructureComponents.log', level=logging.WARNING, filemode='a', format='%(process)d - %(levelname)s - %(message)s')


'''
## STEM OBJECT ##
the Stem object is used to represent RNA secondary structure stems.

Member variable -- data type -- description:
self._label -- String -- the label for the stem as defined in the structure type file.
self._sequence5p -- String -- the 5' portion of the stem sequence.
self._sequence3p -- String -- the 3' portion of the stem sequence.
self._sequenceLen -- Int -- the length of the stem in number of base pairs.
self._sequence5p_index -- (int, int) -- tuple containing the integer value start and stop indices for the 5' portion of the stem sequence.
self._sequence3p_index -- (int, int) -- tuple containing the integer value start and stop indices for the 3' portion of the stem sequence.
self._neighbor5p -- str -- label for 5' neighbor in Structure object
self._neighbor3p -- str -- label for 5' neighbor in Structure object


            5' Sequence

            5' ACGUG 3'
               |||||
            3' UGCAC 5'

            3' Sequence

'''
class Stem:
    # __init__ method for stem object
    def __init__(self, label="", sequence5p="", sequence3p="", sequence5pSpan=(-1, -1), sequence3pSpan=(-1, -1), neighbor5p=('', ''), neighbor3p=('', ''), adjacentBulges=(False, False)):
        self._label = label #sequence label
        self._sequence5p = sequence5p #5' portion of stem
        self._sequence3p = sequence3p #3' portion of stem
        self._sequence = list(zip(list(self._sequence5p), list(self._sequence3p[::-1])))
        self._sequenceLen = (len(sequence5p) + len(sequence3p)) // 2 #sequence length
        self._sequence5pSpan = sequence5pSpan #tuple containing start and stop indices of 5' prime portion of stem
        self._sequence3pSpan = sequence3pSpan #tuple containing start and stop indices of 3' prime portion of stem
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p
        self._adjacentBulges = adjacentBulges

    ###
    ### Internal Methods
    ###

    #define string representation of object
    def __str__(self):
        return f'Stem: {self._label}'

    #define len function operation for Stem objects
    def __len__(self):
        return self._sequenceLen

    #Internal method to set the object _sequence member variable as a list of tuples based on the 5' and 3' sequence variables
    def _setSequence(self):
        if len(self._sequence5p) == len(self._sequence3p):
            self._sequence = list(zip(list(self._sequence5p), list(self._sequence3p[::-1])))

    #internal method to update the sequenceLen member variable when the sequence is changed by the user
    def _setSequenceLen(self):
        self._sequenceLen = (len(self._sequence5p) + len(self._sequence3p)) // 2

    #internal method used during Structure object parsing to track if stem is next to length=1 bulges
    def _addAdjacentBulgeBoolean(self, bulge5p, bulge3p):
        self._adjacentBulges = (bulge5p, bulge3p)

    #Internal method that returns tuple containg booleans for whether or not the stem is adjacent to length=1 bulges
    def _adjacentBulgeBoolean(self):
        return self._adjacentBulges

    #internal method to set the 5' and 3' neighbors for a stem
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    ### User Accesible Methods
    ###

    '''
    Function: Stem.label()
    Description: function returns the label for the stem object. Also allows for user to change label of stems
    Parameters:
            (newLabel=None) -- str -- new label that to update the stem._label member variable
    Return Value:
            str - label for the Stem object
    '''
    def label(self, newLabel=None):
        if newLabel :
            self._label = newLabel
        else:
            return self._label


    '''
    Function: Stem.sequence5p()
    Description: function returns the 5' portion of the stem sequence
    Parameters:
            (newSequence) -- str -- new RNA sequence to define the 5' sequence of the stem.
    Return Value:
            str - The current 5' sequence for the Stem object
    '''
    def sequence5p(self, newSequence=None):
        if (newSequence):  #if new sequence is provided
            if(len(newSequence) == self._sequenceLen): #check that sequence length matchees other 3' sequences
                self._sequence5p = newSequence
                self._setSequence() #reset the self._sequence variable
            else:
                print('Unable to set new 5\' sequence')
        else:
            return self._sequence5p


    '''
    Function: Stem.sequence5p()
    Description: function returns the 3' portion of the stem sequence
    Parameters:
            (newSequence) -- str -- new RNA sequence to define the 3' sequence of the stem.
    Return Value:
            str - The current 3' sequence for the Stem object
    '''
    def sequence3p(self, newSequence=None):
        if (newSequence):  #if new sequence is provided
            if(len(newSequence) == self._sequenceLen): #check that sequence length matchees other 5' sequences
                self._sequence3p = newSequence
                self._setSequence() #reset the self._sequence variable
            else:
                print('Unable to set new 3\' sequence')
        else:
            return self._sequence3p


    '''
    Function: Stem.sequence()
    Description: function returns the stem sequence as a list of tuples containg base pairs. Ex: [('C','G'), ... , ('A', 'U')]
    Parameters:
            (sequence5p=None) -- str -- new RNA sequence to define the 5' sequence of the stem.
            (sequence3p=None) -- str -- new RNA sequence to define the 3' sequence of the stem.
    Return Value:
            list - list of tuples representing the base pair sequence of the stem.
    '''
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


    '''
    Function: Stem.sequenceLen()
    Description: Function returns the length of the Stem object
    Parameters: None
    Return Value:
            int - the integer value length of the stem
    '''
    def sequenceLen(self):
        return self._sequenceLen


    '''
    Function: Stem.span()
    Description: function returns a tuple containing two tuples that contain start and stop indices for the 5' and 3' sequence of the stem
    Parameters: None
    Return Value:
            ((int, int), (int, int)) - a tuple containing the tuple(int, int) start and stop positions for the 5' and 3' stem sequences
    '''
    def span(self):
        return (self._sequence5pSpan, self._sequence3pSpan)


    '''
    Function: Stem.sequence5pSpan()
    Description: function returns the start and stop indices of the 5' portion of the stem in a tuple. Ex: (start, stop)
    Parameters: None
    Return Value:
            (int, int) - a tuple containing the integer value start and stop positions of the 5' portion of the stem
    '''
    def sequence5pSpan(self):
        return self._sequence5pSpan


    '''
    Function: Stem.sequence3pSpan()
    Description: function returns the start and stop indices of the 3' portion of the stem in a tuple. Ex: (start, stop)
    Parameters: None
    Return Value:
            (int, int) - a tuple containing the integer value start and stop positions of the 3' portion of the stem
    '''
    def sequence3pSpan(self):
        return self._sequence3pSpan


    '''
    Function: Stem.neighbors()
    Description: Function returns a tuple containing the labels for the 5' and 3' neighbors of the stem
    Parameters: None
    Return Value:
            ((str, str), (str, str)) - tuple containging tuples containing the string value labels for the neighbors structures of the 5' and 3' portions of the stem sequennce
    '''
    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)


    '''
    Function: Stem.cannonical()
    Description: Function to check if all base pairs in a stem are canonical base pairings
    Parameters: None
    Return Value:
            bool - true or false as to whether or not the stem contains all cannonical base pairings
    '''
    def canonical(self):
        return (self._sequenceLen > 1 and all(pair in CANONICAL_BASE_PAIRS for pair in self._sequence))


    '''
    Function: Stem.energy()
    Description: function calculates the folding free energy change for the stem
    Parameters:
            (strict=True) -- bool -- when true, energy values will only be calculated for cannonical stems/stems with all present energy parameters
            (init=False) -- bool -- when true, the 4.09 Kcal/mol initiation value is inlcuded in energy calculations.
    Return Value:
            float - the calculated energy value for the given stem
    '''
    def energy(self, strict=True, init=False):
        if(self._sequenceLen == 1):
            logging.warning(f'In energy() function for Stem: {self._label}, cannot calculate energy for stem of length 1.')
            return None

        seq = self.sequence() #get stem as list of tuple base pairs

        #check for symmetry
        symmetry = 0
        if self._sequence5p == self._sequence3p:
            symmetry = STEM_SYMMETRY_PENALTY

        #check for AU end penalty
        endPenalty = 0
        if (seq[0] == ('A', 'U') or seq[0] == ('U', 'A') or seq[0] == ('G', 'U') or seq[0] == ('U', 'G')) and (self._adjacentBulgeBoolean()[0] == False):
            endPenalty += STEM_AU_END_PENALTY
        if (seq[-1] == ('A', 'U') or seq[-1] == ('U', 'A') or seq[-1] == ('G', 'U') or seq[-1] == ('U', 'G')) and (self._adjacentBulgeBoolean()[1] == False):
            endPenalty += STEM_AU_END_PENALTY

        #sum up watson crick stacking interactions
        stack = 0
        for i in range(0, self._sequenceLen-1):
            try:
                stack += StackingEnergies[seq[i]][seq[i+1]]
            except KeyError:
                logging.warning(f'In energy() function for Stem: {self._label}, Stacking energy not found for {seq[i]} and {seq[i+1]}.')
                if strict: #default strict mode - only calculate energy for stems with all valid parameters
                    return None
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
self._span -- (int, int) -- tuple containing the integer start and stop indices for the hairpin.
self._closingPair -- (string, string) -- tuple containing two single character strings. The first character corresponds to the 5' base in the closing pair. The second character is the 3' base in the closing pair.
self._closing_span -- (int, int) -- tuple containing two integers. The first integer is the index location of the 5' base in the closing pair. The second integer is the index location of the 3'base in the closing pair.
self._pk -- Int -- The pseudoknot the hairpin is a part of, if any(default value is None)
self._neighbor -- str -- label for the neighboring stem to the hairpin


                  C
                A   G
              G       A
               C     G <- first mismatch = ('C', 'G')
                A - U <- _Closing Pair = ('A', 'U')
                C - G
                G - C
                5'  3'


'''
class Hairpin:
    # __init__ method for stem object
    def __init__(self, label="", sequence="", sequenceSpan=(-1, -1), closingPair=('', ''), closingPairSpan=(-1, -1), pk=None, neighbors=('', '')):
        self._label = label
        self._sequence = sequence
        self._sequenceLen = len(sequence)
        self._span = sequenceSpan
        self._closingPair = closingPair
        self._closingPairSpan = closingPairSpan
        self._pk = pk
        self._neighbors = neighbors


    ###
    ### Internal Methods
    ###

    #define string representation of object
    def __str__(self):
        return f'Hairpin: {self._label}'

    #define len function operation for Hairpin objects
    def __len__(self):
        return self._sequenceLen

    #Internal function to add hairpin neighbors to the object
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbors = (neighbor5p, neighbor3p)


    ###
    ### User Accessible Methods
    ###

    '''
    Function: Hairpin.label()
    Description: Function returns the label for the hairpin object. also allows user to define new label
    Parameters:
            (newLabel=None) -- str -- new label to indentify the hairpin object
    Return Value:
            str - the current label for the given hairpin object
    '''
    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label


    '''
    Function: Hairpin.sequence()
    Description: Function returns the sequence that defines the hairpin structure. Also allows user to define new sequence
    Parameters:
            (newSequence=None) -- str -- new sequence to define the hairpin loop
    Return Value:
            str - the current sequence that defines the hairpin loop
    '''
    def sequence(self, newSequence=None):
        if newSequence:
            self._sequence = newSequence #set new sequence
            self._sequenceLen = len(newSequence) #update sequence length
        else:
            return self._sequence


    '''
    Function: Hairpin.SequenceLen()
    Description: function returns the length of the hairpin
    Parameters: None
    Return Value:
            int - the length of the hairpin sequence
    '''
    def sequenceLen(self):
        return self._sequenceLen


    '''
    Function: Hairpin.span()
    Description: function returns the start and stop indices of the hairpin as a tuple. Ex: (start, stop)
    Parameters: None
    Return Value:
            (int, int) - tuple containing the integer value start and stop positions of the hairpin loop
    '''
    def span(self):
        return self._span


    '''
    Function: Hairpin.closingPair()
    Description: function returns a tuple that contains the closing pair for the hairpin. Ex: (5' closing base, 3' closing base). Also allows user to define new closing pair
    Parameters:
            (newClose) -- (str, str) -- new closing pair for the Hairpin object
    Return Value:
            (str, str) - tuple containing the 5' and 3' closing bases for the Hairpin
    '''
    def closingPair(self, newClose=None):
        if newClose:
            try:
                if newClose[0] and newClose[1]:
                    self._closingPair = newClose
            except:
                print('Please provide and tuple with the opening and closing base pairs for the hairpin.')
        else:
            return self._closingPair


    '''
    Function: Hairpin.closingPairSpan()
    Description: Function returns the index locations of the closing pair bases as a tuple. Ex: (5' closing index, 3' closing index)
    Parameters: None
    Return Value:
            (int, int) - tuple containing the integer value locations of the closing base pairs for the Hairpin
    '''
    def closingPairSpan(self):
        return self._closingPairSpan


    '''
    Function: Hairpin.hairpinPK()
    Description: function returns the pseadoknot label for the hairpin if it exists
    Parameters: None
    Return Value:
            (int) - the pseudoknot that the hairpin is a part of if it exist. Will return none if not part of pseudoknot
    '''
    def hairpinPK(self):
        return self._pk


    '''
    Function: Hairpin.neighbors()
    Description: funtion to get the labels for the StructureComponents adjacent to the hairpin
    Parameters: None
    Return Value:
            (str, str) - tuple containing the labels for the neighboring secondary structures of the 5' and 3' ends of the sequence. Note: for hairpins the neighbors will always be the same.
    '''
    def neighbors(self):
        return self._neighbors


    '''
    Function: Hairpin.cannonical()
    Description: Function to check if the correct parameters are available to calculate the energy of the hairpin
    Parameters: None
    Return Value:
            bool - function return True if all necessary energy values are present for the Hairpin
    '''
    def canonical(self):
        firstMismatch = (self._sequence[0], self._sequence[-1])
        if(self._closingPair not in StackTerminalMismatches) or (firstMismatch not in StackTerminalMismatches[self._closingPair]):
            return False
        elif self._sequenceLen < 3:
            return False
        else:
            return True


    '''
    Function: Hairpin.energy()
    Description: function to calculate folding free energy of hairpin
    Parameters:
            (strict=True) -- bool -- when True, the function will only calculate the energy of the molecule valid energy parameters are present.
    Return Value:
            float - the calculated energy for the hairpin
    '''
    def energy(self, strict=True):
        #check that hairpin is at least 3 nucleotides long
        if self._sequenceLen < 3:
            logging.warning(f'In energy() function for Hairpin: {self._label}, hairpin is less than 3 nucleotides long.')
            return None

        #Check if the hairpin is a special case hairpin with precalculated energy values
        elif (self._closingPair in SpecialHairpins) and (self._sequence in SpecialHairpins[self._closingPair]):
            return SpecialHairpins[self._closingPair][self._sequence]

        #Hairpins of length 3
        elif self._sequenceLen == 3:
            #get hairpin initiation term
            if self._sequenceLen in HairpinInit: #try to get from dictionary
                init = HairpinInit[self._sequenceLen]
            else: #otherwise calculate
                init = HairpinInit[9] + (1.75 * R * T * np.log(float(self._sequenceLen/9.0)))

            #check for all c loop penalty
            if self._sequence.count('C') == self._sequenceLen:
                return init + HAIRPIN_C3

            return init

        #hairpins of 4 nucleotides or greater
        else:
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
                logging.warning(f'In energy() function for Hairpin: {self._label}, terminal mismatch parameters for closing pair: {self._closingPair} and first mismatch: {firstMismatch} not found in Dictionary.')
                if strict:
                    return None #strict mode - only calculate energy for hairpins with valid params
                else:
                    terminalMismatch = 0

            #UU/GA first mismatch bonus
            uu_ga_bonus = 0
            if firstMismatch == ('U', 'U') or firstMismatch == ('G', 'A'):
                uu_ga_bonus = HAIRPIN_UU_GA_FIRST_MISMATCH_BONUS

            #GG first mismatch
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
self._span -- (int, int) -- Tuple containing the integer start and stop indices for the RNA sequene that defines the bulge.
self._closingPair5p -- (string, string) -- Tuple containing 2 single character strings. The first string the the 5' base in 5' closing pair for the bule. The second character is the 3' base in the 5' closing pair.
self._closingPair5pSpan -- (int, int) -- Tuple containing 2 integers. The first integer is the index of the 5' base in 5' closing pair for the bulge. The second integer is the index of the 3' base in the 5' closing pair
self._closingPair3p -- (string, string) -- Tuple containing 2 single character strings. The first string the the 5' base in 3' closing pair for the bule. The second character is the 3' base in the 3' closing pair.
self._closingPair3pSpan -- (int, int) -- Tuple containing 2 integers. The first integer is the index of the 5' base in 3' closing pair for the bule. The second integer is the index of the 3' base in the 3' closing pair
self._pk -- int -- the pseudoknot the bulge is a part of, if any(default value is None)
self._neighbot5p -- str -- label for the 5'neighbor of the bulge
self._neighbot3p -- str -- label for the 3'neighbor of the bulge



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
    def __init__(self, label=None, sequence='', sequenceSpan=(-1, -1), closingPair5p=('', ''), closingPair5pSpan=(-1, -1), closingPair3p=('', ''), closingPair3pSpan=(-1, -1), pk=None, neighbor5p=None, neighbor3p=None):
        self._label = label
        self._sequence = sequence
        self._sequenceLen = len(sequence)
        self._span = sequenceSpan
        self._closingPair5p = closingPair5p
        self._closingPair5pSpan = closingPair5pSpan
        self._closingPair3p = closingPair3p
        self._closingPair3pSpan = closingPair3pSpan
        self._pk = pk
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    ### Internal Methods
    ###

    #defines string representation for object
    def __str__(self):
        return f'Bulge: {self._label}'

    #define len function operation for Bulge objects
    def __len__(self):
        return self._sequenceLen

    #internal method to set the 5' and 3' neighbors for a bulge
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    ### User Accesible Functions
    ###

    '''
    Function: Bulge.label()
    Description: Function returns the label for the bulge object. Also allows user to define new label
    Parameters:
    Return Value:
    '''
    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    '''
    Function: Bulge.sequence()
    Description: Function returns the sequence that defines the bulge structure. Also allows user to define new sequence
    Parameters:
            (newSequence=None) -- str -- str nucleotides sequence to define new bulge sequence
    Return Value:
            str - the current nucleotide sequence that defines the bulge
    '''
    def sequence(self, newSequence=None):
        if newSequence:
            self._sequence = newSequence
            self._sequenceLen = len(newSequence)
        else:
            return self._sequence


    '''
    Function: Bulge.span()
    Description: Function returns the start and stop indices of the bulge as a tuple. Ex: (start, stop)
    Parameters: None
    Return Value:
            (int, int) - tuple containing the start and stop index positions of the bulge loop
    '''
    def span(self):
        return self._span


    '''
    Function: Bulge.sequenceLen()
    Description: Function returns the length of the bulge
    Parameters: None
    Return Value:
            int - the integer value length of the bulge
    '''
    def sequenceLen(self):
        return self._sequenceLen


    '''
    Function: Bulge.closingPair5p()
    Description: Function returns a tuple containg the 5' closing pair for the bulge. Also allows user to define new closing pair
    Parameters:
            (newClose=None) -- (str, str) -- new 5' closing pair for the bulge object
    Return Value:
            (str, str) - a tuple containing the 5' closing base pairs for the bulge object
    '''
    def closingPair5p(self, newClose=None):
        if newClose:
            self._closingPair5p = newClose
        else:
            return self._closingPair5p


    '''
    Function: Bulge.closingPair5pSpan()
    Description: Function returns a tuple containg the indices of the 5' closing pair for the bulge
    Parameters: None
    Return Value:
            (int, int) - tuple containing the index locations of the 5' closing base pair
    '''
    def closingPair5pSpan(self):
        return self._closingPair5pSpan


    '''
    Function: Bulge.closingPair3p(newClose=None)
    Description: Function returns a tuple containg the 3' closing pair for the bulge. Also allows user to define new closing pair
    Parameters:
            (newClose=None) -- (str, str) -- new 3' closing pair for the bulge object
    Return Value:
            (str, str) - a tuple containing the 3' closing base pairs for the bulge object
    '''
    def closingPair3p(self, newClose=None):
        if newClose:
            self._closingPair3p = newClose
        return self._closingPair3p


    '''
    Function: Bulge.closingPair5pSpan()
    Description: Function returns a tuple containg the indices of the 5' closing pair for the bulge
    Parameters: None
    Return Value:
            (int, int) - tuple containing the index locations of the 5' closing base pair
    '''
    def closingPair3pSpan(self):
        return self._closingPair3pSpan


    '''
    Function: Bulge.neighbors()
    Description: function to get the StructureComponents directly adjacent to the bulge
    Parameters: None
    Return Value:
            (str, str) - tuple containg the labels for the 5' and 3' neighbors of the bulge
    '''
    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)


    '''
    Function: Bulge.canonical()
    Description: function to check for valid conditions for calculating bulge energy
    Parameters: None
    Return Value:
            bool - returns True if all valid energy parameters are present for energy calculation
    '''
    def canonical(self):
        if self._sequenceLen == 1:
            if (self._closingPair5p not in StackingEnergies) or (self._closingPair3p not in StackingEnergies[self._closingPair5p]):
                return False
        return True


    '''
    Function: Bulge.energy()
    Description: function calculates the folding free energy change for the bulge
    Parameters:
            (strict=True) -- bool -- when true only energy values for bulges with all valid energy parameters will be calaculated
    Return Value:
            float - the calculated energy of the Bulge
    '''
    def energy(self, strict=True):
        if self._sequenceLen == 1: #bulges of length 1
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

            #check for special C bulge case
            #special C condition = sequence is all 'C' with at least one adjacent 'C'
            specialC = 0 #specialC stores the SPECIAL_C_BULGE value being applied to the particular bulge
            cCount = 0 #cCount stores the number of adjacent C's for the number of states component of the equation
            if self._sequence == 'C' and (self._closingPair5p[0] == 'C' or self._closingPair3p[0] == 'C'):
                specialC = SPECIAL_C_BULGE
                cCount = 1 #number of possible states due to adjacent C's
                if(self._closingPair5p[0] == 'C'):
                    cCount += 1
                if (self._closingPair3p[0] == 'C'):
                    cCount += 1

                return BulgeInit[1] + basePairStack + specialC - (R * T * np.log(cCount))

            #if not special C bulge, return bulge init + basePairStack
            else:
                return BulgeInit[1] + basePairStack

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



              CGC
            AG   CU
            ||   ||
            UC   GA
               A
                 ^3' closing pair
             ^5' closing pair

'''
class InternalLoop:
    # __init__ method for InternalLoop object
    def __init__(self, pLabel=None, label5p=None, label3p=None,  loop5p='', loop3p='', loop5pSpan=(-1, -1), loop3pSpan=(-1, -1), closingPairs=(('', ''), ('', '')), closingPairsSpan=((-1, -1), (-1, -1)), neighbor5p=('', ''), neighbor3p=('', '')):
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
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor5p
        self._strict = True #used for to control energy function

    ###
    ### Internal Methods
    ###

    #defines the string representation of the object
    def __str__(self):
        return f'Inner Loop: {self._parentLabel}'

    #define len function operation for InnerLoop objects
    def __len__(self):
        return self._loopsLen

    #function to update loop lengths upon change
    def _updateLoopLen(self):
        self._loopsLen = (len(self._5pLoop), len(self._3pLoop))

    #internal method to set the 5' and 3' neighbors for a InternalLoop
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    #Function checks that the internal loop has the same 5' closing pair structures
    def _same5pNeighbors(self):
        return (self._neighbor5p[0] == self._neighbor5p[1])

    #Function checks that the internal loop has the same 3' closing pair structures
    def _same3pNeighbors(self):
        return (self._neighbor3p[0] == self._neighbor3p[1])

    ###
    ### User Accessible Methods
    ###


    '''
    Function: InternalLoop.label()
    Description: Function returns the parent label for the inner loop. Also allows user to set new label
    Parameters:
            (newLabel=None) -- str -- new label to identify the InternalLoop object
    Return Value:
            str - the current label for the InternaLoop object
    '''
    def label(self, newLabel=None):
        if newLabel:
            self._parentLabel = newLabel
        else:
            return self._parentLabel


    '''
    Function: InternalLoop.loops()
    Description: Function returns a tuple containing the sequences for the two inner loop subcomponents. Also allows user to define new loop sequences.
    Parameters:
            (loop5p=None) -- str -- new sequence to define the 5' portion of the internal loop
            (loop3p=None) -- str -- new sequence to define the 3' portion of the internal loop
    Return Value:
            (str, str) - tuple containing the 5' and 3' portions of the internal loop sequence
    '''
    def loops(self, loop5p=None, loop3p=None):
        if(loop5p and loop3p):
            self._5pLoop = loop5p
            self._3pLoop = loop3p
            self._updateLoopLen()
        else:
            return (self._5pLoop, self._3pLoop)


    '''
    Function: InternalLoop.loop5p()
    Description: Function that returns the 5' portion of the inner loop. Also allows user to define 5' portion of the loop
    Parameters:
            (loop=None) -- str -- new nucleotide sequence to define the 5' portion of the Internal Loop
    Return Value:
            str - the current sequence that defines the 5' portion of the InternalLoop
    '''
    def loop5p(self, loop=None):
        if(loop):
            self._5pLoop = loop
            self._updateLoopLen()
        else:
            return self._5pLoop


    '''
    Function: InternalLoop.loop3p()
    Description: Function that returns the 3' portion of the inner loop. Also allows user to define 3' portion of the loop
    Parameters:
            (loop=None) -- str -- new nucleotide sequence to define the 3' portion of the Internal Loop
    Return Value:
            str - the current sequence that defines the 3' portion of the InternalLoop
    '''
    def loop3p(self, loop=None):
        if(loop):
            self._3pLoop = loop
            self._updateLoopLen()
        else:
            return self._3pLoop


    '''
    Function: InternalLoop.loopsLen()
    Description: Function returns a tuple containing the the integer value lengths of the two inner loop components
    Parameters: None
    Return Value:
            (int, int) - tuple containing the integer value sequence lengths of the 5' and 3' portions of the Internal Loop
    '''
    def loopsLen(self):
        return self._loopsLen


    '''
    Function: InternalLoop.span()
    Description: Function returns a tuple that contains two tuples containing the integer start and stop positions of the 5' and 3' inner loop components
    Parameters: None
    Return Value:
            ((int, int), (int, int)) - tuple containing 2 tuples that each define the start and stop indices for the 5' and 3' portions of the Internal Loop
    '''
    def span(self):
        return (self._span5p, self._span3p)


    '''
    Function: InternalLoop.closingPairs()
    Description: Function returns a tuple that contains two tuples containing the closing base pairs of the inner loop components
    Parameters: None
    Return Value:
            ((str, str), (str, str)) - tuple containg 3 tuples that define the closing nucleotide base pairs for the 5' and 3' end of the Internal Loop
    '''
    def closingPairs(self):
        return self._closingPairs


    '''
    Function: InternalLoop.closingPairsSpan()
    Description: Function returns a tuple that contains two tuples containing the index locations of the closing base pairs of the inner loop components
    Parameters: None
    Return Value:
            ((int, int), (int, int)) - uple containg 3 tuples that define the index locations of the closing nucleotide base pairs for the 5' and 3' end of the Internal Loop
    '''
    def closingPairsSpan(self):
        return self._closingPairsSpan


    '''
    Function: InternalLoop.neighbors()
    Description: function to get the StructureComponents directly adjacent to the InternalLoop
    Parameters: None
    Return Value:
            ((str, str), (str, str)) - tuple containg 2 tuples that define the neighboring structures to the 5' and 3' ends of the Internal Loop
    '''
    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)


    '''
    Function: InternalLoop.canonical()
    Description: Function to check if valid parameters are available to calculate inner loop energy
    Parameters: None
    Return Value:
            bool - returns True if there is a complete set of parameters for calculating the energy of the internal loop
    '''
    def canonical(self):
        #Check if energy value is present for 1x1 loop
        if len(self._5pLoop) == 1 and len(self._3pLoop) == 1:
            if (self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop) in InnerLoop_1x1_Energies:
                return True
            return False
        #check if energy value is present for 1x2 loop
        elif len(self._5pLoop) == 1 and len(self._3pLoop) == 2:
            if (self._closingPairs[0], self._closingPairs[1], self._5pLoop, self._3pLoop[1], self._3pLoop[0]) in InnerLoop_1x2_Energies:
                return True
            return False
        #check if energy value is present for 2x1 loop
        elif len(self._5pLoop) == 2 and len(self._3pLoop) == 1:
            if ((self._closingPairs[1][1], self._closingPairs[1][0]), (self._closingPairs[0][1], self._closingPairs[0][0]), self._3pLoop, self._5pLoop[1], self._5pLoop[0]) in InnerLoop_1x2_Energies:
                return True
            return False
        #Check if energy value is present for 2x2 loop
        elif len(self._5pLoop) == 2 and len(self._3pLoop) == 2:
            loops = list(zip(list(self._5pLoop), list(self._3pLoop[::-1]))) #convert loop sequences to proper format for dictionary
            if (self._closingPairs[0], self._closingPairs[1], loops[0], loops[1]) in InnerLoop_2x2_Energies:
                return True
            return False
        #Check for valid parameters needed to calculate energy for loops of other lengths
        else:
            if(self._getInnerLoopMismtachEnergy() != None):
                return True
            return False


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
            mismatchEnergy_3x2 += InnerLoopMismatches_2x3[((self._closingPairs[1][1], self._closingPairs[1][0]), mismatch5p)]
        else:
            logging.warning(f'In energy() function for 3x2 InnerLoop: {self._parentLabel}, no mismatch parameter for closing pair: {(self._closingPairs[1][1], self._closingPairs[1][0])} and the 5\' mismatch: {mismatch5p}.')
            if (self._strict):
                return None

        #check for mismatch condition between 3'closing pair and mismatch 2
        if ((self._closingPairs[0][1], self._closingPairs[0][0]), mismatch3p) in InnerLoopMismatches_2x3:
            mismatchEnergy_3x2 += InnerLoopMismatches_2x3[((self._closingPairs[0][1], self._closingPairs[0][0]), mismatch3p)]
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
            mismatchEnergy_2x3 += InnerLoopMismatches_2x3[((self._closingPairs[1][1], self._closingPairs[1][0]), mismatch3p)]
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
self._sequenceLen -- int -- length of the external loop sequence
self._span --tuple(int, int) -- tuple containing the integer start and stop locations for the external loop sequence
self._closingPair5p -- tuple(string, string) -- tuple that contains the 5' closing base pair for the external loop
self._closingPair5pSpan -- tuple(int, int) -- tuple containg the integer index locations for the 5' closing pair
self._closingPair3p -- tuple(string, string) -- tuple that contains the 3' closing base pair for the external loop
self._closingPair3pSpan -- tuple(int, int) -- tuple containg the integer index locations for the 3' closing pair
'''
class ExternalLoop:
    #__init__() method for the external loop object
    def __init__(self, label='', seq='', seqSpan=(-1,-1), closingPair5p=('', ''), closingPair5pSpan=(-1,-1), closingPair3p=('', ''), closingPair3pSpan=(-1, -1), neighbor5p=None, neighbor3p=None):
        self._label = label
        self._sequence = seq
        self._sequenceLen = len(seq)
        self._span = seqSpan
        self._closingPair5p = closingPair5p
        self._closingPair5pSpan = closingPair5pSpan
        self._closingPair3p = closingPair3p
        self._closingPair3pSpan = closingPair3pSpan
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    ### Internal Methods
    ###

    #Defines the string representation of the external loop
    def __str__(self):
        return f'External Loop: {self._label}'

    #define len function operation for ExternalLoop objects
    def __len__(self):
        return self._sequenceLen

    #internal method to set the 5' and 3' neighbors for a external loop
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    ### User Accesible Methods
    ###


    '''
    Function: ExternalLoop.label()
    Description: Function returns the label for the external loop
    Parameters:
            (newLabel=None) -- str -- new label to define the ExternalLoop Object
    Return Value:
            str - the current label for the ExternalLoop object
    '''
    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label


    '''
    Function: ExternalLoop.sequence()
    Description: Function returns the sequence that defines the external loop
    Parameters:
            (newSequence=None) -- str -- new nucleotide sequence to define the ExternalLoop
    Return Value:
            str - current sequence that defines the ExternalLoop
    '''
    def sequence(self):
        return self._sequence


    '''
    Function: ExternalLoop.sequenceLen()
    Description: function returns the length of the external loop sequence
    Parameters: None
    Return Value:
            int - the integer value length of the ExternalLoop sequence
    '''
    def sequenceLen(self):
        return self._sequenceLen


    '''
    Function: ExternalLoop.span()
    Description: Function returns a tuple containing the start and stop index locations for the external loop sequence
    Parameters: None
    Return Value:
            (int, int) - tuple containing the start and stop index locations of the ExternalLoop
    '''
    def span(self):
        return self._span


    '''
    Function: ExternalLoop.neighbors()
    Description: function to get the StructureComponents directly adjacent to the external loop
    Parameters: None
    Return Value:
            (str, str) - tuple containing the labels for the 5' and 3' neighbors of the ExternalLoop
    '''
    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)


'''
ENDS

Member variable -- data type -- description:
self._label -- string -- label for the end objects
self._sequence -- string -- sequence that defines the end objects
self._span -- tuple(int, int) -- tuple containing the integer start and stop locations for the end object
'''
class End:
    #__init__() method for end object
    def __init__(self, label='', sequence='', span=(-1, -1), neighbor=None):
        self._label = label
        self._sequence = sequence
        self._sequenceLen = len(sequence)
        self._span = span
        self._neighbor=None

    ###
    ### Internal Methods
    ###

    #define string representation of end object
    def __str__(self):
        return f'End: {self._label}'

    #define len function operation for End objects
    def __len__(self):
        return self._sequenceLen

    #internal method to set the 5' and 3' neighbors for a end
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    ### User Accesible Methods
    ###


    '''
    Function: End.label()
    Description: Function returns the label for the end object
    Parameters: 
            (newLabel=None) -- str -- new label to define the End object
    Return Value:
            str - the current label for the End object
    '''
    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    
    '''
    Function: End.sequence()
    Description: Function returns the sequence that defines the end object
    Parameters:
            (newSequence=None) -- str -- new nucleotide sequence to define the End object
    Return Value:
            str - the current sequence that defines the End object
    '''
    def sequence(self, newSequence=None):
        if newSequence:
            self._sequence = newSequence
        else:
            return self._sequence

    
    '''
    Function: End.sequenceLen()
    Description: function returns the length of the End sequence
    Parameters: None
    Return Value:
            int - the integer value length of the End object
    '''
    def sequenceLen(self):
        return self._sequenceLen

    
    '''
    Function: End.span()
    Description: Function returns a tuple that contains the integer start and stop index locations for the end object
    Parameters: None
    Return Value:
            (int, int) - tuple containing the integer value start and stop indices of the End object
    '''
    def span(self):
        return self._span

    
    '''
    Function: End.neighbors()
    Description: function to get the StructureComponents directly adjacent to the end
    Parameters: None
    Return Value:
            (str, str) - tuple containing the labels for the StructureComponents neighboring the End in the Structure object
    '''
    def neighbors(self):
        return (self._neighbor5p, self._neighbor3p)


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

    ###
    ### Internal Methods
    ###

    #Defines the string representation of the NCBP object
    def __str__(self):
        return f'NCBP: {self._label}'

    ###
    ### User Accesible Methods
    ###

    
    '''
    Function: NCBP.label()
    Description: Functions returns the label for the NCBP object
    Parameters: 
            (newLabel=None) -- str -- new label to identify the NCBP object
    Return Value:
            str - the current label for the NCBP object
    '''
    def label(self, newLabel=None):
        if newLabel:
            self._label = newLabel
        else:
            return self._label

    
    '''
    Function: NCBP.pair()
    Description: Function returns a tuple containing the two base pairs that define the NCBP object
    Parameters:
            (newPair=None) -- (str, str) -- new tuple to define the NCBP
    Return Value:
            (str, str) - the current base pair that defines the NCBP
    '''
    def pair(self, newPair=None):
        if newPair:
            self._basePair = newPair
        else:
            return self._basePair

    
    '''
    Function: NCBP.span()
    Description: Function returns a tuple containing the integer locations of the base pairs that define the NCBP
    Parameters: None
    Return Value:
            (int, int) - tuple containing the index locations of the NCBP
    '''
    def span(self):
        return self._basePairSpan

    
    '''
    Function: NCBP.parentUnit()
    Description: Function returns a string the identifies the secondary structure that the NCBP occurs in
    Parameters: None
    Return Value:
            str - label for the StructureComponent that contains the NCBP
    '''
    def parentUnit(self):
        return self._parentUnit



'''
MultiLoop

Member variable -- data type -- description:
self._parentLabel -- str -- parent label for all the multiloop subcomponents
self._subunitLabels -- list -- list of all the subunit labels for the multiloop
self._numSubunits -- int -- number of subunits composing the multiloop
self._sequences -- dictionary -- dictionary of the multiloop component sequences. key values are the subunit labels
self._span -- dictionary -- dictionary of the multiloop component spans. key values are the subunit labels
self._closingPairs -- ((str, str), (str, str)) -- tuple containing the 5' and 3' closing base pairs as tuples
self._closingPairsSpan - ((int, int), (int, int)) -- tuple containing the 5' and 3' closing base pair spans as tuples
'''
class MultiLoop:
    #__init__() method for MultiLoop class
    def __init__(self, parentLabel, subunitLabels, sequences, spans, closingPairs, closingPairsSpan):
        self._parentLabel = parentLabel
        self._subunitLabels = subunitLabels
        self._numSubunits = len(sequences)
        self._sequences = sequences
        self._spans = spans
        self._closingPairs = closingPairs
        self._closingPairsSpan = closingPairsSpan

    ###
    ### Internal Methods
    ###

    #define string representation for MultiLoop object
    def __str__(self):
        return f'MultiLoop: {self._parentLabel}'

    #internal method to set the 5' and 3' neighbors for a multiloop
    def _addNeighbors(self, neighbor5p, neighbor3p):
        self._neighbor5p = neighbor5p
        self._neighbor3p = neighbor3p

    ###
    ### User Accesible Methods
    ###

    
    '''
    Function: MultiLoop.label()
    Description: Function to return the parent Label for the MultiLoop object
    Parameters:
            (newLabel=None) -- str -- new label to identify the Multiloop object
    Return Value:
            str - the label for the current MultlLoop object
    '''
    def label(self, newLabel=None):
        if newLabel:
            self._parentLabel = newLabel
        else:
            return self._parentLabel

    
    '''
    Function: MultiLoop.subunitLabels()
    Description: Function to return a list of the MultiLoop object subcomponents
    Parameters: None
    Return Value:
            list - list of the subunit labels for the MultiLoop object
    '''
    def subunitLabels(self):
        return self._subunitLabels

    
    '''
    Function: MultiLoop.numSubunits()
    Description: Function to return the number of subcomponents composing the MultiLoop
    Parameters: None
    Return Value:
            int - the number of subunits that compose the MultiLoop Structure
    '''
    def numSubunits(self):
        return self._numSubunits

    
    '''
    Function: MultiLoop.sequence()
    Description: Function to return dictionary of subunitLabel : Sequence pairs for the MultiLoop object
    Parameters:
            (subunit=None) -- str -- label for the specific subunit being accessed.
    Return Value:
            dict - dictionary of multiloop sequences mapped to their subunit label
            * if a specific subunit label is provided, only that sequence will bw returned
    '''
    def sequence(self, subunit=None):
        if(subunit):
            try:
                sequence = self._sequences[subunit]
                return sequence
            except KeyError:
                return None
        else:
            return self._sequences

    
    '''
    Function: MultiLoop.span()
    Description: Function to return dictionary of subunitLabel : SequenceSpans pairs for the MultiLoop object
    Parameters:
            (subunit=None) -- str -- label for specific subunit being accessed
    Return Value:
            dict - dictionary of tuples defining start and stop indicies for multiloop subunits mapped to their subunit label
            * if a specific subunit label is provided, only that tuple will be returned
    '''
    def span(self, subunit=None):
        if(subunit):
            try:
                span = self._spans[subunit]
                return span
            except KeyError:
                return None
        else:
            return self._spans

    
    '''
    Function: MultiLoop.closingPairs()
    Description: Function to return dictionary of subunitLabel : closingPairs for the MultiLoop object
    Parameters:
            (subunit=None) -- str -- label for specific subunit being accessed
    Return Value:
            dict - dictionary of tuples defining closing base pairs for multiloop subunits mapped to their subunit label
            * if a specific subunit label is provided, only that tuple will be returned
    '''
    def closingPairs(self, subunit=None):
        if(subunit):
            try:
                closingPair = self._closingPairs[subunit]
                return closingPair
            except KeyError:
                return None
        else:
            return self._closingPairs

    
    '''
    Function: MultiLoop.closingPairsSpan()
    Description: Function to return dictionary of subunitLabel : closingPairsSpan for the MultiLoop object
    Parameters:
            (subunit=None) -- str -- label for specific subunit being accessed
    Return Value:
            dict - dictionary of tuples defining closing base pair index locations for multiloop subunits mapped to their subunit label
            * if a specific subunit label is provided, only that tuple will be returned
    '''
    def closingPairsSpan(self, subunit=None):
        if(subunit):
            try:
                closingPairsSpan = self._closingPairsSpan[subunit]
                return closingPairsSpan
            except KeyError:
                return None
        else:
            return self._closingPairsSpan



'''
PSEUDOKNOTS -- UNFINISHED
'''
class PseudoKnot:
    pass
