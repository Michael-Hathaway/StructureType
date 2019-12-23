'''
Filename: LoopInitiationEnergy.py
Author: Michael Hathaway

Description:
Module imported by the StructureType module to be used for energy calculations for bulge, internal loop, and hairpin secondary structures. The content of this file was generated using the parseLoopInitiationEnergies.py script found in the TurnerParameters/scripts directory. The script parses the loopInitiationEnergy.txt file found in the TurnerParameters/parameterTextFiles directory to produce the dictionary below.

Original Source: https://rna.urmc.rochester.edu/NNDB/turner04/loop.txt

Dictionary Usage:
The keys in the dictionary are the lengths of the secondary structure sequences. The values are the initiation energies measured in Kcal/mol
'''

InternalLoopInit = {1: None, 2: None, 3: None, 4: 1.1, 5: 2.0, 6: 2.0, 7: 2.1, 8: 2.3, 9: 2.4, 10: 2.5, 11: 2.6, 12: 2.7, 13: 2.8, 14: 2.9, 15: 2.9, 16: 3.0, 17: 3.1, 18: 3.1, 19: 3.2, 20: 3.3, 21: 3.3, 22: 3.4, 23: 3.4, 24: 3.5, 25: 3.5, 26: 3.5, 27: 3.6, 28: 3.6, 29: 3.7, 30: 3.7}
BulgeInit = {1: 3.8, 2: 2.8, 3: 3.2, 4: 3.6, 5: 4.0, 6: 4.4, 7: 4.6, 8: 4.7, 9: 4.8, 10: 4.9, 11: 5.0, 12: 5.1, 13: 5.2, 14: 5.3, 15: 5.4, 16: 5.4, 17: 5.5, 18: 5.5, 19: 5.6, 20: 5.7, 21: 5.7, 22: 5.8, 23: 5.8, 24: 5.8, 25: 5.9, 26: 5.9, 27: 6.0, 28: 6.0, 29: 6.0, 30: 6.1}
HairpinInit = {1: None, 2: None, 3: 5.4, 4: 5.6, 5: 5.7, 6: 5.4, 7: 6.0, 8: 5.5, 9: 6.4, 10: 6.5, 11: 6.6, 12: 6.7, 13: 6.8, 14: 6.9, 15: 6.9, 16: 7.0, 17: 7.1, 18: 7.1, 19: 7.2, 20: 7.2, 21: 7.3, 22: 7.3, 23: 7.4, 24: 7.4, 25: 7.5, 26: 7.5, 27: 7.5, 28: 7.6, 29: 7.6, 30: 7.7}
