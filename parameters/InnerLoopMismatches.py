'''
Filename: InnerLoopMismatches.py
Author: Michael Hathaway

Description: Module imported by the StructureType module that is used for inner loop energy calculations.
Data Source: https://rna.urmc.rochester.edu/NNDB/turner04/internal-parameters.html
'''

InnerLoopMismatches_2x3 = { #R = purine, Y = pyrimidine
    #(('R', Y), ('A', 'G'))
    (('A', 'U'), ('A', 'G')) : 0,
    (('G', 'U'), ('A', 'G')) : 0,
    (('A', 'C'), ('A', 'G')) : 0,
    (('G', 'C'), ('A', 'G')) : 0,
    #(('Y', 'R'), ('A', 'G'))
    (('U', 'A'), ('A', 'G')) : -0.5,
    (('C', 'A'), ('A', 'G')) : -0.5,
    (('U', 'G'), ('A', 'G')) : -0.5,
    (('C', 'G'), ('A', 'G')) : -0.5,
    #(('R', 'Y'), ('A', 'G'))
    (('A', 'U'), ('G', 'A')) : -1.2,
    (('G', 'U'), ('G', 'A')) : -1.2,
    (('A', 'C'), ('G', 'A')) : -1.2,
    (('G', 'C'), ('G', 'A')) : -1.2,
    #(('Y', 'R'), ('G', 'A'))
    (('U', 'A'), ('G', 'A')) : -1.1,
    (('C', 'A'), ('G', 'A')) : -1.1,
    (('U', 'G'), ('G', 'A')) : -1.1,
    (('C', 'G'), ('G', 'A')) : -1.1,
    #(('R/Y', 'Y/R'), ('G', 'G'))
    (('A', 'U'), ('G', 'G')) : -0.8,
    (('G', 'U'), ('G', 'G')) : -0.8,
    (('A', 'C'), ('G', 'G')) : -0.8,
    (('G', 'C'), ('G', 'G')) : -0.8,
    (('U', 'A'), ('G', 'G')) : -0.8,
    (('C', 'A'), ('G', 'G')) : -0.8,
    (('U', 'G'), ('G', 'G')) : -0.8,
    (('C', 'G'), ('G', 'G')) : -0.8,
    #(('R/Y', 'Y/R'), ('U', 'U'))
    (('A', 'U'), ('U', 'U')) : -0.4,
    (('G', 'U'), ('U', 'U')) : -0.4,
    (('A', 'C'), ('U', 'U')) : -0.4,
    (('G', 'C'), ('U', 'U')) : -0.4,
    (('U', 'A'), ('U', 'U')) : -0.4,
    (('C', 'A'), ('U', 'U')) : -0.4,
    (('U', 'G'), ('U', 'U')) : -0.4,
    (('C', 'G'), ('U', 'U')) : -0.4
}

OtherInnerLoopMismtaches = {
    ('A', 'G') : -0.8,
    ('G', 'A') : -1.0,
    ('G', 'G') : -1.2,
    ('U', 'U') : -0.7,
}
