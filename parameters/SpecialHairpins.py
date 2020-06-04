'''
filename: SpecialHairpins.py
Author: Michael Hathaway

Description: Module imported by the StructureType Module for use in hairpin energy calculations. Hairpin loops of length 3, 4, and 6 nucleotides have stabilities that are poorly fit by models so their stabilities are assigned based on experimental data.

Original source: https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-special-parameters.html

Dictionary Usage: This is a nested dictionary where the first key is the closing pair for the hairpin loop and the second key is the sequence that defines the hairpin. The values are free energy values measured in Kcal/mol.
'''

SpecialHairpins = {
    ('C', 'G') : {
        ('AAC') : 6.8,
        ('UUA') : 6.9,
        ('UACG') : 2.8,
        ('UCCG') : 2.7,
        ('UUCG') : 3.7,
        ('UUUG') : 3.7,
        ('CAAG') : 3.3,
        ('CCAG') : 3.4,
        ('CGAG') : 3.5,
        ('CUAG') : 3.7,
        ('CACG') : 3.7,
        ('CGCG') : 3.6,
        ('CUCG') : 2.5,
        ('UAAG') : 3.6,
        ('UCAG') : 3.7,
        ('UUAG') : 3.5,
        ('UGCG') : 2.8,
        ('AACG') : 5.5,
    },
    ('A', 'U') : {
        ('CAGUGC') : 2.9,
        ('CAGUGA') : 3.6,
        ('CAGUGU') : 1.8,
        ('CAGUAC') : 2.8,
    }
}
