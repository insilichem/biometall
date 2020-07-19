"""
box.py
Module with auxiliar functions
"""

from collections import Counter

def _counterSubset(list1, list2):
    """
    Check if all the elements of list1 are contained in list2.
    
    It counts the quantity of each element (i.e. 3-letter amino acid code) 
    in the list1 (e.g. two 'HIS' and one 'ASP') and checks if the list2 contains 
    at least that quantity of every element.

    Parameters
    ----------
    list1 : list
        List of amino acids in 3-letter code (can be repetitions)
    list2 : list
        List of amino acids in 3 letter code (can be repetitions)

    Returns
    -------
    bool
        True if list2 contains all the elements of list1. False otherwise
    """
    #1. Count how many amino acids of each kind
    c1, c2 = Counter(list1), Counter(list2)
    
    #2. Check that list2 contains at least the same number of every amino acid
    for k, n in c1.items():
        if n > c2[k]:
            return False
    return True