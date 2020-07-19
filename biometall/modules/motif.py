"""
motif.py
Module with functions to hadle combinations of motifs/mutations
"""

import itertools
from .box import _counterSubset

def _check_actual_motif(motif_possibilities, probe_coordinators, name_for_res, 
                            res_for_column):
    """
    Checks if the coordinators of a probe match the motif requested by the user.

    A motif is matched when the possible amino acids that coordinate a probe
    include all the contained in the motif. Only the actual amino acids present 
    in the structure provided by the user are considered. For example, a motif
    [HIS,HIS,ASP/GLU] would be accomplished if the coordinators of the probe 
    contains, at least, either (HIS,HIS,ASP) or (HIS,HIS,GLU).

    The function returns all the possible combinations of amino acids that match
    the motif for the given probe.

    Parameters
    ----------
    motif_possibilities : list of sequences
        All the combinations of amino acids that accomplish the motif
    probe_coordinators : dict
        Contains the amino acids that coordinate the probe. Indexed by number of
        column, the values are the names of the coordinating amino acids.
    name_for_res : dict
        Correspondence between number:chain of residue and residue name
    res_for_column: dict
        Correspondence between column number and number:chain of residue

    Returns
    -------
    list of sequences
        sequences of the amino acids that match the motif for the queried probe
    """
    c = []
    for i in list(probe_coordinators):
        if name_for_res[res_for_column[i]] in probe_coordinators[i]:
            c.append(tuple([name_for_res[res_for_column[i]], i]))
    solutions = set()
    for m in list(itertools.combinations(c, len(motif_possibilities[0]))):
        for possible_motif in motif_possibilities:
            if(_counterSubset(possible_motif, [m_tuple[0] for m_tuple in m])):
                solutions.add(tuple(m_tuple[1] for m_tuple in m))
    return solutions

def _check_possible_mutations(mutated_motif_possibilities, 
                                residues_of_mutated_motif, probe_coordinators):
    """
    Calculates the possible mutations to achieve a motif requested by the user.

    A motif is matched when the possible amino acids that coordinate a probe
    include all the contained in the motif. Both the real and possible mutations
    of the amino acids are considered. For example, a motif [HIS,ASP/GLU] would 
    be accomplished if the coordinators of the probe contain, either in the 
    actual structure or by making a mutation, (HIS,ASP) or (HIS,GLU). 

    The function returns all the mutations that could be made in the residues
    to achieve the coordinating motif for the given probe.

    Parameters
    ----------
    mutated_motif_possibilities : list of sequences
        All the combinations of amino acids that accomplish the motif
    residues_of_mutated_motif : sequence of str
        Sequence of the names of the amino acids involved in the motif
    probe_coordinators : dict
        Contains the amino acids that coordinate the probe. Indexed by number of
        column, the values are the names of the coordinating amino acids and 
        their possible mutations.
    name_for_res : dict
        Correspondence between number:chain of residue and residue name
    res_for_column: dict
        Correspondence between column number and number:chain of residue

    Returns
    -------
    dict
        Possible mutations for a given residue (indexed by column number)
        
    """
    c = []
    for i in list(probe_coordinators):
        c.append([])
        for j in probe_coordinators[i]:
            c[-1].append(tuple([j, i]))
    
    c_possibilities = list(itertools.product(*c))
          
    mutations = {}
    for m in c_possibilities:
        for possible_motif in mutated_motif_possibilities:
            if(_counterSubset(possible_motif, [m_tuple[0] for m_tuple in m])):
                for p in probe_coordinators:
                    for r_name in probe_coordinators[p]:
                        #Exists possibility of mutation
                        if residues_of_mutated_motif.intersection(set([r_name])): 
                            if p not in list(mutations):
                                mutations[p] = set([r_name])
                            else:
                                mutations[p].add(r_name)
                break
        if mutations:
            break
    return mutations