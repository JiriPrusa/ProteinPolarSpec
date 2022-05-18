from __future__ import absolute_import
from __future__ import print_function
__author__ = "Jiri Prusa"
__version__ = "0.1"


class ChainSelector(object):
    """
    ChainSelector description
    """

    def __init__(self, topology):
        self._chains = topology.chains()
        self._atoms = topology.atoms()
        
    
    def getIndexListForChain(self, chainLett):
        # check for atom indexes that coresponds to chain 
        # we are interested in
        reportIndexes = []
        chainIdx = None
        for (chainIndex, chain) in enumerate(self._chains):
                chainName = chr(ord('A')+chainIndex%26)
                if chainName == chainLett:
                    chainIdx = chain
        
        if chainIdx is not None:
            for atom in self._atoms:
                if chainIdx == atom.residue.chain:
                    reportIndexes.append(atom.index)
        else:
            print("Warning!!! Chain %s not found in topology!" % chainLett)
            return 

       #print("Fond chain %s with %i elements." % (self.chainIdx, len(reportIndexes)))
        return reportIndexes
    
       
    def printChains(self):
        for (chainIndex, chain) in enumerate(self._chains):
            chainName = chr(ord('A')+chainIndex%26)
            print("Fond chain %s (%i)." % (chainName, chain))
        return
