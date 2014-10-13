import unittest
from qiime_tools import util as qtu

class SplitPhylogenyTest(unittest.TestCase):
    def test_split_phylogeny(self):
        """
        Testing split_phylogeny() function in util.py.
        """
        
        #test case
        p1 = 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Veillonella; s__denticariosi'
        
        #the following example is provided for testing different case.
        #p1 = 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Peptostreptococcaceae_[XI]; g__Eubacterium_[XI][G-6]; s__minutum'
        #p1 = 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales;'f__Veillonellaceae; g__; s__'
        
        for lvl in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
        	if lvl == 'k':
        		self.assertEquals(qtu.split_phylogeny(p1, 'k'),'k__Bacteria', msg='Error')
        	
        	if lvl == 'p':
        		self.assertEquals(qtu.split_phylogeny(p1, 'p'),'k__Bacteria; p__Firmicutes', msg='Error')
        	
        	if lvl == 'c':
        		self.assertEquals(qtu.split_phylogeny(p1, 'c'), 'k__Bacteria; p__Firmicutes; c__Clostridia', msg='Error')
        	
        	if lvl == 'o':
        		self.assertEquals(qtu.split_phylogeny(p1, 'o'), 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales', msg='Error')
        	
        	if lvl == 'f':
        		self.assertEquals(qtu.split_phylogeny(p1, 'f'), 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae', msg='Error')
        		
        	if lvl == 'g':
        		self.assertEquals(qtu.split_phylogeny(p1, 'g'), 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Veillonella', msg='Error')
        	
        	if lvl == 's':
        		self.assertEquals(qtu.split_phylogeny(p1, 's'), 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Veillonella; s__denticariosi', msg='Error')
        
if __name__ == '__main__':
    unittest.main()