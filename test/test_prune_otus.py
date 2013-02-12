'''
Created on Feb 12, 2013

@author: shareef
'''
import unittest


class PruneOTUsTest(unittest.TestCase):

    def setUp(self):
        self.otus = {'111111': ['SID1_10', 'SID1_22']}


    def tearDown(self):
        pass


    def test_filter_by_sample_pct(self):
        pass
    
    def test_filter_by_sequence_pct(self):
        pass


if __name__ == "__main__":
    unittest.main()