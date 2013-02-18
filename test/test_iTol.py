'''
Created on Feb 16, 2013

@author: Shareef Dabdoub
''' 
import iTol
import unittest

class iTolTest(unittest.TestCase):
    pass
        

class GatherCategoriesTest(iTolTest):

    def setUp(self):
        self.header = ('#SampleID BarcodeSequence LinkerPrimerSequence ' + 
                      'DiseaseState SmokingStatus SampleLocation Description')
        self.header = self.header[1:].split()
        self.categories = ['DiseaseState', 'SmokingStatus', 
                           'SampleLocation=implant']
        self.imap = {'SID1':['SID1','AA','GG','healthy','smoker','implant',''],
                     'SID2':['SID2','AA','GG','pm','smoker','implant',''],
                     'SID3':['SID3','AA','GG','pi','smoker','implant',''],
                     'SID4':['SID4','AA','GG','pm','ns','tooth',''],
                     'SID5':['SID5','AA','GG','healthy','ns','implant',''],
                     'SID6':['SID6','AA','GG','pi','smoker','tooth',''],
                     'SID7':['SID7','AA','GG','pm','ns','implant',''],
                     'SID8':['SID8','AA','GG','pi','ns','implant','']}

    def tearDown(self):
        pass
        
    def test_none_categories(self):
        cats = iTol.gather_categories(self.imap, self.header, None)
        self.assertEqual(len(cats), 1, msg='Result should have only 1 entry')
        self.assertIn('default', cats, msg='key "default" not found in result')
        
    def test_no_categories(self):
        cats = iTol.gather_categories(self.imap, self.header, ['Treatment'])
        self.assertEqual(len(cats), 1, msg='Result should have only 1 entry')
        self.assertIn('default', cats, msg='key "default" not found in result')
        
    def test_expected_input(self):
        cats = iTol.gather_categories(self.imap, self.header, self.categories)
        expected = {'healthy_smoker':iTol.DataCategory({'SID1'},{},{}),
                    'pm_smoker':iTol.DataCategory({'SID2'},{},{}),
                    'pi_smoker':iTol.DataCategory({'SID3'},{},{}),
                    'healthy_ns':iTol.DataCategory({'SID5'},{},{}),
                    'pm_ns':iTol.DataCategory({'SID7'},{},{}),
                    'pi_ns':iTol.DataCategory({'SID8'},{},{})}
        self.assertDictContainsSubset(expected, cats)
        
if __name__ == "__main__":
    unittest.main()