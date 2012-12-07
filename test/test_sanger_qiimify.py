'''
Created on Dec 3, 2012

@author: shareef
'''
import unittest

import sanger_qiimify as sq

class Test(unittest.TestCase):


    def test_generate_barcodes(self):
        num_codes = 50
        bc = sq.generate_barcodes(num_codes)
        self.assertEqual(len(bc), num_codes)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()