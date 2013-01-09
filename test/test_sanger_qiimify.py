'''
Created on Dec 3, 2012

@author: Shareef Dabdoub
'''
import os
import unittest

import sanger_qiimify as sq

class BarcodeTest(unittest.TestCase):
    def test_generate_barcodes(self):
        num_codes = 50
        bc = sq.generate_barcodes(num_codes)
        self.assertEqual(len(bc), num_codes)


class MappingFileTest(unittest.TestCase):
    def setUp(self):
        self.sampleIDs = ['S1', 'S2', 'S3']
        self.barcodes = ['AAACCCTTTGGG', 'CCCTTTGGGAAA', 'TTTGGGAAACCC']
        
        self.mapFN = 'test_write_mapping_file__test_file'
        self.mapF = open(self.mapFN,'w')
        self.sampleMap = sq.write_mapping_file(self.mapF, self.sampleIDs, 
                                               self.barcodes)
        self.addCleanup(os.remove, self.mapFN)
        self.addCleanup(self.mapF.close)
        
    def tearDown(self):
        pass
        
    def parseMapFile(self, mapF):
        """
        Opens a QIIME mapping file and stores the contents in a dictionary
        keyed on SampleID (default) or a user-supplied one. The only 
        required fields are SampleID, BarcodeSequence, LinkerPrimerSequence 
        (in that order), and Description (which must be the final field).
        
        :@param mapFN: Full path to the map file
        
        Example data:
        #SampleID    BarcodeSequence    LinkerPrimerSequence    Treatment    SampleType    DOB    Description
        11.V13    ACGCTCGACA    GTTTGATCCTGGCTCAG    V13    Rat_Oral_Disease    111111    Rat_Oral
        """
        m = {}
        
        with mapF:
            for line in mapF:
                if line.startswith('#') or not line:
                    continue
                line = line.strip().split('\t')
                m[line[0]] = line
                
        return m
    
    def test_write_mapping_file__check_hash(self):
        self.assertIn('S1', self.sampleMap, 
                      'Sample ID not recorded in map')
        
        self.assertEqual(self.sampleMap['S1'].barcode, self.barcodes[0], 
                         'Sample ID/barcode match incorrect in returned map')

    def test_write_mapping_file__check_file_closed(self):        
        self.assertTrue(self.mapF.closed, 'Mapping file never closed')
    
    def test_write_mapping_file__check_num_entries(self):
        with open(self.mapFN, 'rU') as mapF:
            map_len = len(self.parseMapFile(mapF))
            errmsg = ('The number of sample IDs is not preserved in the '
                      'mapping file') 
            self.assertEqual(map_len, len(self.sampleIDs), errmsg)
            
    def test_write_mapping_file__check_num_fields(self):
        with open(self.mapFN, 'rU') as mapF:
            m = self.parseMapFile(mapF)
            item = m[m.keys()[0]]
            errmsg = ('There are fewer than the required 4 fields in the ' 
                      'mapping file')
            self.assertGreaterEqual(len(item), 4, errmsg)
            
    def test_write_mapping_file__check_contents(self):
        with open(self.mapFN, 'rU') as mapF:
            m = self.parseMapFile(mapF)
            self.assertTrue(self.sampleIDs[0] in m, 'Sample ID S1 not found')
            errmsg = 'Sample ID/barcode match incorrect in mapping file'
            self.assertEqual(m[self.sampleIDs[0]][1], self.barcodes[0], errmsg)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()