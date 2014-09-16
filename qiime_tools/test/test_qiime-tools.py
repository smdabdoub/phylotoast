'''
Automated tests for QIIME-tools

@author: Shareef Dabdoub
'''
import unittest

import qiime_tools
from qiime_tools import iTol
from qiime_tools import sanger_qiimify


#################
# iTol.py tests #
#################
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
        expected = {'healthy_smoker':iTol.DataCategory({'SID1'},{}),
                    'pm_smoker':iTol.DataCategory({'SID2'},{}),
                    'pi_smoker':iTol.DataCategory({'SID3'},{}),
                    'healthy_ns':iTol.DataCategory({'SID5'},{}),
                    'pm_ns':iTol.DataCategory({'SID7'},{}),
                    'pi_ns':iTol.DataCategory({'SID8'},{})}
        self.assertDictContainsSubset(expected, cats)



#####################
# sanger_qiimify.py #
#####################
class BarcodeTest(unittest.TestCase):
    def test_generate_barcodes(self):
        num_codes = 50
        bc = sanger_qiimify.generate_barcodes(num_codes)
        self.assertEqual(len(bc), num_codes)


class MappingFileTest(unittest.TestCase):
    def setUp(self):
        self.sampleIDs = ['S1', 'S2', 'S3']
        self.barcodes = ['AAACCCTTTGGG', 'CCCTTTGGGAAA', 'TTTGGGAAACCC']

        self.mapFN = 'test_write_mapping_file__test_file'
        self.mapF = open(self.mapFN,'w')
        self.sampleMap = sanger_qiimify.write_mapping_file(self.mapF, self.sampleIDs,
                                                           self.barcodes)
        self.addCleanup(os.remove, self.mapFN)
        self.addCleanup(self.mapF.close)

    def tearDown(self):
        pass

    def parse_map_file(self, mapF):
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
            map_len = len(self.parse_map_file(mapF))
            errmsg = ('The number of sample IDs is not preserved in the '
                      'mapping file')
            self.assertEqual(map_len, len(self.sampleIDs), errmsg)

    def test_write_mapping_file__check_num_fields(self):
        with open(self.mapFN, 'rU') as mapF:
            m = self.parse_map_file(mapF)
            item = m[m.keys()[0]]
            errmsg = ('There are fewer than the required 4 fields in the '
                      'mapping file')
            self.assertGreaterEqual(len(item), 4, errmsg)

    def test_write_mapping_file__check_contents(self):
        with open(self.mapFN, 'rU') as mapF:
            m = self.parse_map_file(mapF)
            self.assertTrue(self.sampleIDs[0] in m, 'Sample ID S1 not found')
            errmsg = 'Sample ID/barcode match incorrect in mapping file'
            self.assertEqual(m[self.sampleIDs[0]][1], self.barcodes[0], errmsg)



#################
# prune_otus.py #
#################
class PruneOTUsTest(unittest.TestCase):

    def setUp(self):
        self.otus = {'111111': ['SID1_10', 'SID1_22']}


    def tearDown(self):
        pass


    def test_filter_by_sample_pct(self):
        pass

    def test_filter_by_sequence_pct(self):
        pass
