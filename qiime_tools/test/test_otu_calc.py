#!/usr/bin/env python
"""
:Author: Akshay Paropkari

:Contact: paropkari.1@gmail.com

:Date Created: 10/22/2014

:Abstract: Automated Tests for OTU calculations.
"""
import unittest
import json
import tempfile
import sys
from qiime_tools import otu_calc as oc
from qiime_tools import biom_calc as bc


class otu_calc_Test(unittest.TestCase):

    def setUp(self):
        """
        Setting up the test module. Initializing BIOM format file.
        """
        self.biom_text = """{
                     "id":null,
                     "format": "Biological Observation Matrix 0.9.1-dev",
                     "format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
                     "type": "OTU table",
                     "generated_by": "QIIME revision 1.4.0-dev",
                     "date": "2011-12-19T19:00:00",
                     "rows":[
                        {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}},
                        {"id":"GG_OTU_2", "metadata":{"taxonomy":["k__Bacteria", "p__Cyanobacteria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichospermum", "s__"]}},
                        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Archaea", "p__Euryarchaeota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "g__Methanosarcina", "s__"]}},
                        {"id":"GG_OTU_4", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanaerobium", "s__Halanaerobiumsaccharolyticum"]}},
                        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}}
                        ],
                     "columns":[
                        {"id":"Sample1", "metadata":{
                                                 "BarcodeSequence":"CGCTTATCGAGA",
                                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                                 "BODY_SITE":"gut",
                                                 "Description":"human gut"}},
                        {"id":"Sample2", "metadata":{
                                                 "BarcodeSequence":"CATACCAGTAGC",
                                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                                 "BODY_SITE":"gut",
                                                 "Description":"human gut"}},
                        {"id":"Sample3", "metadata":{
                                                 "BarcodeSequence":"CTCTCTACCTGT",
                                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                                 "BODY_SITE":"gut",
                                                 "Description":"human gut"}},
                        {"id":"Sample4", "metadata":{
                                                 "BarcodeSequence":"CTCTCGGCCTGT",
                                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                                 "BODY_SITE":"skin",
                                                 "Description":"human skin"}},
                        {"id":"Sample5", "metadata":{
                                                 "BarcodeSequence":"CTCTCTACCAAT",
                                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                                 "BODY_SITE":"skin",
                                                 "Description":"human skin"}},
                        {"id":"Sample6", "metadata":{
                                                 "BarcodeSequence":"CTAACTACCAAT",
                                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                                 "BODY_SITE":"skin",
                                                 "Description":"human skin"}}
                                ],
                     "matrix_type": "sparse",
                     "matrix_element_type": "int",
                     "shape": [5, 6],
                     "data":[[0,2,1],
                             [1,0,5],
                             [1,1,1],
                             [1,3,2],
                             [1,4,3],
                             [1,5,1],
                             [2,2,1],
                             [2,3,4],
                             [2,5,2],
                             [3,0,2],
                             [3,1,1],
                             [3,2,1],
                             [3,5,1],
                             [4,1,1],
                             [4,2,1]
                            ]
                        }"""
        self.biom = json.loads(self.biom_text)
        self.row = self.biom['rows']

    def test_fuzzy_lookup(self):
        """
        Testing fuzzy_lookup() function of otu_calc.py.

        :return:Returns OK if the test goals were achieved, otherwise
                raises error.
        """
        self.result1 = oc.assign_otu_membership(self.biom)['Sample2']
        self.result = oc.fuzzy_lookup(self.result1, ['Escherichia_spp.'])

        # Obtaining manual result to be compared against.
        hand_calc = 'Escherichia_spp.'

        # Testing the validity of fuzzy_lookup() function.
        self.assertIn(
            hand_calc, str(self.result1),
            msg='Error! Output is not a subset of input.'
        )

    def test_sdi(self):
        """
        Testing sdi() function of otu_calc.py.

        :return:Returns OK if the test goals were achieved, otherwise
                raises error.
        """
        self.result1 = oc.assign_otu_membership(self.biom)
        fset = self.result1['Sample4']
        self.result = oc.sdi(fset)

        # Manual calculation of SDI
        hand_calc = 0.636513937

        # Testing the validity of sdi() function
        self.assertAlmostEqual(
            self.result, hand_calc, places=6,
            msg='Error! SDI was calculated inacccurately.'
        )

    def test_otu_name_biom(self):
        """
        Testing otu_name_biom() function of otu_calc.py.

        :return:Returns OK if the test goals were achieved, otherwise
                raises error.
        """
        self.result = oc.otu_name_biom(self.row[3])
        hand_calc = 'Halanaerobium_Halanaerobiumsaccharolyticum'

        # Testing the validity of otu_name_biom() function
        self.assertEqual(
            self.result, hand_calc,
            msg='Error! otu_name_biom() output does not match manually identified name.'
        )

    def test_otu_name(self):
        """
        Testing otu_name() function of otu_calc.py.

        :return:Returns OK if the test goals were achieved, otherwise
                raises error.
        """
        self.tax = [
            "k__Archaea", "p__Euryarchaeota", "c__Methanomicrobia",
            "o__Methanosarcinales", "f__",
            "g__", "s__"
            ]
        self.result = oc.otu_name(self.tax)
        hand_calc = 'Unclassified_Methanosarcinales'

        # Testing the validity of the otu_name() function
        self.assertEqual(
            self.result, hand_calc,
            msg='Error! The output is not as expected.'
            )

    def test_load_core_file(self):
        """
        Testing load_core_file() function of otu_calc.py

        :return:Returns OK if the test goals were achieved, otherwise
                raises error.
        """
        self.result = oc.load_core_file(
            '/Users/akshay/dev/core_otus_example.txt'
            )
        length = len(self.result)

        # Testing if all core OTU's samples were in the output.
        self.assertEqual(
            length, 17,
            msg='Error! Output does not contain all core OTU in samples.'
        )

    def test_assign_otu_membership(self):
        """
        Testing assign_otu_membership() function of otu_calc.py.

        :return: Returns OK if the test goals were achieved, otherwise
                raises error.
        """
        self.result = oc.assign_otu_membership(self.biom)

        # Obtaining the values to be tested
        result1 = bc.relative_abundance(self.biom, ['Sample1'])
        hand_calc = result1.values()[0].values()
        func_calc = [0.714286, 0.285714]

        # Testing the validity of assign_otu_membership() function
        for hand, func in zip(hand_calc, func_calc):
            self.assertAlmostEqual(
                hand, func, places=5,
                msg='Error! OTU membership calculations are inaccurate!'
            )

    def test_print_membership(self):
        """
        Testing print_membership() function of otu_calc.py. Output of this
        function is written to output.txt file. The contents of output.txt
        file are accessed and used to compare with expected output.

        :return:Returns OK if the test goals were achieved, otherwise
                raises error.
        """
        self.result1 = oc.assign_otu_membership(self.biom)

        # Writing the output of print_membership() to output.txt file
        saveout = sys.stdout
        with open('/Users/akshay/Desktop/docs/output.txt', 'w') as f:
            sys.stdout = f
            self.result = oc.print_membership(self.result1['Sample1'])
        sys.stdout = saveout

        # Obtaining the function output from output.txt file.
        func_calc = []
        with open('/Users/akshay/Desktop/docs/output.txt', 'r') as out:
            while True:
                statement = out.readline()
                func_calc.append(statement)
                if not statement:
                    break

        # Writing manual outputs to compare
        hand_calc = [
            'Halanaerobium_Halanaerobiumsaccharolyticum: 28.57%\n',
            'Dolichospermum_spp.: 71.43%\n'
            ]

        # Testing the validity of print_membership() function
        self.assertListEqual(
            func_calc[0:2], hand_calc,
            msg='Error! Output does not print out the right information correctly.'
        )

    def tearDown(self):
        """
        Tearing down of this unittest framework.
        """
        # Erasing data from output.txt file used in testing print_membership
        # function.
        with open('/Users/akshay/Desktop/docs/output.txt', 'w') as f:
            f.write('')

if __name__ == '__main__':
    unittest.main()
