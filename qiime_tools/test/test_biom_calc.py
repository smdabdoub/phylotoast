#!/usr/bin/env python
"""
:Author: Akshay Paropkari

:Contact: paropkari.1@gmail.com

:Date Created: 10/15/2014

:Abstract: Automated Tests for BIOM calculation.
"""
import unittest
import json
import math
from qiime_tools import biom_calc as bc


class biom_calc_Test(unittest.TestCase):

    def setUp(self):
        """
        Initializing BIOM format file.
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

    def test_relative_abundance(self):
        """
        Testing relative abundance() function of biom.calc.py.

        :return: Returns OK, if testing goal is achieved, otherwise raises
                error.
        """
        sample = 'Sample3'
        self.result = bc.relative_abundance(self.biom, sample)

        # List containing manual calculations
        hand_calc = [1/4.0, 1/4.0, 1/4.0, 1/4.0]

        # Obtaining list of function calculated relative abundance for sample
        result1 = self.result.values()      # result1 is a list
        result2 = result1[0]                # result2 is a dict
        # list containing the calculated relative abundance values
        func_calc = result2.values()

        # Testing the validity of relative_abundance() function.
        for hand, res in zip(hand_calc, func_calc):
            self.assertAlmostEqual(
                hand, res,
                msg='Relative abundances not calculated accurately.'
                )

    def test_mean_otu_pct_abundance(self):
        """
        Testing mean_otu_pct_abundance() function of biom_calc.py.

        :return: Returns OK, if testing goal was achieved, otherwise raises
                error.
        """
        self.rel_a = bc.relative_abundance(self.biom, ['Sample4'])
        self.result = bc.mean_otu_pct_abundance(
            self.rel_a, ['GG_OTU_2', 'GG_OTU_3']
            )

        # Obtaining lists of function calculations and manual hand calculations
        func_calc = self.result.values()
        result1 = self.rel_a.values()       # result1 is a list
        result2 = result1[0]                # result2 is a dict
        # list containing hand calculated relative abundance values
        hand_calc = result2.values()

        # Testing the validity of the calculations of mean_otu_pct_abundance().
        for hand, res in zip(hand_calc, func_calc):
            self.assertAlmostEqual(
                hand*100, res,
                msg='Mean OTU not calculated accurately.'
                )

    def test_MRA(self):
        """
        Testing mean relative abundance calculation, MRA() function
        of biom_calc.py.

        :return: Returns OK, if testing goal was achieved, otherwise
            raises error.
        """
        self.result = bc.MRA(self.biom, 'Sample4')
        self.mean_otu = bc.mean_otu_pct_abundance(
            bc.relative_abundance(self.biom, ['Sample4']),
            ['GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3', 'GG_OUTU4', 'GG_OTU_5']
            )

        # Obtaining lists of function calculations and
        # manual hand calculations
        func_calc = self.result.values()
        hand_calc = self.mean_otu.values()

        # Testing the validity of the calculations of mean_otu_pct_abundance().
        for hand, res in zip(hand_calc, func_calc):
            self.assertAlmostEqual(
                hand, res,
                msg='Mean OTU not calculated accurately.'
                )

    def test_raw_abundance(self):
        """
        Testing raw_abundance() function of biom_calc.py.

        :return: Returns OK, if testing goal is achieved, otherwise raises
                 error.
        """
        self.result = bc.raw_abundance(self.biom, sample_abd=False)
        self.result1 = bc.raw_abundance(self.biom)

        # Lists containing hand and function calculated values.
        hand_calc = [7, 3, 4, 6, 3, 4]
        hand_calc1 = [1, 12, 7, 5, 2]
        hand_calc2 = ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5',
                      'Sample6']
        hand_calc3 = ['GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3', 'GG_OTU_4',
                      'GG_OTU_5']

        # Testing validity of raw_abundance() function.
        self.assertItemsEqual(hand_calc, self.result.values(),
                              msg='Raw abundances not calculated accurately.')
        self.assertItemsEqual(hand_calc1, self.result1.values(),
                              msg='Raw abundances not calculated accurately.')
        self.assertItemsEqual(self.result.keys(), hand_calc2,
                              msg='Abundances not calculated for SampleID\'s')
        self.assertItemsEqual(self.result1.keys(), hand_calc3,
                              msg='Abundances not calculated for OTUID\'s')

    def test_transform_raw_abundance(self):
        """
        Testing transform_raw_abundance() function of biom_calc.py.

        :return: Returns OK if testing goal is achieved, otherwise raises
                 error.
        """
        self.result = bc.transform_raw_abundance(
            self.biom, sample_abd=False
            )
        self.result1 = bc.raw_abundance(self.biom, sample_abd=False)

        # Obtaining manual calculations for comparison testing
        hand_calc = []
        for num in self.result1.values():
            hand_calc.append(math.log10(float(num)))

        # Testing the validity of transform function
        self.assertAlmostEqual(
            self.result.values(), hand_calc, places=10,
            msg='Function did not calculate the transformation accurately.'
        )

    def tearDown(self):
        """
        No particular event to clean, delete or close in testing biom_calc.py.
        """
        pass

if __name__ == '__main__':
    unittest.main()
