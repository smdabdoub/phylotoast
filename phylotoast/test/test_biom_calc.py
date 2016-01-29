#!/usr/bin/env python
"""
:Author: Akshay Paropkari

:Date Created: 10/15/2014

:Abstract: Automated Tests for BIOM calculation.
"""
import sys
import json
import math
import unittest
from phylotoast import biom_calc as bc
import numpy as np
import biom

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
        sampleids = [self.biom["columns"][i]["id"] for i in range(6)]
        otuids = [self.biom["rows"][j]["id"] for j in range(5)]
        data = np.arange(30).reshape(5, 6)
        self.biomf = biom.table.Table(data, otuids, sampleids)

    def test_relative_abundance(self):
        """
        Testing relative abundance() function of biom.calc.py.

        :return: Returns OK, if testing goal is achieved, otherwise raises
                error.
        """
        self.result = bc.relative_abundance(self.biomf)

        # List containing manual calculations
        hand_calc = [0.02857142857, 0.11428571429, 0.2, 0.28571428571,
                     0.37142857143]

        # List containing the calculated relative abundance values
        func_calc = self.result["Sample3"].values()

        # Testing the validity of relative_abundance() function.
        for hand, res in zip(hand_calc, func_calc):
            self.assertAlmostEqual(
                hand, res,
                msg="Relative abundances not calculated accurately."
                )

    def test_mean_otu_pct_abundance(self):
        """
        Testing mean_otu_pct_abundance() function of biom_calc.py.

        :return: Returns OK, if testing goal was achieved, otherwise raises
                error.
        """
        self.rel_a = bc.relative_abundance(self.biomf)

        self.result = bc.mean_otu_pct_abundance(
            self.rel_a, ["GG_OTU_1", "GG_OTU_2"]
            )

        # Obtaining lists of function calculations and manual hand calculations
        func_calc = self.result.values()

        # list containing hand calculated relative abundance values
        hand_calc = [(0 + 0.0153846153846 + 0.0285714285714 + 0.04 + 0.05 +
                      0.0588235294118)/6,
                     (0.1 + 0.107692307692 + 0.114285714286 + 0.12 + 0.125 +
                      0.129411764706)/6]

        # Testing the validity of the calculations of mean_otu_pct_abundance().
        for hand, res in zip(hand_calc, func_calc):
            self.assertAlmostEqual(
                hand*100, res,
                msg="Mean OTU not calculated accurately."
                )

    def test_MRA(self):
        """
        Testing mean relative abundance calculation, MRA() function
        of biom_calc.py.

        :return: Returns OK, if testing goal was achieved, otherwise
            raises error.
        """
        self.result = bc.MRA(self.biomf)
        self.mean_otu = bc.mean_otu_pct_abundance(
            bc.relative_abundance(self.biomf),
            ["GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"]
            )

        # Obtaining lists of function calculations and manual hand calculations
        func_calc = self.result.values()
        hand_calc = self.mean_otu.values()

        # Testing the validity of the calculations of mean_otu_pct_abundance()
        for hand, res in zip(hand_calc, func_calc):
            self.assertAlmostEqual(
                hand, res,
                msg="Mean OTU not calculated accurately."
                )

    def test_raw_abundance(self):
        """
        Testing raw_abundance() function of biom_calc.py.

        :return: Returns OK, if testing goal is achieved, otherwise raises
                 error.
        """
        self.result = bc.raw_abundance(self.biomf, sample_abd=False)
        self.result1 = bc.raw_abundance(self.biomf)
        self.result2 = bc.raw_abundance(self.biomf,
                                        sampleIDs=["Sample2", "Sample5"])
        self.result3 = bc.raw_abundance(self.biomf,
                                        sampleIDs=["Sample1", "Sample4"],
                                        sample_abd=False)

        # Lists containing hand and function calculated values.
        hand_calc = [15, 51, 87, 123, 159]
        hand_calc1 = [60, 65, 70, 75, 80, 85]
        hand_calc2 = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
                      "Sample6"]
        hand_calc3 = ["GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4",
                      "GG_OTU_5"]
        hand_calc4 = [65, 80]
        hand_calc5 = ["Sample2", "Sample5"]
        hand_calc6 = [3, 15, 27, 39, 51]
        hand_calc7 = ["GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4",
                      "GG_OTU_5"]

        # Testing validity of raw_abundance() function.
        self.assertItemsEqual(hand_calc1, sorted(self.result1.values()),
                              msg="Raw abundances not calculated accurately.")
        self.assertItemsEqual(hand_calc, self.result.values(),
                              msg="Raw abundances not calculated accurately.")
        self.assertItemsEqual(self.result.keys(), hand_calc3,
                              msg="Abundances not calculated for SampleID's")
        self.assertItemsEqual(sorted(self.result1.keys()), hand_calc2,
                              msg="Abundances not calculated for OTUID's")
        self.assertItemsEqual(sorted(self.result2.keys()), hand_calc5,
                              msg="Abundances not calculated for SampleID's")
        self.assertItemsEqual(self.result3.keys(), hand_calc7,
                              msg="Abundances not calculated for OTUID's")
        self.assertItemsEqual(sorted(self.result2.values()), hand_calc4,
                              msg="Abundances not calculated for SampleID's")
        self.assertItemsEqual(self.result3.values(), hand_calc6,
                              msg="Abundances not calculated for OTUID's")

    def test_transform_raw_abundance(self):
        """
        Testing transform_raw_abundance() function of biom_calc.py.

        :return: Returns OK if testing goal is achieved, otherwise raises
                 error.
        """
        self.result = bc.transform_raw_abundance(
            self.biomf, sample_abd=False
            )
        self.result1 = bc.raw_abundance(self.biomf, sample_abd=False)

        # Obtaining manual calculations for comparison testing
        hand_calc = [1.17609125906, 1.7075701761, 1.93951925262,
                     2.08990511144, 2.20139712432]

        # Testing the validity of transform function
        for hand, func in zip(hand_calc, self.result.values()):
            self.assertAlmostEqual(
                hand, func,
                msg="Function did not calculate the transformation accurately."
            )

    def test_arcsine_sqrt_transform(self):
        """
        Testing arcsine_sqrt_transform() function of biom_calc.py.

        :return: Returns OK if testing goal is achieved, otherwise raises
                 error.
        """
        self.result1 = bc.relative_abundance(self.biomf)
        self.result2 = bc.arcsine_sqrt_transform(self.result1)

        # Obtaining results to compare.
        hand_calc = [0, 0.32175055439, 0.463647609, 0.57963974036, 0.684719203]
        func_calc = self.result2.values()[3].values()

        # Testing validity of the transforms.
        for hand, func in zip(hand_calc, func_calc):
            self.assertAlmostEqual(
                hand, func, places=7,
                msg="Function did not calculate transformation accurately."
            )

    def tearDown(self):
        """
        No particular event to clean, delete or close in testing biom_calc.py.
        """
        pass

if __name__ == "__main__":
    unittest.main()
