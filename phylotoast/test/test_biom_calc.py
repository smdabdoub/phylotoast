#!/usr/bin/env python
"""
:Author: Akshay Paropkari
:Date Created: 10/15/2014
:Abstract: Automated Tests for BIOM calculation.
"""
import math
import unittest
from phylotoast import biom_calc as bc
from biom import load_table


class biom_calc_Test(unittest.TestCase):
    def setUp(self):
        """
        Initializing BIOM format file.
        """
        self.biomf = load_table("phylotoast/test/test.biom")

    def test_relative_abundance(self):
        """
        Testing relative abundance() function of biom.calc.py.

        :return: Returns OK, if testing goal is achieved, otherwise raises
                error.
        """
        self.result = bc.relative_abundance(self.biomf)

        # List containing manual calculations
        hand_calc = [0.181818182, 0, 0.545454545, 0.272727273, 0]

        # List containing the calculated relative abundance values
        func_calc = self.result["S4"].values()

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
        self.result = bc.mean_otu_pct_abundance(self.rel_a, ["GG_OTU_1", "GG_OTU_2"])

        # Obtaining lists of function calculations and manual hand calculations
        func_calc = self.result.values()

        # list containing hand calculated relative abundance values
        hand_calc = [(0.192307692 + 0.161290323 + 0.111111111 + 0.181818182 +
                      0.086956522 + 0.333333333 + 0.071428571 + 0.230769231 +
                      0 + 0.083333333) / 10,
                     (0.076923077 + 0.258064516 + 0.222222222 + 0 + 0.260869565 +
                      0.148148148 + 0.178571429 + 0.192307692 + 0.111111111 + 0.125) / 10]

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
        self.result2 = bc.raw_abundance(self.biomf, sampleIDs=["S2", "S5"])
        self.result3 = bc.raw_abundance(self.biomf, sampleIDs=["S1", "S4"],
                                        sample_abd=False)

        # Lists containing hand and function calculated values.
        hand_calc = {"GG_OTU_1": 35., "GG_OTU_2": 38., "GG_OTU_3": 54.,
                     "GG_OTU_4": 54., "GG_OTU_5": 42.}
        hand_calc1 = {"S9": 9.0, "S8": 26.0, "S3": 18.0, "S2": 31.0, "S1": 26.0,
                      "S10": 24.0, "S7": 28.0, "S6": 27.0, "S5": 23.0, "S4": 11.0}
        hand_calc2 = {"S2": 31.0, "S5": 23.0}
        hand_calc3 = {"GG_OTU_1": 7., "GG_OTU_2": 2., "GG_OTU_3": 11., "GG_OTU_4": 12.,
                      "GG_OTU_5": 5.}

        # Testing validity of raw_abundance() function.
        self.assertDictEqual(hand_calc, self.result,
                             msg="Raw abundances not calculated accurately.")
        self.assertDictEqual(hand_calc1, self.result1,
                             msg="Raw abundances not calculated accurately.")
        self.assertDictEqual(self.result2, hand_calc2,
                             msg="Abundances not calculated for SampleID's")
        self.assertDictEqual(self.result3, hand_calc3,
                             msg="Abundances not calculated for OTUID's")

    def test_transform_raw_abundance(self):
        """
        Testing transform_raw_abundance() function of biom_calc.py.

        :return: Returns OK if testing goal is achieved, otherwise raises
                 error.
        """
        self.result = bc.transform_raw_abundance(self.biomf, sample_abd=False)
        self.result1 = bc.transform_raw_abundance(self.biomf, fn=math.sqrt)

        # Obtaining manual calculations for comparison testing
        hand_calc = {"GG_OTU_1": 1.544068044, "GG_OTU_2": 1.579783597,
                     "GG_OTU_3": 1.73239376, "GG_OTU_4": 1.73239376,
                     "GG_OTU_5": 1.62324929}
        hand_calc1 = {"S9": 3.0, "S8": 5.09901951, "S3": 4.24264069, "S2": 5.56776436,
                      "S1": 5.09901951, "S10": 4.89897949, "S7": 5.29150262,
                      "S6": 5.19615242, "S5": 4.79583152, "S4": 3.31662479}

        # Testing the validity of transform function
        for hand, func in zip(hand_calc.values(), self.result.values()):
            self.assertAlmostEqual(
                func, hand,
                msg="Raw abundance transformation not computed accurately."
            )
        for hand1, func1 in zip(hand_calc1.values(), self.result1.values()):
            self.assertAlmostEqual(
                func1, hand1,
                msg="Raw abundance transformation not computed accurately."
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
        hand_calc = [0.440510663, 0, 0.830915552, 0.549467245, 0]
        func_calc = self.result2["S4"].values()

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
