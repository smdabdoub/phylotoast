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

    def test_arcsine_sqrt_transform(self):
            """
            Testing arcsine_sqrt_transform() function of biom_calc.py.

            :return: Returns OK if testing goal is achieved, otherwise raises
                     error.
            """
            self.result1 = bc.relative_abundance(self.biomf)
            self.result2 = bc.arcsine_sqrt_transform(self.result1)

            # Obtaining results to compare.
            hand_calc = {"S1": {"GG_OTU_1": 0.453961252, "GG_OTU_2": 0.281034902,
                                "GG_OTU_3": 0.453961252, "GG_OTU_4": 0.629014802,
                                "GG_OTU_5": 0.453961252},
                         "S10": {"GG_OTU_1": 0.292842772, "GG_OTU_2": 0.361367124,
                                 "GG_OTU_3": 0.420534335, "GG_OTU_4": 0.615479709,
                                 "GG_OTU_5": 0.570510448},
                         "S2": {"GG_OTU_1": 0.413273808, "GG_OTU_2": 0.532861869,
                                "GG_OTU_3": 0.532861869, "GG_OTU_4": 0.532861869,
                                "GG_OTU_5": 0.256813917},
                         "S3": {"GG_OTU_1": 0.339836909, "GG_OTU_2": 0.490882678,
                                "GG_OTU_3": 0, "GG_OTU_4": 0.555121168,
                                "GG_OTU_5": 0.673351617},
                         "S4": {"GG_OTU_1": 0.440510663, "GG_OTU_2": 0,
                                "GG_OTU_3": 0.830915552, "GG_OTU_4": 0.549467245,
                                "GG_OTU_5": 0},
                         "S5": {"GG_OTU_1": 0.299334026, "GG_OTU_2": 0.53606149,
                                "GG_OTU_3": 0.584373897, "GG_OTU_4": 0.485049787,
                                "GG_OTU_5": 0.36950894},
                         "S6": {"GG_OTU_1": 0.615479709, "GG_OTU_2": 0.395099667,
                                "GG_OTU_3": 0.575591472, "GG_OTU_4": 0.444859969,
                                "GG_OTU_5": 0.1936583},
                         "S7": {"GG_OTU_1": 0.270549763, "GG_OTU_2": 0.436286927,
                                "GG_OTU_3": 0.387596687, "GG_OTU_4": 0.563942641,
                                "GG_OTU_5": 0.602794553},
                         "S8": {"GG_OTU_1": 0.501093013, "GG_OTU_2": 0.453961252,
                                "GG_OTU_3": 0.588002604, "GG_OTU_4": 0.346579954,
                                "GG_OTU_5": 0.403057074},
                         "S9": {"GG_OTU_1": 0, "GG_OTU_2": 0.339836909,
                                "GG_OTU_3": 0.729727656, "GG_OTU_4": 0,
                                "GG_OTU_5": 0.729727656}}

            # Testing validity of the transforms.
            for sid in sorted(hand_calc.keys()):
                for otuid in sorted(hand_calc[sid].keys()):
                    self.assertAlmostEqual(
                        hand_calc[sid][otuid], self.result2[sid][otuid],
                        msg="Arcsine squareroot transformation was not accurate."
                    )

    def test_mean_otu_pct_abundance(self):
        """
        Testing mean_otu_pct_abundance() function of biom_calc.py.

        :return: Returns OK, if testing goal was achieved, otherwise raises error.
        """
        self.rel_a = bc.relative_abundance(self.biomf)
        self.result = bc.mean_otu_pct_abundance(self.rel_a, ["GG_OTU_1", "GG_OTU_2"])

        # list containing hand calculated relative abundance values
        hand_calc = {"GG_OTU_1": 14.52348298, "GG_OTU_2": 15.73217761,
                     "GG_OTU_3": 26.58131438, "GG_OTU_4": 22.91732137,
                     "GG_OTU_5": 20.24570366}

        # Testing the validity of the calculations of mean_otu_pct_abundance().
        for oid in ["GG_OTU_1", "GG_OTU_2"]:
            self.assertAlmostEqual(
                hand_calc[oid], self.result[oid],
                msg="Mean OTU percent abundance not calculated accurately."
            )

    def test_MRA(self):
        """
        Testing mean relative abundance calculation, MRA() function of biom_calc.py.

        :return: Returns OK, if testing goal was achieved, otherwise raises error.
        """
        self.result = bc.MRA(self.biomf)

        # Obtaining lists of function calculations and manual hand calculations
        hand_calc = {"GG_OTU_1": 14.52348298, "GG_OTU_2": 15.73217761,
                     "GG_OTU_3": 26.58131438, "GG_OTU_4": 22.91732137,
                     "GG_OTU_5": 20.24570366}

        # Testing the validity of the calculations of mean_otu_pct_abundance()
        for oid in hand_calc.keys():
            self.assertAlmostEqual(
                hand_calc[oid], self.result[oid],
                msg="Mean relative OTU percent abundance (MRA) not calculated accurately."
            )

    def test_raw_abundance(self):
        """
        Testing raw_abundance() function of biom_calc.py.

        :return: Returns OK, if testing goal is achieved, otherwise raises error.
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

    def test_relative_abundance(self):
        """
        Testing relative abundance() function of biom.calc.py.

        :return: Returns OK, if testing goal is achieved, otherwise raises error.
        """
        self.result = bc.relative_abundance(self.biomf)

        # List containing manual calculations
        hand_calc = {"S1": {"GG_OTU_1": 0.192307692, "GG_OTU_2": 0.076923077,
                            "GG_OTU_3": 0.192307692, "GG_OTU_4": 0.346153846,
                            "GG_OTU_5": 0.192307692},
                     "S10": {"GG_OTU_1": 0.083333333, "GG_OTU_2": 0.125,
                             "GG_OTU_3": 0.166666667, "GG_OTU_4": 0.333333333,
                             "GG_OTU_5": 0.291666667},
                     "S2": {"GG_OTU_1": 0.161290323, "GG_OTU_2": 0.258064516,
                            "GG_OTU_3": 0.258064516, "GG_OTU_4": 0.258064516,
                            "GG_OTU_5": 0.064516129},
                     "S3": {"GG_OTU_1": 0.111111111, "GG_OTU_2": 0.222222222,
                            "GG_OTU_3": 0.0, "GG_OTU_4": 0.277777778,
                            "GG_OTU_5": 0.388888889},
                     "S4": {"GG_OTU_1": 0.181818182, "GG_OTU_2": 0.0,
                            "GG_OTU_3": 0.545454545, "GG_OTU_4": 0.272727273,
                            "GG_OTU_5": 0.0},
                     "S5": {"GG_OTU_1": 0.086956522, "GG_OTU_2": 0.260869565,
                            "GG_OTU_3": 0.304347826, "GG_OTU_4": 0.217391304,
                            "GG_OTU_5": 0.130434783},
                     "S6": {"GG_OTU_1": 0.333333333, "GG_OTU_2": 0.148148148,
                            "GG_OTU_3": 0.296296296, "GG_OTU_4": 0.185185185,
                            "GG_OTU_5": 0.037037037},
                     "S7": {"GG_OTU_1": 0.071428571, "GG_OTU_2": 0.178571429,
                            "GG_OTU_3": 0.142857143, "GG_OTU_4": 0.285714286,
                            "GG_OTU_5": 0.321428571},
                     "S8": {"GG_OTU_1": 0.230769231, "GG_OTU_2": 0.192307692,
                            "GG_OTU_3": 0.307692308, "GG_OTU_4": 0.115384615,
                            "GG_OTU_5": 0.153846154},
                     "S9": {"GG_OTU_1": 0.0, "GG_OTU_2": 0.111111111,
                            "GG_OTU_3": 0.444444444, "GG_OTU_4": 0.0,
                            "GG_OTU_5": 0.444444444}}

        # Testing the validity of relative_abundance() function.
        for sid in sorted(hand_calc.keys()):
            for otuid in sorted(hand_calc[sid].keys()):
                self.assertAlmostEqual(
                    hand_calc[sid][otuid], self.result[sid][otuid],
                    msg="Relative abundances not calculated accurately."
                )

    def test_transform_raw_abundance(self):
        """
        Testing transform_raw_abundance() function of biom_calc.py.

        :return: Returns OK if testing goal is achieved, otherwise raises error.
        """
        self.result = bc.transform_raw_abundance(self.biomf, sample_abd=False)
        self.result1 = bc.transform_raw_abundance(self.biomf, fn=math.sqrt)

        # Obtaining manual calculations for comparison testing
        hand_calc = {"GG_OTU_1": 1.544068044, "GG_OTU_2": 1.579783597,
                     "GG_OTU_3": 1.73239376, "GG_OTU_4": 1.73239376,
                     "GG_OTU_5": 1.62324929}
        hand_calc1 = {"S1": 5.099019514, "S2": 5.567764363, "S3": 4.242640687,
                      "S4": 3.31662479, "S5": 4.795831523, "S6": 5.196152423,
                      "S7": 5.291502622, "S8": 5.099019514, "S9": 3, "S10": 4.898979486}

        # Testing the validity of transform function
        for oid in hand_calc.keys():
            self.assertAlmostEqual(
                hand_calc[oid], self.result[oid],
                msg="Raw abundance transformation not computed accurately (Test1)."
            )
        for sid in hand_calc1.keys():
            self.assertAlmostEqual(
                hand_calc1[sid], self.result1[sid],
                msg="Raw abundance transformation not computed accurately (Test2)."
            )

    def tearDown(self):
        """
        No particular event to clean, delete or close in testing biom_calc.py.
        """
        pass

if __name__ == "__main__":
    unittest.main()
