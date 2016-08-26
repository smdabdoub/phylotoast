#!/usr/bin/env python
"""
:Author: Akshay Paropkari
:Date Created: 10/22/2014
:Abstract: Automated Tests for OTU calculations.
"""
import biom
import unittest
from phylotoast import otu_calc as oc


class otu_calc_Test(unittest.TestCase):

    def setUp(self):
        """
        Setting up the test module. Initializing BIOM format file.
        """
        pass

    def test_otu_name(self):
        """
        Testing the otu_name() function of otu_calc.py.

        :return: Returns OK if the test goals were achieved, otherwise
                 raises error.
        """
        self.taxa = {
            "Unclassified_Methanosarcinales":
            ["k__Archaea", "p__Euryarchaeota", "c__Methanomicrobia",
             "o__Methanosarcinales", "f__",
             "g__", "s__concilii"],
            "Campylobacter_gracilis":
            ["k__Bacteria", "p__Proteobacteria", "c__Epsilonproteobacteria",
             "o__Campylobacterales", "f__Campylobacteraceae", "g__Campylobacter",
             "s__gracilis"],
            "Escherichia_spp.":
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"],
            "Fusobacterium_nucleatum":
            ["k__Bacteria", "p__Fusobacteria", "c__Fusobacteria", "o__",
             "f__", "g__Fusobacterium", "s__nucleatum"],
            "Fusobacterium_spp.":
            ["k__Bacteria", "p__Fusobacteria", "c__Fusobacteria", "o__",
            "f__", "g__Fusobacterium", "s__"]
        }

        for expected, test in self.taxa.items():
            self.result = oc.otu_name(test)

            # Testing the validity of the otu_name() function
            self.assertEqual(
                self.result, expected,
                msg="Error!\nExpected result: {}.\notu_name() result: {}".
                    format(expected, self.result)
            )

    def test_load_core_file(self):
        """
        Testing the load_core_file() function of otu_calc.py
        :return: Returns OK if the test goals were achieved, otherwise
                 raises error.
        """
        result = oc.load_core_file("phylotoast/test/test_core.txt")
        hand_calc = {"Actinomyces_spp.", "Campylobacter_spp.", "Capnocytophaga_spp.",
                     "Catonella_spp.", "Corynebacterium_spp.", "Dialister_spp.",
                     "Eikenella_spp.", "Filifactor_spp.", "Fusobacterium_spp.",
                     "Gemella_spp.", "Granulicatella_spp.", "Kingella_spp.",
                     "Leptotrichia_spp.", "Megasphaera_spp.", "Parvimonas_spp.",
                     "Prevotella_melaninogenica", "Prevotella_spp.", "Selenomonas_noxia",
                     "Selenomonas_spp.", "Streptococcus_anginosus", "Streptococcus_equi",
                     "Streptococcus_infantis", "Streptococcus_spp.",
                     "Unclassified_Lachnospiraceae", "Unclassified_TM7-3",
                     "Unclassified_[Mogibacteriaceae]", "Veillonella_dispar",
                     "Veillonella_parvula", "Veillonella_spp."}

        # Testing if all core OTU's samples were in the output.
        self.assertSetEqual(
            result, hand_calc,
            msg="Error! Genus-species names not calculated as expected."
        )

    def test_assign_otu_membership(self):
        """
        Testing assign_otu_membership() function of otu_calc.py.

        :return: Returns OK if the test goals were achieved, otherwise
                 raises error.
        """
        self.biomf = biom.load_table("phylotoast/test/test.biom")
        self.result = oc.assign_otu_membership(self.biomf)

        # Obtaining the values to be tested
        hand_calc = {"S9": ["GG_OTU_2", "GG_OTU_3", "GG_OTU_5"],
                     "S3": ["GG_OTU_1", "GG_OTU_2", "GG_OTU_4", "GG_OTU_5"],
                     "S6": ["GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"]}

        # Testing the validity of assign_otu_membership() function
        for sid in ["S3", "S6", "S9"]:
            self.assertListEqual(
                sorted(hand_calc[sid]), sorted(self.result[sid]),
                msg="Error! OTU membership calculations are inaccurate!"
            )

    def tearDown(self):
        """
        Tearing down of this unittest framework.
        """
        pass

if __name__ == "__main__":
    unittest.main()
