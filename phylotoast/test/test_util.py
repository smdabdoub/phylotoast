"""
:Author: Akshay Paropkari

:Date Created: 10/13/2014

:Abstract: Automated tests for util.py functions.
"""
import os
import unittest
import tempfile
from collections import namedtuple
from phylotoast import util as ut


class util_Test(unittest.TestCase):

    def setUp(self):
        """
        Setting up files, or data for testing purposes.
        """
        self.map_header, self.map_data = ut.parse_map_file("phylotoast/test/test_mapping_file.txt")

    def test_parseFASTA(self):
        """
        Testing parseFASTA function.

        :return: Returns OK is test goals were achieved, otherwise raises
                 error.
        """
        FASTARecord = namedtuple("FASTA_Record", "id descr data")
        parseFASTA_result = ut.parseFASTA("phylotoast/test/test_FASTA.fna")
        manually_parsed = [FASTARecord(id="PIDF154_1", descr="HU82XDC01DBOHO orig_bc=ACAGGTCG new_bc=ACAGGTCG bc_diffs=0", data="AGTGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGGAGATTAAGTAGCTTGCTATTTAATCTTAGTGGCGCACGGGTGAGTAATATATAGCTAATCTGCCCTACACTAGAGGACAACAGTTGGAAACGACTGCTAATACTCTATACTCCTTCTTTACATAAGTTAAGTCGGGAAAGTTTTTCGGTGTAGGATGAGGCTATATCGTATCAGCTAGTGGTAGGTAACGGCCTACCAAGGCTATGACGCGTAACTGGTCTGAGAGGATGATCAGTCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGGAATATTGCTCAAATGGGGGGAAAACCCTGAAAGCAGCAACGCCGCGTGGAGGATGACACTTTTCGGA"),
                           FASTARecord(id="PIDTA158_2", descr="HU82XDC01A3N0T orig_bc=ACCGCAGG new_bc=ACCGCAGG bc_diffs=0", data="GATGAACGCTAGCGATAGGCTTAACACATGCAAGTCGAGGGCATCACGAATTAGCAATAGTTTGGTGGCGACCGGCGCACGGGTGCGTAACACGTATACAACCTACCTTCAATTGGGGAATAACCTGGAGAAATTTGGACTAATACCCCATAGTAAACGGGAGAGGCATTCTTTTTTGTTTAAAGATTTATTGATTGGAGATGGGTATGCGTAGGATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAACGATCCTTAGGGGTT"),
                           FASTARecord(id="PIDF160_3", descr="HU82XDC01DTNIU orig_bc=ACCGTAGA new_bc=ACCGTAGA bc_diffs=0", data="GATGAACGCTGACAGAATGCTTAACACATGCAAGTCTACTTGAACTTCGGTTTGGGTGGCGGACGGGTGAGTAACGCGTAAAGAACTTGCCTCACAGTTAGGGACAACATTTGGAAACGAATGCTAATACCTGATATTATGATTTTAGGGCATCCTAAGATTATGAAAGCTATATGCGCTGTGAGAGAGCTTTGCGTCCCATTAGCTAGTTGGAGAGGTAACGGCTCACCAAGGCGATGATGGGTAGCCGGCCTGAGAGGGTGAACGGGCCACAAGGGGACTGAGACACGGCCCTTACTCCTACGGGAGGCAGCAGTGGGGAATATTGGGACAATGGAACCAAAAGTCTGATCCAGCAATTCTGTGTGCACGATG"),
                           FASTARecord(id="PIDTA.TB168_4", descr="HU82XDC01ETBU0 orig_bc=GCGCAACG new_bc=GCGCAACG bc_diffs=0", data="GATGAACGCTGACAGAATGCTTAACACATGCAAGTCAACTTGAACTTCGGTTTGGGTGGCGGACGGGTGAGTAACGCGTAAAGAACTTGCCTCACAGCTAGGGACAACATTTGGAAACGAATGCTAATACCTGATATTATGATTATATGGCATCGTATAATTATGAAAGCTATATGCGCTGTGAGAGAGCTTTGCGTCCCATTAGCTAGTTGGAGAGGTAACGGCTCACCAAGGCGATGATGGGTAGCCGGCCTGAGAGGGTGATCGGCCACAAGGGGACTGAGACACGGCCCTTACTCCTACGGGAGGCAGCAGTGGGGGAATATTGGGACAATGGGACCGAGAGTCTGATCCAGCAACTCTGTGTGCACGAT"),
                           FASTARecord(id="PIDTA.TB140_5", descr="HU82XDC01AVWB9 orig_bc=ACTGGAGA new_bc=ACTGGAGA bc_diffs=0", data="GATGAACGCTGACAGAATGCTTAACACATGCAAGTCAACTTGAATTTGGGTTTTAACTTAGATTTGGGTGGCGGACGGGTGAGTAACGCGTAAAGAACTTGCCTCACAGCTAGGGACAACATTTAGAAATGAATGCTAATACCTGATATTATGATTTTAAGGCATCTTAGAATTATGAAAGCTATAAGCACTGTGAGAGAGCTTTGCGTCCCATTAGCTAGTTGGAGAGGTAACAGCTCACCAAGGC")]
        for rec1, rec2 in zip(parseFASTA_result, manually_parsed):
            self.assertEqual(
                rec1, rec2,
                msg="FASTA records not parsed as expected."
            )

    def test_storeFASTA(self):
        """
        Testing storeFASTA function.

        :return: Returns OK if test goals were achieved, otherwise raises
                 error.
        """
        FASTARecord = namedtuple("FASTA_Record", "id descr data")
        storeFASTA_result = ut.storeFASTA("phylotoast/test/test_FASTA.fna")
        manually_parsed = [FASTARecord(id="PIDF154_1", descr="HU82XDC01DBOHO orig_bc=ACAGGTCG new_bc=ACAGGTCG bc_diffs=0", data="AGTGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGGAGATTAAGTAGCTTGCTATTTAATCTTAGTGGCGCACGGGTGAGTAATATATAGCTAATCTGCCCTACACTAGAGGACAACAGTTGGAAACGACTGCTAATACTCTATACTCCTTCTTTACATAAGTTAAGTCGGGAAAGTTTTTCGGTGTAGGATGAGGCTATATCGTATCAGCTAGTGGTAGGTAACGGCCTACCAAGGCTATGACGCGTAACTGGTCTGAGAGGATGATCAGTCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGGAATATTGCTCAAATGGGGGGAAAACCCTGAAAGCAGCAACGCCGCGTGGAGGATGACACTTTTCGGA"),
                           FASTARecord(id="PIDTA158_2", descr="HU82XDC01A3N0T orig_bc=ACCGCAGG new_bc=ACCGCAGG bc_diffs=0", data="GATGAACGCTAGCGATAGGCTTAACACATGCAAGTCGAGGGCATCACGAATTAGCAATAGTTTGGTGGCGACCGGCGCACGGGTGCGTAACACGTATACAACCTACCTTCAATTGGGGAATAACCTGGAGAAATTTGGACTAATACCCCATAGTAAACGGGAGAGGCATTCTTTTTTGTTTAAAGATTTATTGATTGGAGATGGGTATGCGTAGGATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAACGATCCTTAGGGGTT"),
                           FASTARecord(id="PIDF160_3", descr="HU82XDC01DTNIU orig_bc=ACCGTAGA new_bc=ACCGTAGA bc_diffs=0", data="GATGAACGCTGACAGAATGCTTAACACATGCAAGTCTACTTGAACTTCGGTTTGGGTGGCGGACGGGTGAGTAACGCGTAAAGAACTTGCCTCACAGTTAGGGACAACATTTGGAAACGAATGCTAATACCTGATATTATGATTTTAGGGCATCCTAAGATTATGAAAGCTATATGCGCTGTGAGAGAGCTTTGCGTCCCATTAGCTAGTTGGAGAGGTAACGGCTCACCAAGGCGATGATGGGTAGCCGGCCTGAGAGGGTGAACGGGCCACAAGGGGACTGAGACACGGCCCTTACTCCTACGGGAGGCAGCAGTGGGGAATATTGGGACAATGGAACCAAAAGTCTGATCCAGCAATTCTGTGTGCACGATG"),
                           FASTARecord(id="PIDTA.TB168_4", descr="HU82XDC01ETBU0 orig_bc=GCGCAACG new_bc=GCGCAACG bc_diffs=0", data="GATGAACGCTGACAGAATGCTTAACACATGCAAGTCAACTTGAACTTCGGTTTGGGTGGCGGACGGGTGAGTAACGCGTAAAGAACTTGCCTCACAGCTAGGGACAACATTTGGAAACGAATGCTAATACCTGATATTATGATTATATGGCATCGTATAATTATGAAAGCTATATGCGCTGTGAGAGAGCTTTGCGTCCCATTAGCTAGTTGGAGAGGTAACGGCTCACCAAGGCGATGATGGGTAGCCGGCCTGAGAGGGTGATCGGCCACAAGGGGACTGAGACACGGCCCTTACTCCTACGGGAGGCAGCAGTGGGGGAATATTGGGACAATGGGACCGAGAGTCTGATCCAGCAACTCTGTGTGCACGAT"),
                           FASTARecord(id="PIDTA.TB140_5", descr="HU82XDC01AVWB9 orig_bc=ACTGGAGA new_bc=ACTGGAGA bc_diffs=0", data="GATGAACGCTGACAGAATGCTTAACACATGCAAGTCAACTTGAATTTGGGTTTTAACTTAGATTTGGGTGGCGGACGGGTGAGTAACGCGTAAAGAACTTGCCTCACAGCTAGGGACAACATTTAGAAATGAATGCTAATACCTGATATTATGATTTTAAGGCATCTTAGAATTATGAAAGCTATAAGCACTGTGAGAGAGCTTTGCGTCCCATTAGCTAGTTGGAGAGGTAACAGCTCACCAAGGC")]
        for rec3, rec4 in zip(storeFASTA_result, manually_parsed):
            self.assertEqual(
                rec3, rec4,
                msg="FASTA records not parsed as expected."
            )

    def test_parse_map_file(self):
        """
        Testing parse_map_file function of util.py

        :return: Returns OK if test goals were achieved, otherwise raises
                 error.
        """
        # Obtaining lists to compare results
        header = ["#SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Treatment",
                  "Color", "Smoking", "Gender", "DOB", "Description"]
        data = {"PC.636": ["PC.636", "ACGGTGAGTGTC", "YATGCTGCCTCCCGTAGGAGT", "Fast",
                           "#0000CC", "Current_Smoker", "Female", "20080116",
                           "Fasting_I.D._636"],
                "PC.355": ["PC.355", "AACTCGTCGATG", "YATGCTGCCTCCCGTAGGAGT", "Control",
                           "#008000", "Current_Smoker", "Male", "20061218",
                           "Control_I.D._355"],
                "PC.607": ["PC.607", "AACTGTGCGTAC", "YATGCTGCCTCCCGTAGGAGT", "Fast",
                           "#0000CC", "Never_Smoker", "Male", "20071112",
                           "Fasting_I.D._607"],
                "PC.634": ["PC.634", "ACAGAGTCGGCT", "YATGCTGCCTCCCGTAGGAGT", "Fast",
                           "#0000CC", "Current_Smoker", "Female", "20080116",
                           "Fasting_I.D._634"],
                "PC.635": ["PC.635", "ACCGCAGAGTCA", "YATGCTGCCTCCCGTAGGAGT", "Fast",
                           "#0000CC", "Current_Smoker", "Male", "20080116",
                           "Fasting_I.D._635"],
                "PC.593": ["PC.593", "AGCAGCACTTGT", "YATGCTGCCTCCCGTAGGAGT", "Control",
                           "#008000", "Never_Smoker", "Female", "20071210",
                           "Control_I.D._593"],
                "PC.356": ["PC.356", "ACAGACCACTCA", "YATGCTGCCTCCCGTAGGAGT", "Control",
                           "#008000", "Current_Smoker", "Female", "20061126",
                           "Control_I.D._356"],
                "PC.481": ["PC.481", "ACCAGCGACTAG", "YATGCTGCCTCCCGTAGGAGT", "Control",
                           "#008000", "Never_Smoker", "Male", "20070314",
                           "Control_I.D._481"],
                "PC.354": ["PC.354", "AGCACGAGCCTA", "YATGCTGCCTCCCGTAGGAGT", "Control",
                           "#008000", "Never_Smoker", "Female", "20061218",
                           "Control_I.D._354"]}
        # Testing the validity
        self.assertEqual(
            self.map_header, header,
            msg="The mapping of sampleID's and the row information was not performed accurately."
        )
        self.assertEqual(
            self.map_data, data,
            msg="The parse_map_file() output was not keyed accurately on respective SampleID."
        )

    def test_write_map_file(self):
        """
        Testing write_map_file function of util.py.

        :return: Returns OK if test goals were achieved, otherwise raises
                 error.
        """
        mapfile = tempfile.NamedTemporaryFile(delete=False)
        ut.write_map_file(mapfile.name, self.map_data.values(), self.map_header)
        mapfile.seek(0)
        func_calc = mapfile.read().split("\n")[:-1]

        # Obtaining original file to compare with results.
        with open("phylotoast/test/test_mapping_file.txt") as inmapf:
            in_data = inmapf.read()
        hand_calc = in_data.split("\r")

        # Testing the validity of the function.
        for a, b in zip(hand_calc, func_calc):
            self.assertListEqual(
                a.split("\t"), b.split("\t"),
                msg="Mapping file data not written accurately to temp file."
                )

        # Deleting temporary files.
        os.unlink(mapfile.name)

    def test_parse_taxonomy_table(self):
        """
        Testing parse_taxonomy_table function.

        :return: Returns OK if test goals were achieved, otherwise raises
                 error.
        """
        taxa_data = ut.parse_taxonomy_table("phylotoast/test/test_taxa.txt")

        # Testing the validity of the function.
        hand_calc = {"018AP132": "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Neisseriaceae; g__Neisseria; s__HOT.018",
                     "057BE024": "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__HOT.057",
                     "083BS091": "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae_[XIVa]; g__Lachnoanaerobaculum; s__HOT.083",
                     "105_3039": "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Peptostreptococcaceae_[XI]; g__Eubacterium_[XI][G-1]; s__infirmum",
                     "122_8622": "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Megasphaera; s__micronuciformis",
                     "130Snoxi": "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Selenomonas; s__noxia",
                     "139EW076": "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Selenomonas; s__dianae",
                     "151_K168": "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Selenomonas; s__sputigena",
                     "214DE081": "k__Bacteria; p__Fusobacteria; c__Fusobacteria; o__Fusobacteriales; f__Leptotrichiaceae; g__Leptotrichia; s__shahii",
                     "220FB074": "k__Bacteria; p__Fusobacteria; c__Fusobacteria; o__Fusobacteriales; f__Leptotrichiaceae; g__Leptotrichiaceae_[G-1]; s__HOT.220",
                     "222_7816": "k__Bacteria; p__Fusobacteria; c__Fusobacteria; o__Fusobacteriales; f__Leptotrichiaceae; g__Leptotrichia; s__wadei"}
        for ids in hand_calc:
            self.assertEqual(
                taxa_data[ids], hand_calc[ids],
                msg="Taxonomy file was not accurately parsed into (OTU, taxonomy) dict."
                )

    def test_split_phylogeny(self):
        """
        Testing split_phylogeny() function of util.py.

        :return: Returns OK for successful run of the test, otherwise raises
                 error.
        """
        p1 = "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Veillonella; s__denticariosi"

        for lvl in ["k", "p", "c", "o", "f", "g", "s"]:
            if lvl == "k":
                self.assertEqual(
                    ut.split_phylogeny(p1, "k"), "k__Bacteria",
                    msg="Error. Identification failed at level 'k'."
                    )

            if lvl == "p":
                self.assertEqual(
                    ut.split_phylogeny(p1, "p"), "k__Bacteria; p__Firmicutes",
                    msg="Error. Identification failed at level 'p'."
                    )

            if lvl == "c":
                self.assertEqual(
                    ut.split_phylogeny(p1, "c"),
                    "k__Bacteria; p__Firmicutes; c__Clostridia",
                    msg="Error. Identification failed at level 'c'."
                    )

            if lvl == "o":
                self.assertEqual(
                    ut.split_phylogeny(p1, "o"),
                    "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales",
                    msg="Error. Identification failed at level 'o'."
                    )

            if lvl == "f":
                self.assertEqual(
                    ut.split_phylogeny(p1, "f"),
                    "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae",
                    msg="Error. Identification failed at level 'f'."
                    )

            if lvl == "g":
                self.assertEqual(
                    ut.split_phylogeny(p1, "g"),
                    "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Veillonella",
                    msg="Error. Identification failed at level 'g'."
                    )

            if lvl == "s":
                self.assertEqual(
                    ut.split_phylogeny(p1, "s"),
                    "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Veillonella; s__denticariosi",
                    msg="Error. Identification failed at level 's'."
                    )

    def test_ensure_dir(self):
        """
        Testing ensure_dir() function.

        :return: Returns OK if test goals were achieved, otherwise raises
                 error.
        """
        dir_path = "phylotoast/test/test_dir/"
        ut.ensure_dir(dir_path)
        self.assertTrue(os.path.isdir(dir_path),
                        msg="Directory could not be created.")
        os.rmdir(dir_path)
        self.assertFalse(os.path.isdir(dir_path), msg="test_dir could not be deleted.")

    def test_file_handle(self):
        """
        Testing file_handle() function of util.py.

        :return: Returns OK for successful run of the test, otherwise raises
                 error.
        """
        open_file_handle = ut.file_handle("phylotoast/test/test_mapping_file.txt")

        # Test function return value
        self.assertIsInstance(open_file_handle, file,
                              msg="Input did not return an open file handle.")
        with self.assertRaisesRegexp(ValueError, "Input file is closed."):
            f = open("phylotoast/test/test_mapping_file.txt", "rU")
            f.close()
            ut.file_handle(f)

    def test_gather_categories(self):
        """
        Testing gather_category function from iTol.py. If successful, the
        function will be moved to util.py.

        :return: Returns OK if test goals were achieved, otherwise raises
                error.
        """
        DataCategory = namedtuple("DataCategory", "sids results")
        result = ut.gather_categories(self.map_data, self.map_header)
        result1 = ut.gather_categories(self.map_data, self.map_header,
                                       ["Treatment"])  # one category given
        result2 = ut.gather_categories(self.map_data, self.map_header,
                                       ["Smoking=Control"])  # incorrect condition given
        result3 = ut.gather_categories(self.map_data, self.map_header,
                                       ["Treatment=Fast"])  # correct condition given
        result4 = ut.gather_categories(self.map_data, self.map_header,
                                       ["Treatment", "Smoking"])  # 2 categories given
        result5 = ut.gather_categories(self.map_data, self.map_header,
                                       ["Treatment", "Smoking=Never_Smoker"])  # 1 category 1 condition
        result6 = ut.gather_categories(self.map_data, self.map_header,
                                       ["Treatment=Control", "Smoking=Current_Smoker"])  # 1 category 1 condition
        result7 = ut.gather_categories(self.map_data, self.map_header,
                                       ["Smoking=Current_Smoker", "Smoking=Never_Smoker"])  # 2 conditions given
        result8 = ut.gather_categories(self.map_data, self.map_header,
                                       ["Smoking=Never_Smoker", "Treatment",
                                        "Gender=Female"])  # more than 2 categories - mix

        # Testing if the function calculates without any categories.
        manual = {"default": DataCategory({"PC.354", "PC.355", "PC.356", "PC.481",
                                           "PC.593", "PC.607", "PC.634", "PC.635",
                                           "PC.636"}, {})}
        self.assertDictEqual(
            result, manual,
            msg="With no category or condition given, gather_categories() did not return "
                "all SampleIDs as expected."
        )

        # Testing if the function accurately calculates for one category
        manual1 = {"Control": DataCategory({"PC.355", "PC.356", "PC.354",
                                            "PC.481", "PC.593"}, {}),
                   "Fast": DataCategory({"PC.634", "PC.635", "PC.636", "PC.607"}, {})}
        self.assertDictEqual(
            result1, manual1,
            msg="With one category given, gather_categories() did not return per "
                "category SampleIDs as expected."
        )

        # Testing if the function accurately calculates for incorrect condition given
        self.assertDictEqual(
            result2, manual,
            msg="With incorrect condition given, gather_categories() did not return "
                "all SampleIDs by default as expected."
        )

        # Testing if the function accurately calculates for correct one condition given
        manual3 = {"Fast": DataCategory({"PC.634", "PC.635", "PC.636", "PC.607"}, {})}
        self.assertDictEqual(
            result3, manual3,
            msg="With one correct condition given, gather_categories() did not return "
                "SampleIDs for the condition given, as expected."
        )

        # Testing if the function accurately calculates for correct one condition given
        manual4 = {"Control_Current_Smoker": DataCategory({"PC.355", "PC.356"}, {}),
                   "Control_Never_Smoker": DataCategory({"PC.354", "PC.481", "PC.593"}, {}),
                   "Fast_Current_Smoker": DataCategory({"PC.634", "PC.635", "PC.636"}, {}),
                   "Fast_Never_Smoker": DataCategory({"PC.607"}, {})}
        self.assertDictEqual(
            result4, manual4,
            msg="With multiple categories given, gather_categories() did not return "
                "SampleIDs for all category combinations, as expected."
        )

        # Testing if the function accurately calculates for one category and condition
        manual5 = {"Control_Never_Smoker": DataCategory({"PC.354", "PC.481", "PC.593"}, {}),
                   "Fast_Never_Smoker": DataCategory({"PC.607"}, {})}
        self.assertDictEqual(
            result5, manual5,
            msg="With one category and one condition given, gather_categories() did not "
                "return SampleIDs for all category-condition combinations, as expected."
        )

        # Testing if the function accurately calculates for one category and condition
        manual6 = {"Control_Current_Smoker": DataCategory({"PC.355", "PC.356"}, {})}
        self.assertDictEqual(
            result6, manual6,
            msg="With two specific conditions given, gather_categories() did not "
                "return SampleIDs for all condition combinations, as expected."
        )

        # Testing if the function accurately calculates for one category and condition
        manual7 = {"Current_Smoker": DataCategory({"PC.355", "PC.356", "PC.634",
                                                   "PC.635", "PC.636"}, {}),
                   "Never_Smoker": DataCategory({"PC.354", "PC.481", "PC.593", "PC.607"}, {})}
        self.assertDictEqual(
            result7, manual7,
            msg="With two conditions from same category given, gather_categories() did "
                "not return SampleIDs for all condition combinations, as expected."
        )

        # Testing if function accurately categorizes SampleIDs for multiple categories
        manual8 = {"Control_Never_Smoker_Female": DataCategory({"PC.354", "PC.593"}, {})}
        self.assertDictEqual(
            result8, manual8,
            msg="With two or more conditions/categories given, gather_categories() did "
                "not return SampleIDs for all condition combinations, as expected."
        )

    def test_parse_unifrac(self):
        """
        Testing parse_unifrac function of util.py.

        :return: Returns OK if test goals were achieved, otherwise raises
                 error.
        """
        with self.assertRaisesRegexp(ValueError, "File format not supported/recognized. "
                                     "Please check input unifrac file."):
            ut.parse_unifrac("phylotoast/test/test_taxa.txt")

        # Testing parse method for older unifrac format >= QIIME 1.9
        new_result = ut.parse_unifrac("phylotoast/test/test_new_unifrac.txt")
        new_pcd_keys = ["S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20"]
        self.assertEqual(
            new_result["pcd"].keys(), new_pcd_keys,
            msg="(QIIME 1.9) Principal Coordinate Vectors not calculated for all SampleIDs."
        )
        new_pcd = {"S11": [-0.157411912, 0.003113048, -0.033995362, 0.005761881, -0.08984718, 1.11E-06, 0.027054912, -0.020760899, 0.013421892, -0.025369248, 0, 0],
                   "S12": [-0.208758691, 0.119322024, 0.050022374, 0.057134613, -0.071335134, -0.002672489, 0.033457743, 0.015571627, 0.011477451, 0.03488542, 0, 0],
                   "S13": [0.068062494, 0.05191316, -0.053871698, -0.089607635, 0.040471386, 0.015394245, 0.001131232, 0.013414867, -0.051237119, -0.005849389, 0, 0],
                   "S14": [-0.051047792, 0.15554913, 0.051688267, -0.091735979, -0.001046725, 0.049856361, -0.039844435, 0.011459861, 0.006870138, 0.076553853, 0, 0],
                   "S15": [0.020092707, -0.136548397, -0.064225328, 0.001162091, 0.003647047, 0.036781766, 0.006616754, 0.028236267, 0.058819874, -0.007981049, 0, 0],
                   "S16": [-0.209895198, -0.172771574, -0.011870011, -0.021167774, 0.139398297, -0.04502447, 0.118228288, 0.079750916, 0.011682316, -0.040032913, 0, 0],
                   "S17": [-0.091537115, -0.13574047, -0.057444957, 0.019893426, 0.090763351, -0.005562631, 0.050826084, -0.043032858, 0.08985785, 0.006663828, 0, 0],
                   "S18": [-0.110229531, -0.432593524, -0.01330442, 0.105612832, 0.006608581, 0.070763799, -0.103434704, 0.039120951, -0.025255827, 0.038543059, 0, 0],
                   "S19": [-0.060261412, -0.03463619, -0.091585076, -0.027507828, -0.052794257, -0.003637563, 0.087998352, -0.08520105, -0.032339037, 0.041615824, 0, 0],
                   "S20": [0.25677973, -0.002601114, -0.006285835, 0.04308064, 0.020366331, 0.023860169, 0.021074805, 0.032318581, -0.037580605, -0.037847306, 0, 0]}
        new_eigenval = [1.036722355, 0.69841792, 0.214030904, 0.136343588, 0.131048126, 0.086219923, 0.076792575, 0.058771137, 0.04978709, 0.036961875, 0.036426793, 0]
        new_varexp = [38.7399495, 26.0982845, 7.9978466, 5.0948489, 4.8969696, 3.2218419, 2.8695634, 2.1961433, 1.86043, 1.3811809, 1.3611861, 0]
        for a1, b1 in zip(new_result["eigvals"], new_eigenval):
            self.assertAlmostEqual(
                a1, b1,
                msg="(QIIME 1.9) Eigenvalues not parsed accurately with their corresponding key."
            )
        for a2, b2 in zip(new_result["varexp"], new_varexp):
            self.assertAlmostEqual(
                a2, b2,
                msg="(QIIME 1.9) Variance values not parsed accurately to their corresponding key."
            )
        for sid in new_pcd:
            for a3, b3 in zip(new_result["pcd"][sid], new_pcd[sid]):
                self.assertAlmostEqual(
                    a3, b3,
                    msg="(QIIME 1.9) PCD points not parsed accurately to correspond to their keys."
                )

        # Testing parse method for older unifrac format <= QIIME 1.8
        old_result = ut.parse_unifrac("phylotoast/test/test_old_unifrac.txt")
        old_pcd_keys = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"]
        self.assertEqual(
            old_result["pcd"].keys(), old_pcd_keys,
            msg="(QIIME 1.8) Principal Coordinate Vectors not calculated for all SampleIDs."
        )
        old_pcd = {"S1": [0.100836831, 0.007753384, -0.057748508, 0.091966901, 0.020047889, 0.074339039, 0.043220932, -0.027168825, 0.066974054, 0.007435508],
                   "S2": [0.199962746, 0.073960653, -0.032422978, -0.042610315, -0.022374308, 0.003714962, -0.051645417, 0.056297627, -0.042463295, 0.043060609],
                   "S3": [0.047565709, 0.168165361, -0.074130952, 0.02487798, -0.079852925, 0.000323633, -0.003968622, 0.008167576, 0.004530349, -0.056961177],
                   "S4": [0.140443295, -0.039475886, -0.025525758, 0.042153004, 0.110942957, 0.012739223, 0.033180866, 0.011832021, 0.068494426, -0.047639907],
                   "S5": [-0.136757116, -0.120031189, 0.017439573, 0.039528692, 0.039809067, -0.034443631, 0.093800166, 0.212262866, -0.01425893, -0.044323662],
                   "S6": [-0.248686792, 0.082467735, -0.040308376, 0.118185703, 0.135928483, 0.043041123, -0.039510796, -0.035302159, -0.055888817, 0.013733698],
                   "S7": [0.013351414, 0.090151459, 0.056540284, 0.025991866, -0.029388885, -0.032180352, -0.094075979, -0.013512448, 0.070055931, -0.049158479],
                   "S8": [0.030827257, 0.185374719, 0.066041818, -0.091731703, 0.028183239, 0.016665481, -0.056988868, 0.011934123, 0.069909254, -0.037223541],
                   "S9": [0.089360294, 0.091572201, 0.093558961, -0.105185682, 0.006075005, -0.01497293, -0.080077528, 0.038141749, 0.036383198, -0.027530755],
                   "S10": [0.134964652, 0.019636571, -0.0054345, 0.060611834, 0.000711488, 0.063693172, 0.074841712, -0.013591079, 0.068696143, 0.003720656]}
        old_eigenval = [0.938839916, 0.376712138, 0.229056828, 0.196605771, 0.168680373, 0.15734919, 0.14605245, 0.127838108, 0.118759911, 0.093651517]
        old_varexp = [26.84931482, 10.77336255, 6.550657673, 5.622609514, 4.823987944, 4.499934304, 4.176865661, 3.655964721, 3.396342854, 2.678283087]
        for a4, b4 in zip(old_result["eigvals"], old_eigenval):
            self.assertAlmostEqual(
                a4, b4,
                msg="(QIIME 1.8) Eigenvalues not parsed accurately with their corresponding key."
            )
        for a5, b5 in zip(old_result["varexp"], old_varexp):
            self.assertAlmostEqual(
                a5, b5,
                msg="(QIIME 1.8) Variance values not parsed accurately to their corresponding key."
            )
        for sid in old_pcd:
            for a6, b6 in zip(old_result["pcd"][sid], old_pcd[sid]):
                self.assertAlmostEqual(
                    a6, b6,
                    msg="(QIIME 1.8) PCD points not parsed accurately to correspond to their keys."
                )

    def test_color_mapping(self):
        """
        Testing the color-group mapping for obtaining colors for visualizations
        from mapping file.

        :return:Returns OK if test goals were achieved, otherwise raises
                 error.
        """
        colormap1 = ut.color_mapping(self.map_data, self.map_header, "Treatment", "Color")
        self.assertEqual({"Control": "#008000", "Fast": "#0000CC"},
                         colormap1, msg="Color-group mapping not computed "
                         "accurately. Please check category and color "
                         "columns.")
        colormap2 = ut.color_mapping(self.map_data, self.map_header, "Treatment")
        self.assertEqual({"Control": "#8DD3C7", "Fast": "#FFFFB3"},
                         colormap2, msg="With no color column given, the "
                         "color-group mapping not computed accurately.")

    def tearDown(self):
        """
        Closing all files or directories that are created or open for testing.
        """
        pass

if __name__ == "__main__":
    unittest.main()
