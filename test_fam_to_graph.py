import unittest
import json
import sys
import os
import tempfile
import subprocess
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from fam_to_graph import GraphMaker


class TestPanacondaBaseline(unittest.TestCase):
    def setUp(self):
        # Create a mock patric feature tab file
        self.mock_features = """#genome\torganism\tcontig\tfeature\tstart\tend\tgroup
genome_1\tMockus_bacterium\tcontig_1\tfig|1.1\t100\t200\tgroupA
genome_1\tMockus_bacterium\tcontig_1\tfig|1.2\t250\t350\tgroupB
genome_1\tMockus_bacterium\tcontig_1\tfig|1.3\t400\t500\tgroupC
genome_2\tMockus_bacterium\tcontig_1\tfig|2.1\t110\t210\tgroupA
genome_2\tMockus_bacterium\tcontig_1\tfig|2.2\t260\t360\tgroupB
genome_2\tMockus_bacterium\tcontig_1\tfig|2.3\t410\t510\tgroupC"""

        self.temp_dir = tempfile.TemporaryDirectory()
        self.feature_file = os.path.join(self.temp_dir.name, "mock_features.txt")

        with open(self.feature_file, "w") as f:
            f.write(self.mock_features)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_graph_creation_and_attributes(self):
        """Test that the graph generates properly with required frontend attributes."""
        gmaker = GraphMaker(
            feature_files=[self.feature_file],
            feature_selector="genome_id",
            file_type="patric_tab",
            context="genome",
            ksize=3,
            break_conflict=False,
            label_function=False,
            diversity="genus",
            minSeq=1,
            alpha="pgfam_id"
        )

        gmaker.processFeatures()
        gmaker.RF_to_PG()
        gmaker.checkPGGraph()
        gmaker.checkRFGraph()
        gmaker.calcStatistics()
        gmaker.compute_graph_metrics()
        gmaker.detect_and_annotate_superbubbles()
        gmaker.tag_structural_variants()
        gmaker.flag_bridges()
        gmaker.annotate_major_blocks_tfs()

        undirected = gmaker.get_undirected_pg_graph()
        gmaker.serialize_graph_for_gexf(undirected)

        # 1. Verify Nodes and Edges exist
        self.assertGreater(undirected.number_of_nodes(), 0, "Graph should have nodes")

        # 2. Check for frontend-required Node attributes
        for node, data in undirected.nodes(data=True):
            self.assertIn("features", data, "Node missing 'features' attribute")
            self.assertIn("diversity", data, "Node missing 'diversity' attribute")

            # 3. Verify JSON Formatting is valid
            parsed_features = json.loads(data["features"])
            self.assertIsInstance(parsed_features, dict, "Features JSON should parse to a dict")

        # 4. Check for frontend-required Edge attributes
        for u, v, data in undirected.edges(data=True):
            self.assertIn("genomes", data, "Edge missing 'genomes' attribute")
            self.assertIn("sequences", data, "Edge missing 'sequences' attribute")

    def test_max_instance_overlap_fix(self):
        """Verify maxInstanceOverlap updates max_intersect and returns only maximum overlap keys."""
        gmaker = GraphMaker(
            feature_files=[self.feature_file],
            feature_selector="genome_id",
            file_type="patric_tab",
            context="genome",
            ksize=3
        )
        res = gmaker.maxInstanceOverlap("1.2.3.1r", ["1.9.9.1r", "1.2.3.1r", "1.2.8.1r"])
        self.assertEqual(res, ["1.2.3.1r"])

    def test_parallel_paralog_guidance_and_reports(self):
        """Verify parallel chromosomes/contigs entering a shared repeat stay guided to separate paralog paths."""
        # Two genomes, each with two chromosomes sharing a 4-gene repeat (R1-R2-R3-R4) flanked by distinct anchors
        features = """#genome\torganism\tcontig\tfeature\tstart\tend\tgroup
g1\tOrg1\tchr1\tf1_1\t100\t200\tA1
g1\tOrg1\tchr1\tf1_2\t200\t300\tA2
g1\tOrg1\tchr1\tf1_3\t300\t400\tA3
g1\tOrg1\tchr1\tf1_4\t400\t500\tR1
g1\tOrg1\tchr1\tf1_5\t500\t600\tR2
g1\tOrg1\tchr1\tf1_6\t600\t700\tR3
g1\tOrg1\tchr1\tf1_7\t700\t800\tR4
g1\tOrg1\tchr1\tf1_8\t800\t900\tB1
g1\tOrg1\tchr1\tf1_9\t900\t1000\tB2
g1\tOrg1\tchr1\tf1_10\t1000\t1100\tB3
g1\tOrg1\tchr2\tf2_1\t100\t200\tC1
g1\tOrg1\tchr2\tf2_2\t200\t300\tC2
g1\tOrg1\tchr2\tf2_3\t300\t400\tC3
g1\tOrg1\tchr2\tf2_4\t400\t500\tR1
g1\tOrg1\tchr2\tf2_5\t500\t600\tR2
g1\tOrg1\tchr2\tf2_6\t600\t700\tR3
g1\tOrg1\tchr2\tf2_7\t700\t800\tR4
g1\tOrg1\tchr2\tf2_8\t800\t900\tD1
g1\tOrg1\tchr2\tf2_9\t900\t1000\tD2
g1\tOrg1\tchr2\tf2_10\t1000\t1100\tD3
g2\tOrg2\tchr1\tf3_1\t100\t200\tA1
g2\tOrg2\tchr1\tf3_2\t200\t300\tA2
g2\tOrg2\tchr1\tf3_3\t300\t400\tA3
g2\tOrg2\tchr1\tf3_4\t400\t500\tR1
g2\tOrg2\tchr1\tf3_5\t500\t600\tR2
g2\tOrg2\tchr1\tf3_6\t600\t700\tR3
g2\tOrg2\tchr1\tf3_7\t700\t800\tR4
g2\tOrg2\tchr1\tf3_8\t800\t900\tB1
g2\tOrg2\tchr1\tf3_9\t900\t1000\tB2
g2\tOrg2\tchr1\tf3_10\t1000\t1100\tB3
g2\tOrg2\tchr2\tf4_1\t100\t200\tC1
g2\tOrg2\tchr2\tf4_2\t200\t300\tC2
g2\tOrg2\tchr2\tf4_3\t300\t400\tC3
g2\tOrg2\tchr2\tf4_4\t400\t500\tR1
g2\tOrg2\tchr2\tf4_5\t500\t600\tR2
g2\tOrg2\tchr2\tf4_6\t600\t700\tR3
g2\tOrg2\tchr2\tf4_7\t700\t800\tR4
g2\tOrg2\tchr2\tf4_8\t800\t900\tD1
g2\tOrg2\tchr2\tf4_9\t900\t1000\tD2
g2\tOrg2\tchr2\tf4_10\t1000\t1100\tD3"""
        repeat_file = os.path.join(self.temp_dir.name, "repeat_features.txt")
        out_gexf = os.path.join(self.temp_dir.name, "repeat_out.gexf")
        with open(repeat_file, "w") as f:
            f.write(features)

        cmd = [sys.executable, os.path.join(os.path.dirname(__file__), "src", "fam_to_graph.py"), "--output", out_gexf, repeat_file]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"CLI failed: {res.stderr}")
        self.assertTrue(os.path.exists(out_gexf))
        for report_suffix in ["_synteny_similarity.tsv", "_mobilized_gene_families.tsv", "_superbubbles.tsv", "_structural_variants.tsv", "_scaffold_recommendations.tsv"]:
            self.assertTrue(os.path.exists(os.path.join(self.temp_dir.name, "repeat_out" + report_suffix)))

        # Verify that the two parallel chromosomes did not overwrite each other's repeat assignments:
        # 6 unique nodes on chr1 (A1-A3, B1-B3) + 6 unique nodes on chr2 (C1-C3, D1-D3) + 4*2=8 repeat nodes = 20 nodes
        tree = ET.parse(out_gexf)
        root = tree.getroot()
        summary_el = root.find(".//{*}summary")
        self.assertIsNotNone(summary_el)
        summary = json.loads(summary_el.text)
        self.assertEqual(summary["total_nodes"], 20)
        self.assertEqual(summary["copy_number_variants"], 4)


if __name__ == '__main__':
    unittest.main()

