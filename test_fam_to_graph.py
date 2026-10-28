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

    def test_inversion_and_shift_repeat_ambiguity(self):
        """Verify opposite-direction contigs restore reference strand, true DNA inversions stay contiguous across k-mers, and shift splits preserve repeat_ambiguity."""
        # g1 & g2 forward: A1(+) A2(+) A3(+) B1(+) B2(+) B3(+) B4(+) C1(+) C2(+) C3(+)
        # g3 has B1..B4 inverted in both gene order and DNA strand: A1(+) A2(+) A3(+) B4(-) B3(-) B2(-) B1(-) C1(+) C2(+) C3(+)
        # g4 is the entire contig reverse-complemented (opposite assembly direction, NOT a rearrangement relative to g1/g2):
        #   C3(-) C2(-) C1(-) B4(-) B3(-) B2(-) B1(-) A3(-) A2(-) A1(-)
        features = """#genome\torganism\tcontig\tfeature\tstart\tend\tstrand\tgroup
g1\tOrg1\tc1\tf1_1\t100\t200\t+\tA1
g1\tOrg1\tc1\tf1_2\t200\t300\t+\tA2
g1\tOrg1\tc1\tf1_3\t300\t400\t+\tA3
g1\tOrg1\tc1\tf1_4\t400\t500\t+\tB1
g1\tOrg1\tc1\tf1_5\t500\t600\t+\tB2
g1\tOrg1\tc1\tf1_6\t600\t700\t+\tB3
g1\tOrg1\tc1\tf1_7\t700\t800\t+\tB4
g1\tOrg1\tc1\tf1_8\t800\t900\t+\tC1
g1\tOrg1\tc1\tf1_9\t900\t1000\t+\tC2
g1\tOrg1\tc1\tf1_10\t1000\t1100\t+\tC3
g2\tOrg2\tc1\tf2_1\t100\t200\t+\tA1
g2\tOrg2\tc1\tf2_2\t200\t300\t+\tA2
g2\tOrg2\tc1\tf2_3\t300\t400\t+\tA3
g2\tOrg2\tc1\tf2_4\t400\t500\t+\tB1
g2\tOrg2\tc1\tf2_5\t500\t600\t+\tB2
g2\tOrg2\tc1\tf2_6\t600\t700\t+\tB3
g2\tOrg2\tc1\tf2_7\t700\t800\t+\tB4
g2\tOrg2\tc1\tf2_8\t800\t900\t+\tC1
g2\tOrg2\tc1\tf2_9\t900\t1000\t+\tC2
g2\tOrg2\tc1\tf2_10\t1000\t1100\t+\tC3
g3\tOrg3\tc1\tf3_1\t100\t200\t+\tA1
g3\tOrg3\tc1\tf3_2\t200\t300\t+\tA2
g3\tOrg3\tc1\tf3_3\t300\t400\t+\tA3
g3\tOrg3\tc1\tf3_4\t400\t500\t-\tB4
g3\tOrg3\tc1\tf3_5\t500\t600\t-\tB3
g3\tOrg3\tc1\tf3_6\t600\t700\t-\tB2
g3\tOrg3\tc1\tf3_7\t700\t800\t-\tB1
g3\tOrg3\tc1\tf3_8\t800\t900\t+\tC1
g3\tOrg3\tc1\tf3_9\t900\t1000\t+\tC2
g3\tOrg3\tc1\tf3_10\t1000\t1100\t+\tC3
g4\tOrg4\tc1\tf4_1\t100\t200\t-\tC3
g4\tOrg4\tc1\tf4_2\t200\t300\t-\tC2
g4\tOrg4\tc1\tf4_3\t300\t400\t-\tC1
g4\tOrg4\tc1\tf4_4\t400\t500\t-\tB4
g4\tOrg4\tc1\tf4_5\t500\t600\t-\tB3
g4\tOrg4\tc1\tf4_6\t600\t700\t-\tB2
g4\tOrg4\tc1\tf4_7\t700\t800\t-\tB1
g4\tOrg4\tc1\tf4_8\t800\t900\t-\tA3
g4\tOrg4\tc1\tf4_9\t900\t1000\t-\tA2
g4\tOrg4\tc1\tf4_10\t1000\t1100\t-\tA1"""
        inv_file = os.path.join(self.temp_dir.name, "inv_features.txt")
        with open(inv_file, "w") as f:
            f.write(features)

        gmaker = GraphMaker(
            feature_files=[inv_file],
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
        gmaker.calcStatistics()
        gmaker.compute_graph_metrics()
        inv_events, _ = gmaker.tag_structural_variants()

        # Exactly 1 contiguous inverted component (B1-B2-B3-B4 plus boundary edges), and A1-A2-A3 / C1-C2-C3 are NOT inverted
        self.assertEqual(inv_events, 1)
        fam_to_node = {d["label"]: n for n, d in gmaker.pg_graph.nodes(data=True)}
        for non_inv in ["A1", "A2", "A3", "C1", "C2", "C3"]:
            nid = fam_to_node[non_inv]
            self.assertFalse(gmaker.pg_graph.nodes[nid].get("inverted_block", False), f"{non_inv} should not be marked inverted_block")
        for inv_fam in ["B1", "B2", "B3", "B4"]:
            nid = fam_to_node[inv_fam]
            self.assertTrue(gmaker.pg_graph.nodes[nid].get("inverted_block", False), f"{inv_fam} should be marked inverted_block")
        # Add a localized repeat R1-R2-R3 duplicated in g2 so shift conflict splits R1-R2-R3 into original + new_pg nodes,
        # and verify BOTH original and split nodes have repeat_ambiguity == True.
        shift_features = """#genome\torganism\tcontig\tfeature\tstart\tend\tstrand\tgroup
g1\tOrg1\tc1\ts1_1\t100\t200\t+\tA1
g1\tOrg1\tc1\ts1_2\t200\t300\t+\tA2
g1\tOrg1\tc1\ts1_3\t300\t400\t+\tA3
g1\tOrg1\tc1\ts1_4\t400\t500\t+\tR1
g1\tOrg1\tc1\ts1_5\t500\t600\t+\tR2
g1\tOrg1\tc1\ts1_6\t600\t700\t+\tR3
g1\tOrg1\tc1\ts1_7\t700\t800\t+\tC1
g1\tOrg1\tc1\ts1_8\t800\t900\t+\tC2
g1\tOrg1\tc1\ts1_9\t900\t1000\t+\tC3
g2\tOrg2\tc1\ts2_1\t100\t200\t+\tA1
g2\tOrg2\tc1\ts2_2\t200\t300\t+\tA2
g2\tOrg2\tc1\ts2_3\t300\t400\t+\tA3
g2\tOrg2\tc1\ts2_4\t400\t500\t+\tR1
g2\tOrg2\tc1\ts2_5\t500\t600\t+\tR2
g2\tOrg2\tc1\ts2_6\t600\t700\t+\tR3
g2\tOrg2\tc1\ts2_7\t700\t800\t+\tR1
g2\tOrg2\tc1\ts2_8\t800\t900\t+\tR2
g2\tOrg2\tc1\ts2_9\t900\t1000\t+\tR3
g2\tOrg2\tc1\ts2_10\t1000\t1100\t+\tC1
g2\tOrg2\tc1\ts2_11\t1100\t1200\t+\tC2
g2\tOrg2\tc1\ts2_12\t1200\t1300\t+\tC3"""
        shift_file = os.path.join(self.temp_dir.name, "shift_features.txt")
        with open(shift_file, "w") as f:
            f.write(shift_features)
        gm_shift = GraphMaker(
            feature_files=[shift_file],
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
        gm_shift.processFeatures()
        gm_shift.RF_to_PG()
        gm_shift.calcStatistics()
        gm_shift.compute_graph_metrics()
        r_nodes = [n for n, d in gm_shift.pg_graph.nodes(data=True) if d.get("label") in {"R1", "R2", "R3"}]
        self.assertEqual(len(r_nodes), 6, "Expected 3 original + 3 split PG nodes for R1..R3")
        for rn in r_nodes:
            self.assertTrue(gm_shift.pg_graph.nodes[rn].get("repeat_ambiguity", False), f"Node {rn} ({gm_shift.pg_graph.nodes[rn].get('label')}) missing repeat_ambiguity")
            self.assertGreater(gm_shift.pg_graph.nodes[rn].get("cnv_cluster_id", 0), 0)

        # Verify isolated paralogs in different k-mer contexts on different contigs do NOT get flagged as repeat_ambiguity or repeat_ambiguity_fragmentation
        iso_paralog_features = """#genome\torganism\tcontig\tfeature\tstart\tend\tstrand\tgroup
g1\tOrg1\tchr1\tp1_1\t100\t200\t+\tA1
g1\tOrg1\tchr1\tp1_2\t200\t300\t+\tA2
g1\tOrg1\tchr1\tp1_3\t300\t400\t+\tA3
g1\tOrg1\tchr1\tp1_4\t400\t500\t+\tPARA
g1\tOrg1\tchr1\tp1_5\t500\t600\t+\tB1
g1\tOrg1\tchr1\tp1_6\t600\t700\t+\tB2
g1\tOrg1\tchr1\tp1_7\t700\t800\t+\tB3
g1\tOrg1\tchr2\tp2_1\t100\t200\t+\tC1
g1\tOrg1\tchr2\tp2_2\t200\t300\t+\tC2
g1\tOrg1\tchr2\tp2_3\t300\t400\t+\tC3
g1\tOrg1\tchr2\tp2_4\t400\t500\t+\tPARA
g1\tOrg1\tchr2\tp2_5\t500\t600\t+\tD1
g1\tOrg1\tchr2\tp2_6\t600\t700\t+\tD2
g1\tOrg1\tchr2\tp2_7\t700\t800\t+\tD3
g2\tOrg2\tchr1\tp3_1\t100\t200\t+\tA1
g2\tOrg2\tchr1\tp3_2\t200\t300\t+\tA2
g2\tOrg2\tchr1\tp3_3\t300\t400\t+\tA3
g2\tOrg2\tchr1\tp3_4\t400\t500\t+\tPARA
g2\tOrg2\tchr1\tp3_5\t500\t600\t+\tB1
g2\tOrg2\tchr1\tp3_6\t600\t700\t+\tB2
g2\tOrg2\tchr1\tp3_7\t700\t800\t+\tB3
g2\tOrg2\tchr2\tp4_1\t100\t200\t+\tC1
g2\tOrg2\tchr2\tp4_2\t200\t300\t+\tC2
g2\tOrg2\tchr2\tp4_3\t300\t400\t+\tC3
g2\tOrg2\tchr2\tp4_4\t400\t500\t+\tPARA
g2\tOrg2\tchr2\tp4_5\t500\t600\t+\tD1
g2\tOrg2\tchr2\tp4_6\t600\t700\t+\tD2
g2\tOrg2\tchr2\tp4_7\t700\t800\t+\tD3"""
        iso_file = os.path.join(self.temp_dir.name, "iso_paralog.txt")
        with open(iso_file, "w") as f:
            f.write(iso_paralog_features)
        gm_iso = GraphMaker(
            feature_files=[iso_file],
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
        gm_iso.processFeatures()
        gm_iso.RF_to_PG()
        gm_iso.calcStatistics()
        gm_iso.compute_graph_metrics()
        para_nodes = [n for n, d in gm_iso.pg_graph.nodes(data=True) if d.get("family") == "PARA"]
        self.assertEqual(len(para_nodes), 2, "Expected 2 distinct PG-nodes for isolated paralog PARA in different k-mer contexts")
        for pn in para_nodes:
            self.assertEqual(gm_iso.pg_graph.nodes[pn].get("cnv_cluster_id", 0), 0)
            self.assertFalse(gm_iso.pg_graph.nodes[pn].get("repeat_ambiguity", False))
            self.assertFalse(gm_iso.pg_graph.nodes[pn].get("repeat_ambiguity_fragmentation", False))


if __name__ == '__main__':
    unittest.main()

