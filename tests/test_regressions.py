"""Focused regression tests for the reviewed LiftOverZ fixes.

Run from the repository root with:
    python -m unittest discover -s tests -v

A tiny pyliftover stub is used so these safety/format tests can also run in
minimal environments without the external conversion library installed.
"""
from __future__ import annotations

import gzip
import importlib.util
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Import the application's pure-Python validation and format logic without
# requiring a real chain reader. Individual tests supply fake coordinate maps.
if "pyliftover" not in sys.modules:
    stub = types.ModuleType("pyliftover")
    class UnusedLiftOver:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("The pyliftover stub is not a real chain reader")
    stub.LiftOver = UnusedLiftOver
    sys.modules["pyliftover"] = stub

spec = importlib.util.spec_from_file_location("liftover_gui_tested", ROOT / "liftover_gui.py")
app = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = app
assert spec.loader is not None
spec.loader.exec_module(app)

patch_spec = importlib.util.spec_from_file_location("patch_chains_tested", ROOT / "patch_chains.py")
patcher = importlib.util.module_from_spec(patch_spec)
sys.modules[patch_spec.name] = patcher
assert patch_spec.loader is not None
patch_spec.loader.exec_module(patcher)


class FakeLiftOver:
    def __init__(self, mapper):
        self.mapper = mapper

    def convert_coordinate(self, chrom, pos):
        return self.mapper(chrom, pos)


class RegressionTests(unittest.TestCase):
    def test_unmapped_sidecar_cannot_overwrite_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "result.vcf"
            input_path = Path(str(output) + ".unmapped")
            input_path.write_text("important input\n", encoding="utf-8")
            err = app.check_paths(str(input_path), str(output), cmd_type="vcf", save_unmapped=True)
            self.assertEqual(err, "err_output_conflict")
            self.assertEqual(input_path.read_text(encoding="utf-8"), "important input\n")

    def test_output_cannot_overwrite_chain_or_map(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "input.bed"
            chain_path = Path(tmp) / "chain.gz"
            input_path.write_text("chr1\t0\t1\n", encoding="utf-8")
            chain_path.write_bytes(b"chain placeholder")
            self.assertEqual(
                app.check_paths(str(input_path), str(chain_path), protected_paths=[str(chain_path)]),
                "err_output_conflict",
            )

    def test_wig_interval_with_changed_span_is_unmapped(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.wig"
            output = Path(tmp) / "output.wig"
            source.write_text("variableStep chrom=chrM span=5\n1 2.5\n", encoding="utf-8")

            def mapper(chrom, pos):
                # Source bases 0..4 map to target bases 10..13: length changes 5 -> 4.
                return [("chr1", 10 + pos if pos < 3 else 13, "+")]

            logs = []
            engine = app.LiftOverEngine(
                FakeLiftOver(mapper), None, {}, True, True, logs.append, lambda key, **kw: key,
            )
            stats = engine.run("wig", str(source), str(output))
            self.assertEqual(stats["mapped"], 0)
            self.assertEqual(stats["unmapped"], 1)
            self.assertTrue(Path(str(output) + ".unmapped").exists())
            self.assertNotIn("2.5", output.read_text(encoding="utf-8"))

    def test_stale_assembly_metadata_is_removed_from_vcf(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.vcf"
            output = Path(tmp) / "output.vcf"
            source.write_text(
                "##fileformat=VCFv4.2\n"
                "##assembly=hg19\n"
                "##reference=old_reference.fa\n"
                "##contig=<ID=chr1,length=100>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t1\t.\tA\tC\t.\tPASS\t.\n",
                encoding="utf-8",
            )
            def mapper(chrom, pos):
                return [(chrom, pos, "+")]
            engine = app.LiftOverEngine(
                FakeLiftOver(mapper), None, {}, True, False, lambda msg: None,
                lambda key, **kw: key,
            )
            engine.run("vcf", str(source), str(output))
            text = output.read_text(encoding="utf-8")
            self.assertNotIn("##assembly=hg19", text)
            self.assertNotIn("##reference=old_reference.fa", text)
            self.assertNotIn("##contig=<ID=chr1,length=100>", text)
            self.assertIn("##contig=<ID=chr1>", text)

    def test_chain_patch_uses_exact_source_contig_name_not_mito_alias_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            chain_path = Path(tmp) / "test.chain.gz"
            text = (
                # Source chrM chain to a nuclear contig (NUMT): this is replaced.
                "chain 100 chrM 16569 + 0 10 chr1 248956422 + 10 20 1\n"
                "10\n\n"
                # chrMT is a separate named source contig, even with a similar size.
                "chain 100 chrMT 16569 + 0 10 chrMT 16569 + 0 10 2\n"
                "10\n\n"
                "chain 100 chr1 248956422 + 0 10 chr1 248956422 + 0 10 3\n"
                "10\n\n"
            )
            with gzip.open(chain_path, "wt", encoding="utf-8") as fh:
                fh.write(text)
            removed = patcher.patch_chain_file(str(chain_path), patcher.PATCHES["hg38ToHg19.over.chain.gz"])
            self.assertEqual(removed, 1)
            with gzip.open(chain_path, "rt", encoding="utf-8") as fh:
                result = fh.read()
            self.assertNotIn("chain 100 chrM 16569 + 0 10 chr1", result)
            self.assertIn("chain 100 chrMT 16569 + 0 10 chrMT", result)
            self.assertIn("chain 100 chr1 248956422", result)
            self.assertIn("chain 100000000 chrM 16569", result)
            patcher._verify_gzip(str(chain_path))

    def test_chrmt_source_does_not_fall_back_to_legacy_chrm(self):
        class Chain:
            source_name = "chrM"
            source_size = 16571
        class ChainFile:
            chains = [Chain()]
        class NamedLiftOver(FakeLiftOver):
            chain_file = ChainFile()

        lo = NamedLiftOver(lambda chrom, pos: [("chrM", pos, "+")] if chrom == "chrM" else None)
        engine = app.LiftOverEngine(lo, None, {}, True, False, lambda msg: None,
                                    lambda key, **kw: key)
        # chrMT represents the current 16,569-bp contig; don't route it to
        # the explicitly known legacy 16,571-bp chrM just because auto_chr is on.
        self.assertIsNone(engine._lift_coord("chrMT", 100))
        self.assertEqual(engine._lift_coord("chrM", 100), ("chrM", 100, "+"))

    def test_mitochondrial_output_keeps_actual_target_contig_name(self):
        self.assertEqual(app.LiftOverEngine._format_output_chrom("chrM", "chrMT"), "chrMT")
        self.assertEqual(app.LiftOverEngine._format_output_chrom("MT", "chrM"), "chrM")


    def test_hg19_chrmt_gets_explicit_rcrs_identity_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            chain_path = Path(tmp) / "hg19ToHg38.over.chain.gz"
            text = (
                # Legacy hg19 chrM is patched separately (length 16,571).
                "chain 100 chrM 16571 + 0 10 chr1 248956422 + 10 20 10\n"
                "10\n\n"
                # Existing chrMT NUMT/fragmented alignment is replaced by the
                # explicit, known rCRS-to-rCRS identity block.
                "chain 100 chrMT 16569 + 0 10 chr2 242193529 + 10 20 11\n"
                "10\n\n"
                "chain 100 chr1 248956422 + 0 10 chr1 248956422 + 0 10 12\n"
                "10\n\n"
            )
            with gzip.open(chain_path, "wt", encoding="utf-8") as fh:
                fh.write(text)
            removed_legacy = patcher.patch_chain_file(
                str(chain_path), patcher.PATCHES["hg19ToHg38.over.chain.gz"]
            )
            self.assertEqual(removed_legacy, 1)
            identity_block = patcher.EXTRA_PATCHES["hg19ToHg38.over.chain.gz"]
            removed_chrmt = patcher.patch_chain_file(
                str(chain_path), identity_block, require_match=False
            )
            self.assertEqual(removed_chrmt, 1)
            with gzip.open(chain_path, "rt", encoding="utf-8") as fh:
                result = fh.read()
            self.assertIn("chain 100000000 chrM 16571", result)
            self.assertIn("chain 100000001 chrMT 16569 + 0 16569 chrM 16569", result)
            self.assertNotIn("chain 100 chrMT 16569 + 0 10 chr2", result)
            patcher._verify_gzip(str(chain_path))

    def test_unprefixed_mt_alias_prefers_chrmt_when_both_source_contigs_exist(self):
        class LegacyChain:
            source_name = "chrM"
            source_size = 16571
        class CurrentChain:
            source_name = "chrMT"
            source_size = 16569
        class ChainFile:
            chains = [LegacyChain(), CurrentChain()]
        class NamedLiftOver(FakeLiftOver):
            chain_file = ChainFile()

        calls = []
        def mapper(chrom, pos):
            calls.append(chrom)
            if chrom == "chrMT":
                return [("chrM", pos, "+")]
            if chrom == "chrM":
                return [("chrM_legacy", pos, "+")]
            return None

        engine = app.LiftOverEngine(
            NamedLiftOver(mapper), None, {}, True, False, lambda msg: None,
            lambda key, **kw: key,
        )
        self.assertEqual(engine._lift_coord("MT", 12), ("chrM", 12, "+"))
        self.assertIn("chrMT", calls)
        self.assertNotIn("chrM", calls)

    def test_fasta_lookup_does_not_use_legacy_chrm_for_chrmt(self):
        fasta = {"chrM": "A" * 16571}
        self.assertIsNone(app.get_fasta_seq(fasta, "chrMT", 0, 1))
        self.assertEqual(app.get_fasta_seq(fasta, "chrM", 0, 1), "A")
        self.assertIsNone(app.get_fasta_seq({"chrMT": "C" * 16569}, "chrM", 0, 1))
        self.assertIsNone(app.lookup_length({"chrM": 16571}, "chrMT"))
        self.assertIsNone(app.lookup_length({"chrMT": 16569}, "chrM"))
        self.assertEqual(app.lookup_length({"chrM": 16571, "chrMT": 16569}, "chrMT"), 16569)

    def test_synthetic_human_chain_patch_does_not_force_hs1_identity(self):
        self.assertNotIn("hg38ToHs1.over.chain.gz", patcher.PATCHES)
        self.assertNotIn("hs1ToHg38.over.chain.gz", patcher.PATCHES)
        self.assertNotIn("hg19ToHs1.over.chain.gz", patcher.PATCHES)
        self.assertNotIn("hs1ToHg19.over.chain.gz", patcher.PATCHES)


if __name__ == "__main__":
    unittest.main()
