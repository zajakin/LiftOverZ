"""Проверка допущений о реальном API pyliftover (выполняется только если библиотека установлена).

Нужна, потому что логика chrM/chrMT опирается на lo.chain_file.chains[*].source_name / source_size;
если в какой-то версии pyliftover эти поля изменятся, тест упадёт в CI, а не молча отключит защиту.
"""
from __future__ import annotations

import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

try:
    import pyliftover
    from pyliftover import LiftOver
    HAVE_REAL = getattr(pyliftover, "__file__", None) is not None
except ImportError:
    HAVE_REAL = False

if HAVE_REAL:
    spec = importlib.util.spec_from_file_location("liftover_gui_api_tested", ROOT / "liftover_gui.py")
    app = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = app
    spec.loader.exec_module(app)

CHAIN = (
    "chain 100000000 chrM 16571 + 0 16571 chrM 16569 + 0 16569 1\n3106 2 0\n13463\n\n"
    "chain 100000001 chrMT 16569 + 0 16569 chrM 16569 + 0 16569 2\n16569\n\n"
    "chain 1000 chr1 1000 + 0 1000 chr1 1002 + 0 1002 3\n500 0 2\n500\n\n"
)


@unittest.skipUnless(HAVE_REAL, "real pyliftover is not installed")
class RealPyliftoverApiTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.chain = Path(self.tmp.name) / "t.over.chain.gz"
        with gzip.open(self.chain, "wt", encoding="utf-8") as fh:
            fh.write(CHAIN)
        self.lo = LiftOver(str(self.chain))

    def tearDown(self):
        self.tmp.cleanup()

    def test_convert_coordinate_returns_list_of_tuples(self):
        res = self.lo.convert_coordinate("chr1", 600)
        self.assertEqual(res[0][0], "chr1")
        self.assertEqual(res[0][1], 602)
        self.assertEqual(res[0][2], "+")
        self.assertIsNone(self.lo.convert_coordinate("no_such_contig", 1))

    def test_chain_source_sizes_are_readable(self):
        sizes = app._chain_source_sizes(self.lo)
        self.assertIsNotNone(sizes, "lo.chain_file.chains[*].source_name/source_size are not available")
        self.assertEqual(sizes.get("chrM"), 16571)
        self.assertEqual(sizes.get("chrMT"), 16569)

    def test_mito_routing_with_real_chain(self):
        engine = app.LiftOverEngine(self.lo, None, {}, True, False, lambda m: None,
                                    lambda key, **kw: key, chain_file=str(self.chain))
        # legacy hg19 chrM: shifted by +2 after offset 3106
        self.assertEqual(engine._lift_coord("chrM", 4000)[:2], ("chrM", 4002))
        # unprefixed MT is rCRS: routed to the chrMT identity chain, no shift
        self.assertEqual(engine._lift_coord("MT", 4000)[:2], ("chrM", 4000))


if __name__ == "__main__":
    unittest.main()
