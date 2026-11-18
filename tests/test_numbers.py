"""Typing 1, 2, 3 ... after a numbered list loads that line."""
import io
import contextlib
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pymol import cmd                                          # noqa: E402
from prospector.core import cache, http, report, search        # noqa: E402
from prospector.pymol_layer import commands as base            # noqa: E402
from tests import fixture_transport                            # noqa: E402


class NumbersTest(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying
        base.register()

    def tearDown(self):
        http.transport = http._urllib_transport
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def run_line(self, line):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            cmd.do(line)
        return buffer.getvalue()

    def test_number_before_any_list_explains(self):
        del base._choices[:]
        self.assertIn("There is no list yet", self.run_line("1"))

    def test_ambiguous_list_numbers_load(self):
        result = search.find("insulin")
        if result.kind != "ambiguous":
            self.skipTest("no recording of an ambiguous search")
        picks = report.choices(result)
        self.assertTrue(picks[0].startswith("grab "))
        text = report.render(result, "grab")
        self.assertIn("just type its number", text)
        self.assertIn("->  grab ", text)
        base._set_choices(result)
        ran = []
        original = cmd.do
        try:
            cmd.do = lambda line, *a, **k: ran.append(line) if line == picks[1] \
                else original(line, *a, **k)
            self.run_line("2")
        finally:
            cmd.do = original
        self.assertEqual(ran, [picks[1]])

    def test_number_out_of_range(self):
        base._choices[:] = ["grab 1UBQ", "grab 4HHB"]
        self.assertIn("from 1 to 2", self.run_line("7"))


if __name__ == "__main__":
    unittest.main()
