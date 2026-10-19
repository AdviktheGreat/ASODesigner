"""Tests for ViennaRNA wrappers (RNAfold, RNAplfold, RNAcofold, RNAduplex).

All tests mock subprocess calls so they run without ViennaRNA installed.
"""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from aso_designer.utils.exceptions import ViennaRNAError
from aso_designer.viennarna import (
    ViennaRNAVersion,
    _parse_rnacofold_output,
    _parse_rnaduplex_output,
    _parse_rnafold_output,
    _parse_rnaplfold_lunp,
    _parse_rnaup_stdout,
    _parse_version,
    run_rnacofold,
    run_rnaduplex,
    run_rnafold,
    run_rnaplfold,
    run_rnaup,
    verify_viennarna,
)

# ---------------------------------------------------------------------------
# Version parsing
# ---------------------------------------------------------------------------


class TestParseVersion:
    def test_standard_format(self):
        v = _parse_version("RNAfold 2.6.4")
        assert v is not None
        assert v.major == 2
        assert v.minor == 6
        assert v.patch == 4

    def test_viennarna_prefix(self):
        v = _parse_version("ViennaRNA 2.5.1")
        assert v is not None
        assert v.major == 2
        assert v.minor == 5
        assert v.patch == 1

    def test_bare_version(self):
        v = _parse_version("2.4.6")
        assert v is not None
        assert v.major == 2
        assert v.minor == 4
        assert v.patch == 6

    def test_no_version(self):
        assert _parse_version("no version here") is None

    def test_version_comparison(self):
        v = ViennaRNAVersion(2, 6, 4)
        assert v >= (2, 4, 0)
        assert v >= (2, 6, 4)
        assert not v >= (3, 0, 0)

    def test_version_string(self):
        v = ViennaRNAVersion(2, 6, 4)
        assert v.version_string == "2.6.4"


# ---------------------------------------------------------------------------
# verify_viennarna
# ---------------------------------------------------------------------------


class TestVerifyViennarna:
    @patch("aso_designer.viennarna.shutil.which")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_success(self, mock_run, mock_which):
        mock_which.return_value = "/usr/local/bin/RNAfold"
        mock_run.return_value = MagicMock(stdout="RNAfold 2.6.4\n", stderr="")
        v = verify_viennarna()
        assert v.major == 2
        assert v.minor == 6
        assert v.patch == 4

    @patch("aso_designer.viennarna.shutil.which")
    def test_not_found(self, mock_which):
        mock_which.return_value = None
        with pytest.raises(ViennaRNAError, match="RNAfold not found"):
            verify_viennarna()

    @patch("aso_designer.viennarna.shutil.which")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_old_version(self, mock_run, mock_which):
        mock_which.return_value = "/usr/local/bin/RNAfold"
        mock_run.return_value = MagicMock(stdout="RNAfold 2.3.0\n", stderr="")
        with pytest.raises(ViennaRNAError, match="required"):
            verify_viennarna()

    @patch("aso_designer.viennarna.shutil.which")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_unparseable_version(self, mock_run, mock_which):
        mock_which.return_value = "/usr/local/bin/RNAfold"
        mock_run.return_value = MagicMock(stdout="garbage\n", stderr="")
        with pytest.raises(ViennaRNAError, match="Could not parse"):
            verify_viennarna()

    @patch("aso_designer.viennarna.shutil.which")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_subprocess_failure(self, mock_run, mock_which):
        mock_which.return_value = "/usr/local/bin/RNAfold"
        mock_run.side_effect = subprocess.SubprocessError("boom")
        with pytest.raises(ViennaRNAError, match="Failed to execute"):
            verify_viennarna()


# ---------------------------------------------------------------------------
# RNAfold output parsing
# ---------------------------------------------------------------------------


class TestParseRnafoldOutput:
    def test_standard_output(self):
        mock = ">test\nAUGCAUGC\n((((...)))) (-2.30)\n"
        r = _parse_rnafold_output("AUGCAUGC", mock)
        assert r.sequence == "AUGCAUGC"
        assert r.structure == "((((...))))"
        assert r.mfe == -2.30

    def test_dna_input_converted(self):
        # run_rnafold converts DNA to RNA, so the sequence in result is RNA
        mock = ">test\nAUGCAUGC\n((((...)))) (-1.50)\n"
        r = _parse_rnafold_output("AUGCAUGC", mock)
        assert r.sequence == "AUGCAUGC"
        assert r.mfe == -1.50

    def test_long_sequence(self):
        seq = "AUGCAUGCAUGCAUGCAUGC"
        struct = "((((....))))....((..))"
        mock = f">test\n{seq}\n{struct} (-5.40)\n"
        r = _parse_rnafold_output(seq, mock)
        assert r.structure == struct
        assert r.mfe == -5.40

    def test_positive_mfe(self):
        mock = ">test\nAUG\n... (0.00)\n"
        r = _parse_rnafold_output("AUG", mock)
        assert r.mfe == 0.0

    def test_regex_fallback(self):
        # Malformed but MFE still extractable
        mock = "some junk\nenergy = (-3.20)\n"
        r = _parse_rnafold_output("AUGCAUGC", mock)
        assert r.mfe == -3.20

    def test_parse_failure(self):
        with pytest.raises(ViennaRNAError, match="Could not parse"):
            _parse_rnafold_output("AUGCAUGC", "no useful data here")


# ---------------------------------------------------------------------------
# RNAfold wrapper (mocked subprocess)
# ---------------------------------------------------------------------------


class TestRunRnafold:
    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_success(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAfold"
        mock_run.return_value = MagicMock(
            stdout=">test\nAUGCAUGC\n((((...)))) (-2.30)\n", stderr=""
        )
        r = run_rnafold("ATGCATGC")
        assert r.sequence == "AUGCAUGC"
        assert r.structure == "((((...))))"
        assert r.mfe == -2.30

    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_failure(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAfold"
        mock_run.side_effect = subprocess.SubprocessError("boom")
        with pytest.raises(ViennaRNAError, match="execution failed"):
            run_rnafold("ATGCATGC")


# ---------------------------------------------------------------------------
# RNAplfold output parsing
# ---------------------------------------------------------------------------


class TestParseRnaplfoldLunp:
    def test_standard_parse(self, tmp_path):
        lunp = tmp_path / "test_lunp"
        lunp.write_text(
            "#unpaired probabilities\n"
            "#i   p_u   u\n"
            "1    0.95  1\n"
            "2    0.87  1\n"
            "3    0.12  1\n"
            "4    0.45  1\n"
        )
        r = _parse_rnaplfold_lunp(lunp)
        assert r.positions == [1, 2, 3, 4]
        assert r.unpaired_probs == [0.95, 0.87, 0.12, 0.45]

    def test_empty_file(self, tmp_path):
        lunp = tmp_path / "test_lunp"
        lunp.write_text("#header only\n")
        with pytest.raises(ViennaRNAError, match="No unpaired probabilities"):
            _parse_rnaplfold_lunp(lunp)

    def test_skip_comments(self, tmp_path):
        lunp = tmp_path / "test_lunp"
        lunp.write_text("#comment1\n#comment2\n1    0.50  1\n#mid comment\n2    0.30  1\n")
        r = _parse_rnaplfold_lunp(lunp)
        assert r.positions == [1, 2]
        assert r.unpaired_probs == [0.50, 0.30]


# ---------------------------------------------------------------------------
# RNAplfold wrapper (mocked subprocess + temp dir)
# ---------------------------------------------------------------------------


class TestRunRnaplfold:
    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_success(self, mock_run, mock_path, tmp_path):
        mock_path.return_value = "/usr/local/bin/RNAplfold"

        # Create the expected lunp file in the temp dir
        def side_effect(*args, **kwargs):
            cwd = kwargs.get("cwd", tmp_path)
            lunp = Path(cwd) / "asodesigner_lunp"
            lunp.write_text("#unpaired probabilities\n1    0.90  1\n2    0.10  1\n")

        mock_run.side_effect = side_effect
        r = run_rnaplfold("ATGCATGC", window_size=80, max_span=50)
        assert r.positions == [1, 2]
        assert r.unpaired_probs == [0.90, 0.10]

    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_no_lunp_file(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAplfold"
        mock_run.return_value = MagicMock(stdout="", stderr="")
        with pytest.raises(ViennaRNAError, match="did not produce"):
            run_rnaplfold("ATGCATGC")

    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_failure(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAplfold"
        mock_run.side_effect = subprocess.SubprocessError("boom")
        with pytest.raises(ViennaRNAError, match="execution failed"):
            run_rnaplfold("ATGCATGC")


# ---------------------------------------------------------------------------
# RNAcofold output parsing
# ---------------------------------------------------------------------------


class TestParseRnacofoldOutput:
    def test_standard_output(self):
        mock = ">seq1\nAUGCAUGC&AUGCAUGC\n((((...))))&((((...)))) (-4.60)\n"
        r = _parse_rnacofold_output("AUGCAUGC", "AUGCAUGC", mock)
        assert r.sequence1 == "AUGCAUGC"
        assert r.sequence2 == "AUGCAUGC"
        assert r.structure == "((((...))))&((((...))))"
        assert r.mfe == -4.60

    def test_regex_fallback(self):
        mock = "junk\nenergy = (-7.20)\n"
        r = _parse_rnacofold_output("AUG", "GCA", mock)
        assert r.mfe == -7.20

    def test_parse_failure(self):
        with pytest.raises(ViennaRNAError, match="Could not parse"):
            _parse_rnacofold_output("AUG", "GCA", "no useful data")


class TestRunRnacofold:
    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_success(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAcofold"
        mock_run.return_value = MagicMock(
            stdout=">seq1\nAUGCAUGC&AUGCAUGC\n((((...))))&((((...)))) (-4.60)\n",
            stderr="",
        )
        r = run_rnacofold("ATGCATGC", "ATGCATGC")
        assert r.mfe == -4.60

    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_failure(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAcofold"
        mock_run.side_effect = subprocess.SubprocessError("boom")
        with pytest.raises(ViennaRNAError, match="execution failed"):
            run_rnacofold("ATGC", "ATGC")


# ---------------------------------------------------------------------------
# RNAduplex output parsing
# ---------------------------------------------------------------------------


class TestParseRnaduplexOutput:
    def test_standard_output(self):
        mock = ".(((((((((...)))))))))...  &...((((((....))))))...  (-12.40)\n1,15  : 30,16\n"
        r = _parse_rnaduplex_output("AUGCAUGC", "AUGCAUGC", mock)
        assert r.mfe == -12.40
        assert r.start1 == 1
        assert r.end1 == 15
        assert r.start2 == 30
        assert r.end2 == 16

    def test_parse_failure(self):
        with pytest.raises(ViennaRNAError, match="Could not parse"):
            _parse_rnaduplex_output("AUG", "GCA", "no useful data")


class TestRunRnaduplex:
    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_success(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAduplex"
        mock_run.return_value = MagicMock(
            stdout=".(((((((((...)))))))))...  &...((((((....))))))...  (-12.40)\n1,15  : 30,16\n",
            stderr="",
        )
        r = run_rnaduplex("ATGCATGC", "ATGCATGC")
        assert r.mfe == -12.40
        assert r.start1 == 1
        assert r.end1 == 15

    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_failure(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAduplex"
        mock_run.side_effect = subprocess.SubprocessError("boom")
        with pytest.raises(ViennaRNAError, match="execution failed"):
            run_rnaduplex("ATGC", "ATGC")


# ---------------------------------------------------------------------------
# RNAup parser and wrapper (mocked subprocess)
# ---------------------------------------------------------------------------

RNAUP_STDOUT = """RNAup output in file: asodesigner_unup.out
>t_1
AUGCUAGC
[region 4]
   1    0.98    0.90    -0.1053
   2    0.96    0.85    0.1625
   3    0.40    0.10    2.3026
"""


class TestParseRnaupStdout:
    def test_standard_layout(self):
        r = _parse_rnaup_stdout(RNAUP_STDOUT, "AUGCUAGC")
        assert r.sequence == "AUGCUAGC"
        assert r.positions == [1, 2, 3]
        assert r.unpaired_probs == [0.98, 0.96, 0.40]
        assert r.region_probs == [0.90, 0.85, 0.10]
        assert r.opening_dg == [-0.1053, 0.1625, 2.3026]
        assert r.region_length == 4

    def test_missing_region_block(self):
        with pytest.raises(ViennaRNAError, match="Could not parse RNAup"):
            _parse_rnaup_stdout(">t_1\nAUGCUAGC\n", "AUGCUAGC")

    def test_bad_row(self):
        bad = RNAUP_STDOUT.replace("2    0.96    0.85    0.1625", "2    0.96")
        with pytest.raises(ViennaRNAError, match="output row"):
            _parse_rnaup_stdout(bad, "AUGCUAGC")


class TestRunRnaup:
    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_success_parses_stdout(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAup"
        mock_run.return_value = MagicMock(stdout=RNAUP_STDOUT, stderr="")
        r = run_rnaup("ATGCTAGC", unpaired_region=4)
        assert r.region_length == 4
        assert r.positions == [1, 2, 3]
        assert r.opening_dg == [-0.1053, 0.1625, 2.3026]
        # T converted to U in the sequence handed to RNAup
        stdin = mock_run.call_args.kwargs["input"]
        assert "AUGCUAGC" in stdin

    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_passes_unpaired_region(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAup"
        mock_run.return_value = MagicMock(stdout=RNAUP_STDOUT, stderr="")
        run_rnaup("ATGC", unpaired_region=6)
        args = mock_run.call_args.args[0]
        assert "-u" in args
        assert args[args.index("-u") + 1] == "6"
        assert "--id" in args

    @patch("aso_designer.viennarna._get_tool_path")
    @patch("aso_designer.viennarna.subprocess.run")
    def test_run_failure(self, mock_run, mock_path):
        mock_path.return_value = "/usr/local/bin/RNAup"
        mock_run.side_effect = subprocess.SubprocessError("boom")
        with pytest.raises(ViennaRNAError, match="execution failed"):
            run_rnaup("ATGC")
