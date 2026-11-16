"""Tests for molecular-fragment analysis."""

from chemoinformatics_workbench.standardization import analyze_fragments


def test_single_molecule_has_one_fragment():
	"""A connected neutral molecule should contain one fragment."""

	result = analyze_fragments("COc1cc(C=O)ccc1O")

	assert result.is_valid is True
	assert result.fragment_count == 1
	assert result.is_multifragment is False
	assert result.largest_fragment_smiles == "COc1cc(C=O)ccc1O"
	assert result.total_formal_charge == 0
	assert result.error is None


def test_sodium_acetate_has_two_fragments():
	"""A sodium salt should retain and report both disconnected fragments."""

	result = analyze_fragments("[Na+].CC(=O)[O-]")

	assert result.is_valid is True
	assert result.fragment_count == 2
	assert result.is_multifragment is True
	assert set(result.fragment_smiles) == {"[Na+]", "CC(=O)[O-]"}
	assert result.largest_fragment_smiles == "CC(=O)[O-]"
	assert result.total_formal_charge == 0


def test_zwitterion_remains_one_fragment():
	"""Internal charges do not make a connected molecule multi-fragment."""

	result = analyze_fragments("[NH3+]CC(=O)[O-]")

	assert result.is_valid is True
	assert result.fragment_count == 1
	assert result.is_multifragment is False
	assert result.total_formal_charge == 0


def test_largest_organic_fragment_is_identified():
	"""The fragment with the most heavy atoms should be identified."""

	result = analyze_fragments("CCO.CCCCC")

	assert result.is_valid is True
	assert result.fragment_count == 2
	assert set(result.fragment_smiles) == {"CCO", "CCCCC"}
	assert result.largest_fragment_smiles == "CCCCC"


def test_invalid_smiles_returns_failed_result():
	"""A malformed structure should produce a traceable failure."""

	result = analyze_fragments("C1CC")

	assert result.is_valid is False
	assert result.fragment_count == 0
	assert result.fragment_smiles == ()
	assert result.largest_fragment_smiles is None
	assert result.error is not None


def test_empty_smiles_returns_specific_error():
	"""An empty molecular representation should be rejected explicitly."""

	result = analyze_fragments("   ")

	assert result.is_valid is False
	assert result.fragment_count == 0
	assert result.error == "The SMILES string is empty."
