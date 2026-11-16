"""Tests for molecular charge normalization."""

from chemoinformatics_workbench.normalization import normalize_charge


def test_neutral_molecule_remains_unchanged():
	"""A neutral molecular structure should not be altered."""

	result = normalize_charge("CCO")

	assert result.is_valid is True
	assert result.status == "unchanged"
	assert result.normalized_smiles == "CCO"
	assert result.original_formal_charge == 0
	assert result.normalized_formal_charge == 0
	assert result.was_changed is False
	assert result.error is None


def test_acetate_generates_neutral_acetic_acid_candidate():
	"""A protonatable anion should generate a neutral candidate."""

	result = normalize_charge("CC(=O)[O-]")

	assert result.is_valid is True
	assert result.status == "normalized"
	assert result.normalized_smiles == "CC(=O)O"
	assert result.original_formal_charge == -1
	assert result.normalized_formal_charge == 0
	assert result.was_changed is True
	assert result.error is None


def test_glycine_zwitterion_generates_neutral_candidate():
	"""Opposite internal charges may be normalized without changing net charge."""

	result = normalize_charge("[NH3+]CC(=O)[O-]")

	assert result.is_valid is True
	assert result.status == "normalized"
	assert result.normalized_smiles == "NCC(=O)O"
	assert result.original_formal_charge == 0
	assert result.normalized_formal_charge == 0
	assert result.was_changed is True


def test_permanent_cation_retains_residual_charge():
	"""A quaternary ammonium charge should not be artificially removed."""

	result = normalize_charge("C[N+](C)(C)C")

	assert result.is_valid is True
	assert result.status == "residual_charge"
	assert result.normalized_smiles == "C[N+](C)(C)C"
	assert result.original_formal_charge == 1
	assert result.normalized_formal_charge == 1
	assert result.was_changed is False
	assert result.error is None


def test_multifragment_structure_requires_parent_selection_first():
	"""Charge normalization should not silently remove counterions."""

	result = normalize_charge("[Na+].CC(=O)[O-]")

	assert result.is_valid is True
	assert result.status == "multifragment_rejected"
	assert result.normalized_smiles is None
	assert result.was_changed is False
	assert result.error is not None


def test_invalid_smiles_cannot_be_normalized():
	"""A malformed molecular structure should produce a traceable failure."""

	result = normalize_charge("C1CC")

	assert result.is_valid is False
	assert result.status == "invalid"
	assert result.normalized_smiles is None
	assert result.error is not None


def test_empty_smiles_cannot_be_normalized():
	"""An empty molecular representation should be rejected explicitly."""

	result = normalize_charge("   ")

	assert result.is_valid is False
	assert result.status == "invalid"
	assert result.normalized_smiles is None
	assert result.error == "The SMILES string is empty."
