"""Automated tests for molecular SMILES validation."""

from chemoinformatics_workbench.validation import validate_smiles


def test_valid_smiles_is_canonicalized():
	"""A valid non-canonical SMILES should return its canonical form."""

	result = validate_smiles("C(C)O")

	assert result.is_valid is True
	assert result.canonical_smiles == "CCO"
	assert result.error is None


def test_chiral_smiles_preserves_stereochemistry():
	"""A defined stereocentre should remain present after canonicalization."""

	ibuprofen = "CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1"
	result = validate_smiles(ibuprofen)

	assert result.is_valid is True
	assert result.canonical_smiles is not None
	assert "@" in result.canonical_smiles
	assert result.error is None


def test_charged_smiles_is_accepted():
	"""A chemically valid zwitterionic representation should be accepted."""

	glycine_zwitterion = "[NH3+]CC(=O)[O-]"
	result = validate_smiles(glycine_zwitterion)

	assert result.is_valid is True
	assert result.canonical_smiles == glycine_zwitterion
	assert result.error is None


def test_unclosed_ring_is_rejected():
	"""A SMILES containing an unclosed ring should be rejected."""

	result = validate_smiles("C1CC")

	assert result.is_valid is False
	assert result.canonical_smiles is None
	assert result.error == "RDKit could not parse or sanitize this structure."


def test_empty_smiles_is_rejected():
	"""An empty or whitespace-only value should be rejected explicitly."""

	result = validate_smiles("   ")

	assert result.is_valid is False
	assert result.canonical_smiles is None
	assert result.error == "The SMILES string is empty."


def test_non_text_value_is_rejected():
	"""A non-string input should fail cleanly instead of crashing."""

	result = validate_smiles(123)

	assert result.is_valid is False
	assert result.canonical_smiles is None
	assert result.error == "The SMILES value must be a string."
