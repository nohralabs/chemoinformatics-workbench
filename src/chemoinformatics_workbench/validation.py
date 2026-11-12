"""Tools for validating and canonicalizing molecular SMILES."""

from dataclasses import dataclass

from rdkit import Chem


@dataclass(frozen=True)
class ValidationResult:
	"""Structured result returned after validating one SMILES string."""

	original_smiles: str
	is_valid: bool
	canonical_smiles: str | None
	error: str | None


def validate_smiles(smiles: str) -> ValidationResult:
	"""Validate a SMILES string and return its canonical representation."""

	if not isinstance(smiles, str):
		return ValidationResult(
			original_smiles=str(smiles),
			is_valid=False,
			canonical_smiles=None,
			error="The SMILES value must be a string.",
		)

	cleaned_smiles = smiles.strip()

	if not cleaned_smiles:
		return ValidationResult(
			original_smiles=smiles,
			is_valid=False,
			canonical_smiles=None,
			error="The SMILES string is empty.",
		)

	molecule = Chem.MolFromSmiles(cleaned_smiles)

	if molecule is None:
		return ValidationResult(
			original_smiles=smiles,
			is_valid=False,
			canonical_smiles=None,
			error="RDKit could not parse or sanitize this structure.",
		)

	return ValidationResult(
		original_smiles=smiles,
		is_valid=True,
		canonical_smiles=Chem.MolToSmiles(molecule, canonical=True),
		error=None,
	)
