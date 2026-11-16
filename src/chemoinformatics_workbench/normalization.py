"""Tools for generating charge-normalized molecular candidates."""

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


@dataclass(frozen=True)
class ChargeNormalizationResult:
	"""Structured result returned after charge normalization."""

	original_smiles: str
	is_valid: bool
	status: str
	normalized_smiles: str | None
	original_formal_charge: int | None
	normalized_formal_charge: int | None
	was_changed: bool
	error: str | None


def _calculate_formal_charge(molecule: Chem.Mol) -> int:
	"""Calculate the total formal charge of an RDKit molecule."""

	return sum(
		atom.GetFormalCharge()
		for atom in molecule.GetAtoms()
	)


def normalize_charge(smiles: str) -> ChargeNormalizationResult:
	"""Generate a charge-normalized candidate for one connected structure."""

	if not isinstance(smiles, str):
		return ChargeNormalizationResult(
			original_smiles=str(smiles),
			is_valid=False,
			status="invalid",
			normalized_smiles=None,
			original_formal_charge=None,
			normalized_formal_charge=None,
			was_changed=False,
			error="The SMILES value must be text.",
		)

	cleaned_smiles = smiles.strip()

	if not cleaned_smiles:
		return ChargeNormalizationResult(
			original_smiles=smiles,
			is_valid=False,
			status="invalid",
			normalized_smiles=None,
			original_formal_charge=None,
			normalized_formal_charge=None,
			was_changed=False,
			error="The SMILES string is empty.",
		)

	molecule = Chem.MolFromSmiles(cleaned_smiles)

	if molecule is None:
		return ChargeNormalizationResult(
			original_smiles=smiles,
			is_valid=False,
			status="invalid",
			normalized_smiles=None,
			original_formal_charge=None,
			normalized_formal_charge=None,
			was_changed=False,
			error="RDKit could not parse or sanitize this structure.",
		)

	if len(Chem.GetMolFrags(molecule)) > 1:
		return ChargeNormalizationResult(
			original_smiles=smiles,
			is_valid=True,
			status="multifragment_rejected",
			normalized_smiles=None,
			original_formal_charge=_calculate_formal_charge(molecule),
			normalized_formal_charge=None,
			was_changed=False,
			error=(
				"Parent selection is required before charge "
				"normalization of a multi-fragment structure."
			),
		)

	canonical_original = Chem.MolToSmiles(
		molecule,
		canonical=True,
		isomericSmiles=True,
	)

	original_formal_charge = _calculate_formal_charge(molecule)

	uncharger = rdMolStandardize.Uncharger()
	normalized_molecule = uncharger.uncharge(molecule)

	normalized_smiles = Chem.MolToSmiles(
		normalized_molecule,
		canonical=True,
		isomericSmiles=True,
	)

	normalized_formal_charge = _calculate_formal_charge(
		normalized_molecule
	)

	was_changed = normalized_smiles != canonical_original

	if normalized_formal_charge != 0:
		status = "residual_charge"
	elif was_changed:
		status = "normalized"
	else:
		status = "unchanged"

	return ChargeNormalizationResult(
		original_smiles=smiles,
		is_valid=True,
		status=status,
		normalized_smiles=normalized_smiles,
		original_formal_charge=original_formal_charge,
		normalized_formal_charge=normalized_formal_charge,
		was_changed=was_changed,
		error=None,
	)
