"""Tools for analysing disconnected molecular fragments."""

from dataclasses import dataclass

from rdkit import Chem


@dataclass(frozen=True)
class FragmentAnalysisResult:
	"""Structured result returned after molecular-fragment analysis."""

	original_smiles: str
	is_valid: bool
	canonical_smiles: str | None
	fragment_count: int
	is_multifragment: bool
	fragment_smiles: tuple[str, ...]
	largest_fragment_smiles: str | None
	total_formal_charge: int | None
	error: str | None


def analyze_fragments(smiles: str) -> FragmentAnalysisResult:
	"""Analyse the disconnected fragments in one SMILES string."""

	if not isinstance(smiles, str):
		return FragmentAnalysisResult(
			original_smiles=str(smiles),
			is_valid=False,
			canonical_smiles=None,
			fragment_count=0,
			is_multifragment=False,
			fragment_smiles=(),
			largest_fragment_smiles=None,
			total_formal_charge=None,
			error="The SMILES value must be text.",
		)

	cleaned_smiles = smiles.strip()

	if not cleaned_smiles:
		return FragmentAnalysisResult(
			original_smiles=smiles,
			is_valid=False,
			canonical_smiles=None,
			fragment_count=0,
			is_multifragment=False,
			fragment_smiles=(),
			largest_fragment_smiles=None,
			total_formal_charge=None,
			error="The SMILES string is empty.",
		)

	molecule = Chem.MolFromSmiles(cleaned_smiles)

	if molecule is None:
		return FragmentAnalysisResult(
			original_smiles=smiles,
			is_valid=False,
			canonical_smiles=None,
			fragment_count=0,
			is_multifragment=False,
			fragment_smiles=(),
			largest_fragment_smiles=None,
			total_formal_charge=None,
			error="RDKit could not parse or sanitize this structure.",
		)

	fragments = Chem.GetMolFrags(
		molecule,
		asMols=True,
		sanitizeFrags=True,
	)

	canonical_fragments = tuple(
		Chem.MolToSmiles(
			fragment,
			canonical=True,
			isomericSmiles=True,
		)
		for fragment in fragments
	)

	largest_fragment = max(
		fragments,
		key=lambda fragment: fragment.GetNumHeavyAtoms(),
	)

	largest_fragment_smiles = Chem.MolToSmiles(
		largest_fragment,
		canonical=True,
		isomericSmiles=True,
	)

	total_formal_charge = sum(
		atom.GetFormalCharge()
		for atom in molecule.GetAtoms()
	)

	return FragmentAnalysisResult(
		original_smiles=smiles,
		is_valid=True,
		canonical_smiles=Chem.MolToSmiles(
			molecule,
			canonical=True,
			isomericSmiles=True,
		),
		fragment_count=len(fragments),
		is_multifragment=len(fragments) > 1,
		fragment_smiles=canonical_fragments,
		largest_fragment_smiles=largest_fragment_smiles,
		total_formal_charge=total_formal_charge,
		error=None,
	)
