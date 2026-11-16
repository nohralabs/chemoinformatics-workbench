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
@dataclass(frozen=True)
class ParentSelectionResult:
    """Structured result returned after parent-component selection."""

    original_smiles: str
    is_valid: bool
    status: str
    parent_smiles: str | None
    excluded_fragment_smiles: tuple[str, ...]
    error: str | None


def select_parent_fragment(smiles: str) -> ParentSelectionResult:
    """Select a parent component using a conservative carbon-based policy."""

    fragment_analysis = analyze_fragments(smiles)

    if not fragment_analysis.is_valid:
        return ParentSelectionResult(
            original_smiles=str(smiles),
            is_valid=False,
            status="invalid",
            parent_smiles=None,
            excluded_fragment_smiles=(),
            error=fragment_analysis.error,
        )

    molecule = Chem.MolFromSmiles(
        fragment_analysis.canonical_smiles
    )

    fragments = Chem.GetMolFrags(
        molecule,
        asMols=True,
        sanitizeFrags=True,
    )

    fragment_records = [
        {
            "molecule": fragment,
            "smiles": Chem.MolToSmiles(
                fragment,
                canonical=True,
                isomericSmiles=True,
            ),
            "heavy_atom_count": fragment.GetNumHeavyAtoms(),
            "contains_carbon": any(
                atom.GetAtomicNum() == 6
                for atom in fragment.GetAtoms()
            ),
        }
        for fragment in fragments
    ]

    if len(fragment_records) == 1:
        return ParentSelectionResult(
            original_smiles=smiles,
            is_valid=True,
            status="single_fragment",
            parent_smiles=fragment_records[0]["smiles"],
            excluded_fragment_smiles=(),
            error=None,
        )

    carbon_fragments = [
        record
        for record in fragment_records
        if record["contains_carbon"]
    ]

    if not carbon_fragments:
        return ParentSelectionResult(
            original_smiles=smiles,
            is_valid=True,
            status="no_carbon_parent",
            parent_smiles=None,
            excluded_fragment_smiles=(),
            error=(
                "No carbon-containing fragment is available "
                "for parent selection."
            ),
        )

    largest_heavy_atom_count = max(
        record["heavy_atom_count"]
        for record in carbon_fragments
    )

    largest_candidates = [
        record
        for record in carbon_fragments
        if record["heavy_atom_count"] == largest_heavy_atom_count
    ]

    if len(largest_candidates) != 1:
        return ParentSelectionResult(
            original_smiles=smiles,
            is_valid=True,
            status="ambiguous_parent",
            parent_smiles=None,
            excluded_fragment_smiles=(),
            error=(
                "Multiple carbon-containing fragments share "
                "the largest heavy-atom count."
            ),
        )

    selected_parent = largest_candidates[0]

    excluded_fragments = tuple(
        record["smiles"]
        for record in fragment_records
        if record is not selected_parent
    )

    return ParentSelectionResult(
        original_smiles=smiles,
        is_valid=True,
        status="parent_selected",
        parent_smiles=selected_parent["smiles"],
        excluded_fragment_smiles=excluded_fragments,
        error=None,
    )