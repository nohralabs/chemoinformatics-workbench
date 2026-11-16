"""Automated tests for dataset-level molecular validation."""

import pandas as pd
import pytest

from chemoinformatics_workbench.pipeline import validate_molecular_dataset
from chemoinformatics_workbench.schema import MolecularDataSchemaError


def make_molecular_dataset() -> pd.DataFrame:
	"""Create a mixed dataset containing valid and invalid structures."""

	return pd.DataFrame(
		{
			"record_id": ["MOL-001", "MOL-002", "MOL-003", "MOL-004"],
			"name": [
				"Ethanol",
				"Chiral ibuprofen",
				"Glycine zwitterion",
				"Invalid ring example",
			],
			"smiles": [
				"C(C)O",
				"CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1",
				"[NH3+]CC(=O)[O-]",
				"C1CC",
			],
			"source": ["Synthetic test data"] * 4,
			"domain": [
				"control",
				"pharmaceutical",
				"charged structure",
				"negative control",
			],
		},
		index=[10, 20, 30, 40],
	)


def test_pipeline_adds_validation_columns():
	"""The pipeline should add all expected validation outputs."""

	data = make_molecular_dataset()

	validated = validate_molecular_dataset(data)

	assert "is_valid" in validated.columns
	assert "canonical_smiles" in validated.columns
	assert "validation_error" in validated.columns


def test_pipeline_returns_expected_chemical_results():
    """Valid structures should pass and the malformed ring should fail."""

    data = make_molecular_dataset()

    validated = validate_molecular_dataset(data)

    assert validated.loc[10, "canonical_smiles"] == "CCO"
    assert bool(validated.loc[20, "is_valid"]) is True
    assert "@" in validated.loc[20, "canonical_smiles"]
    assert bool(validated.loc[30, "is_valid"]) is True
    assert bool(validated.loc[40, "is_valid"]) is False
    assert pd.isna(validated.loc[40, "canonical_smiles"])
    assert validated.loc[40, "validation_error"] is not None

def test_pipeline_preserves_metadata_and_index():
	"""Record identifiers, provenance, domains and row indices should remain intact."""

	data = make_molecular_dataset()

	validated = validate_molecular_dataset(data)

	assert validated["record_id"].tolist() == data["record_id"].tolist()
	assert validated["source"].tolist() == data["source"].tolist()
	assert validated["domain"].tolist() == data["domain"].tolist()
	assert validated.index.tolist() == [10, 20, 30, 40]


def test_pipeline_does_not_modify_original_dataset():
	"""The raw input table should remain unchanged after validation."""

	data = make_molecular_dataset()
	original_copy = data.copy(deep=True)

	validate_molecular_dataset(data)

	pd.testing.assert_frame_equal(data, original_copy)


def test_pipeline_rejects_invalid_schema():
	"""The pipeline should stop before chemistry processing if the schema is invalid."""

	data = make_molecular_dataset().drop(columns=["source"])

	with pytest.raises(
		MolecularDataSchemaError,
		match="Missing required columns: source",
	):
		validate_molecular_dataset(data)
def test_pipeline_adds_fragment_analysis_columns():
    """The pipeline should add all expected fragment-analysis outputs."""

    data = make_molecular_dataset()

    processed = validate_molecular_dataset(data)

    expected_columns = {
        "fragment_count",
        "is_multifragment",
        "fragment_smiles",
        "largest_fragment_smiles",
        "total_formal_charge",
        "fragment_analysis_error",
    }

    assert expected_columns.issubset(processed.columns)


def test_pipeline_reports_multifragment_salt():
    """A salt should remain intact while its fragments are reported."""

    data = make_molecular_dataset().iloc[[0]].copy()
    data.loc[10, "name"] = "Sodium acetate"
    data.loc[10, "smiles"] = "[Na+].CC(=O)[O-]"
    data.loc[10, "domain"] = "multifragment salt"

    processed = validate_molecular_dataset(data)

    assert processed.loc[10, "fragment_count"] == 2
    assert bool(processed.loc[10, "is_multifragment"]) is True
    assert set(processed.loc[10, "fragment_smiles"].split("|")) == {
        "[Na+]",
        "CC(=O)[O-]",
    }
    assert processed.loc[10, "largest_fragment_smiles"] == "CC(=O)[O-]"
    assert processed.loc[10, "total_formal_charge"] == 0
    assert pd.isna(processed.loc[10, "fragment_analysis_error"])