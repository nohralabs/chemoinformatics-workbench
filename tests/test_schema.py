"""Automated tests for the molecular-data schema."""

import pandas as pd
import pytest

from chemoinformatics_workbench.schema import (
	MolecularDataSchemaError,
	validate_input_schema,
)


def make_valid_dataset() -> pd.DataFrame:
	"""Create a small valid molecular dataset for schema testing."""

	return pd.DataFrame(
		{
			"record_id": ["MOL-001", "MOL-002"],
			"name": ["Ethanol", "Propylene carbonate"],
			"smiles": ["CCO", "CC1COC(=O)O1"],
			"source": ["Synthetic test data", "Synthetic test data"],
			"domain": ["control", "cyclic carbonate"],
		}
	)


def test_valid_dataset_passes_schema_validation():
	"""A complete dataset with unique identifiers should be accepted."""

	data = make_valid_dataset()

	validate_input_schema(data)


def test_non_dataframe_is_rejected():
	"""An object that is not a pandas DataFrame should be rejected."""

	with pytest.raises(
		TypeError,
		match="The molecular dataset must be a pandas DataFrame.",
	):
		validate_input_schema(["CCO", "CC1COC(=O)O1"])


def test_missing_required_column_is_rejected():
	"""A dataset missing a required column should be rejected."""

	data = make_valid_dataset().drop(columns=["source"])

	with pytest.raises(
		MolecularDataSchemaError,
		match="Missing required columns: source",
	):
		validate_input_schema(data)


def test_empty_dataset_is_rejected():
	"""A dataset with the correct columns but no records should be rejected."""

	data = make_valid_dataset().iloc[0:0]

	with pytest.raises(
		MolecularDataSchemaError,
		match="The molecular dataset contains no records.",
	):
		validate_input_schema(data)


def test_blank_record_id_is_rejected():
	"""Every molecular record should have a non-empty identifier."""

	data = make_valid_dataset()
	data.loc[0, "record_id"] = "   "

	with pytest.raises(
		MolecularDataSchemaError,
		match="Every molecular record must have a non-empty record_id.",
	):
		validate_input_schema(data)


def test_duplicate_record_ids_are_rejected():
	"""Molecular record identifiers should be unique."""

	data = make_valid_dataset()
	data.loc[1, "record_id"] = "MOL-001"

	with pytest.raises(
		MolecularDataSchemaError,
		match="Duplicate record_id values: MOL-001",
	):
		validate_input_schema(data)


def test_missing_provenance_is_rejected():
	"""Every record should retain its source information."""

	data = make_valid_dataset()
	data.loc[0, "source"] = None

	with pytest.raises(
		MolecularDataSchemaError,
		match="Column 'source' contains missing values.",
	):
		validate_input_schema(data)
