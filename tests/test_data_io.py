"""Automated tests for molecular CSV input and output."""

from pathlib import Path

import pandas as pd
import pytest

from chemoinformatics_workbench.data_io import process_csv


def write_test_csv(path: Path) -> None:
	"""Write a temporary molecular CSV used by the tests."""

	path.write_text(
		"record_id,name,smiles,source,domain\n"
		"MOL-001,Ethanol,C(C)O,Synthetic test data,control\n"
		"MOL-002,Invalid ring,C1CC,Synthetic negative control,invalid structure\n"
		"MOL-003,Missing structure,,Synthetic negative control,missing structure\n",
		encoding="utf-8",
	)


def test_process_csv_creates_validated_output(tmp_path):
	"""Processing a CSV should create a validated file and return its data."""

	input_path = tmp_path / "raw" / "molecules.csv"
	output_path = tmp_path / "processed" / "validated_molecules.csv"

	input_path.parent.mkdir(parents=True)
	write_test_csv(input_path)

	result = process_csv(input_path, output_path)

	assert output_path.is_file()
	assert len(result) == 3
	assert result.loc[0, "canonical_smiles"] == "CCO"
	assert bool(result.loc[0, "is_valid"]) is True
	assert bool(result.loc[1, "is_valid"]) is False
	assert bool(result.loc[2, "is_valid"]) is False


def test_blank_smiles_receives_specific_error(tmp_path):
	"""An empty CSV cell should remain empty and receive a clear error."""

	input_path = tmp_path / "molecules.csv"
	output_path = tmp_path / "validated.csv"

	write_test_csv(input_path)

	result = process_csv(input_path, output_path)

	assert result.loc[2, "smiles"] == ""
	assert result.loc[2, "validation_error"] == "The SMILES string is empty."


def test_saved_csv_contains_validation_columns(tmp_path):
	"""The exported CSV should contain the calculated validation fields."""

	input_path = tmp_path / "molecules.csv"
	output_path = tmp_path / "validated.csv"

	write_test_csv(input_path)
	process_csv(input_path, output_path)

	saved_data = pd.read_csv(
		output_path,
		dtype=str,
		keep_default_na=False,
	)

	assert "is_valid" in saved_data.columns
	assert "canonical_smiles" in saved_data.columns
	assert "validation_error" in saved_data.columns


def test_input_and_output_paths_must_differ(tmp_path):
	"""The raw input file should never be overwritten."""

	input_path = tmp_path / "molecules.csv"
	write_test_csv(input_path)

	with pytest.raises(
		ValueError,
		match="Input and output CSV paths must be different.",
	):
		process_csv(input_path, input_path)


def test_missing_input_file_is_rejected(tmp_path):
	"""A missing source file should produce a clear error."""

	missing_path = tmp_path / "missing.csv"
	output_path = tmp_path / "validated.csv"

	with pytest.raises(
		FileNotFoundError,
		match="Input CSV file not found:",
	):
		process_csv(missing_path, output_path)
