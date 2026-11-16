"""Schema definitions and quality checks for molecular datasets."""

import pandas as pd


REQUIRED_INPUT_COLUMNS = (
	"record_id",
	"name",
	"smiles",
	"source",
	"domain",
)

VALIDATION_OUTPUT_COLUMNS = (
    "is_valid",
    "canonical_smiles",
    "validation_error",
)

FRAGMENT_ANALYSIS_OUTPUT_COLUMNS = (
    "fragment_count",
    "is_multifragment",
    "fragment_smiles",
    "largest_fragment_smiles",
    "total_formal_charge",
    "fragment_analysis_error",
)


class MolecularDataSchemaError(ValueError):

	"""Raised when a molecular dataset does not satisfy the required schema."""


def validate_input_schema(data: pd.DataFrame) -> None:
	"""Check that a molecular dataset has the required structure and provenance."""

	if not isinstance(data, pd.DataFrame):
		raise TypeError("The molecular dataset must be a pandas DataFrame.")

	missing_columns = [
		column for column in REQUIRED_INPUT_COLUMNS if column not in data.columns
	]

	if missing_columns:
		raise MolecularDataSchemaError(
			f"Missing required columns: {', '.join(missing_columns)}"
		)

	if data.empty:
		raise MolecularDataSchemaError("The molecular dataset contains no records.")

	blank_record_ids = (
		data["record_id"].isna()
		| data["record_id"].astype(str).str.strip().eq("")
	)

	if blank_record_ids.any():
		raise MolecularDataSchemaError(
			"Every molecular record must have a non-empty record_id."
		)

	duplicated_record_ids = data["record_id"].duplicated(keep=False)

	if duplicated_record_ids.any():
		duplicate_values = sorted(
			data.loc[duplicated_record_ids, "record_id"].astype(str).unique()
		)
		raise MolecularDataSchemaError(
			f"Duplicate record_id values: {', '.join(duplicate_values)}"
		)

	for provenance_column in ("source", "domain"):
		blank_values = (
			data[provenance_column].isna()
			| data[provenance_column].astype(str).str.strip().eq("")
		)

		if blank_values.any():
			raise MolecularDataSchemaError(
				f"Column '{provenance_column}' contains missing values."
			)
