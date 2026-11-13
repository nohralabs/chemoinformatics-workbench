"""Dataset-level molecular validation pipelines."""

import pandas as pd

from chemoinformatics_workbench.schema import (
	VALIDATION_OUTPUT_COLUMNS,
	validate_input_schema,
)
from chemoinformatics_workbench.validation import validate_smiles


def validate_molecular_dataset(data: pd.DataFrame) -> pd.DataFrame:
	"""Validate every SMILES while preserving the original molecular records."""

	validate_input_schema(data)

	validated_data = data.copy(deep=True)

	validation_results = [
		validate_smiles(smiles) for smiles in validated_data["smiles"]
	]

	validation_table = pd.DataFrame(
		{
			"is_valid": [result.is_valid for result in validation_results],
			"canonical_smiles": [
				result.canonical_smiles for result in validation_results
			],
			"validation_error": [result.error for result in validation_results],
		},
		index=validated_data.index,
	)

	validated_data[list(VALIDATION_OUTPUT_COLUMNS)] = validation_table

	return validated_data
