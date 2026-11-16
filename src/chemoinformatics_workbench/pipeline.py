"""Dataset-level molecular processing pipelines."""

import pandas as pd

from chemoinformatics_workbench.schema import (
    FRAGMENT_ANALYSIS_OUTPUT_COLUMNS,
    VALIDATION_OUTPUT_COLUMNS,
    validate_input_schema,
)
from chemoinformatics_workbench.standardization import analyze_fragments
from chemoinformatics_workbench.validation import validate_smiles


def validate_molecular_dataset(data: pd.DataFrame) -> pd.DataFrame:
    """Validate and characterize every molecular record."""

    validate_input_schema(data)

    processed_data = data.copy(deep=True)

    validation_results = [
        validate_smiles(smiles)
        for smiles in processed_data["smiles"]
    ]

    validation_table = pd.DataFrame(
        {
            "is_valid": [
                result.is_valid
                for result in validation_results
            ],
            "canonical_smiles": [
                result.canonical_smiles
                for result in validation_results
            ],
            "validation_error": [
                result.error
                for result in validation_results
            ],
        },
        index=processed_data.index,
    )

    processed_data[list(VALIDATION_OUTPUT_COLUMNS)] = validation_table

    fragment_results = [
        analyze_fragments(smiles)
        for smiles in processed_data["smiles"]
    ]

    fragment_table = pd.DataFrame(
        {
            "fragment_count": [
                result.fragment_count
                for result in fragment_results
            ],
            "is_multifragment": [
                result.is_multifragment
                for result in fragment_results
            ],
            "fragment_smiles": [
                "|".join(result.fragment_smiles)
                for result in fragment_results
            ],
            "largest_fragment_smiles": [
                result.largest_fragment_smiles
                for result in fragment_results
            ],
            "total_formal_charge": pd.array(
    [
        result.total_formal_charge
        for result in fragment_results
    ],
    dtype="Int64",
),
            "fragment_analysis_error": [
                result.error
                for result in fragment_results
            ],
        },
        index=processed_data.index,
    )

    processed_data[
        list(FRAGMENT_ANALYSIS_OUTPUT_COLUMNS)
    ] = fragment_table

    return processed_data