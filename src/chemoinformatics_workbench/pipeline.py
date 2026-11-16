"""Dataset-level molecular processing pipelines."""

import pandas as pd

import chemoinformatics_workbench.schema
from chemoinformatics_workbench.standardization import (
    analyze_fragments,
    select_parent_fragment,
)
from chemoinformatics_workbench.validation import validate_smiles


def validate_molecular_dataset(data: pd.DataFrame) -> pd.DataFrame:
    """Validate and characterize every molecular record."""

    chemoinformatics_workbench.schema.validate_input_schema(data)

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

    processed_data[
        list(
            chemoinformatics_workbench.schema
            .VALIDATION_OUTPUT_COLUMNS
        )
    ] = validation_table

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
        list(
            chemoinformatics_workbench.schema
            .FRAGMENT_ANALYSIS_OUTPUT_COLUMNS
        )
    ] = fragment_table

    parent_results = [
        select_parent_fragment(smiles)
        for smiles in processed_data["smiles"]
    ]

    parent_table = pd.DataFrame(
        {
            "parent_selection_status": [
                result.status
                for result in parent_results
            ],
            "parent_smiles": [
                result.parent_smiles
                for result in parent_results
            ],
            "excluded_fragment_smiles": [
                "|".join(result.excluded_fragment_smiles)
                for result in parent_results
            ],
            "parent_selection_error": [
                result.error
                for result in parent_results
            ],
        },
        index=processed_data.index,
    )

    processed_data[
        list(
            chemoinformatics_workbench.schema
            .PARENT_SELECTION_OUTPUT_COLUMNS
        )
    ] = parent_table

    return processed_data