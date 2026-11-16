"""Dataset-level molecular processing pipelines."""

import pandas as pd

import chemoinformatics_workbench.schema
from chemoinformatics_workbench.standardization import (
    analyze_fragments,
    select_parent_fragment,
)
from chemoinformatics_workbench.validation import validate_smiles
from chemoinformatics_workbench.normalization import normalize_charge

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
    normalization_results = [
        (
            normalize_charge(result.parent_smiles)
            if result.parent_smiles is not None
            else None
        )
        for result in parent_results
    ]

    normalization_table = pd.DataFrame(
        {
            "charge_normalization_status": [
                (
                    result.status
                    if result is not None
                    else "not_attempted"
                )
                for result in normalization_results
            ],
            "charge_normalized_smiles": [
                (
                    result.normalized_smiles
                    if result is not None
                    else None
                )
                for result in normalization_results
            ],
            "original_parent_formal_charge": pd.array(
                [
                    (
                        result.original_formal_charge
                        if result is not None
                        else None
                    )
                    for result in normalization_results
                ],
                dtype="Int64",
            ),
            "normalized_formal_charge": pd.array(
                [
                    (
                        result.normalized_formal_charge
                        if result is not None
                        else None
                    )
                    for result in normalization_results
                ],
                dtype="Int64",
            ),
            "charge_was_changed": [
                (
                    result.was_changed
                    if result is not None
                    else False
                )
                for result in normalization_results
            ],
            "charge_normalization_error": [
                (
                    result.error
                    if result is not None
                    else (
                        "Charge normalization was not attempted "
                        "because parent selection did not produce "
                        "a molecular structure."
                    )
                )
                for result in normalization_results
            ],
        },
        index=processed_data.index,
    )

    processed_data[
        list(
            chemoinformatics_workbench.schema
            .CHARGE_NORMALIZATION_OUTPUT_COLUMNS
        )
    ] = normalization_table
    
    return processed_data 