"""Input and output utilities for molecular datasets."""

from pathlib import Path

import pandas as pd

from chemoinformatics_workbench.pipeline import validate_molecular_dataset


def process_csv(
	input_path: str | Path,
	output_path: str | Path,
) -> pd.DataFrame:
	"""Read, validate and export a molecular dataset without altering the raw file."""

	input_path = Path(input_path)
	output_path = Path(output_path)

	if not input_path.is_file():
		raise FileNotFoundError(f"Input CSV file not found: {input_path}")

	if input_path.resolve() == output_path.resolve():
		raise ValueError("Input and output CSV paths must be different.")

	raw_data = pd.read_csv(
		input_path,
		dtype=str,
		keep_default_na=False,
		encoding="utf-8",
	)

	validated_data = validate_molecular_dataset(raw_data)

	output_path.parent.mkdir(
		parents=True,
		exist_ok=True,
	)

	validated_data.to_csv(
		output_path,
		index=False,
		encoding="utf-8",
	)

	return validated_data
