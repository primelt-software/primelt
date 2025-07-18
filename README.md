Certainly! Here is a `README.md` text you can copy and paste into your repository:

# PRIMELT

**PRIMELT** is a Python package and command-line tool for performing primary magma calculations with PRIMELT3-P (Herzberg et al. 2023; G-cubed) and PRIMARSMELT (Hernández-Montenegro et al., 2024; JGR: Planets) calculations on geochemical datasets of terrestrial and Martian Rocks.

## Features

- Run PRIMELT3-P or PRIMARSMELT calculations from the command line or Python.
- Input and output in CSV and Excel formats.
- Customizable MgO and FeO source values.

## Requirements

- Python 3.8 or higher
- pip (Python package installer)

The following Python packages will be installed automatically if missing:
- pandas
- matplotlib
- numpy
- openpyxl
- scipy

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/primelt-software/primelt
   cd primelt
   ```

2. **Install the package:**
   ```bash
   pip install .
   ```

## Usage

### Command Line

After installation, use the primelt command:

**PRIMELT3-P:**
```bash
primelt primelt3 path/to/input_file.csv --mgo 38.12 --feo 8.02 --output path/to/output_dir
```

**PRIMARSMELT:**
```bash
primelt primarsmelt path/to/input_file.csv --mgo 30.2 --feo 17.9 --output path/to/output_dir
```

- `primelt3` or `primarsmelt`: Calculation mode (required)
- `path/to/input_file.csv`: Path to your input CSV file (required)
- `--mgo`: MgO value for source (optional, default: 38.12 for primelt3, 30.2 for primarsmelt)
- `--feo`: FeO value for source (optional, default: 8.02 for primelt3, 17.9 for primarsmelt)
- `--output`: Output directory (optional, default: input file directory)

## Input File Format

- Input files should be in CSV format with columns matching the expected geochemical data.
- See the `examples/` folder or documentation for sample input files.

## Output

- Results are saved as CSV and Excel files in the specified output directory.
- Output files include `_afm.csv`, `_batch.csv`, and `_all.xlsx` with different result tables.

## License

MIT License

