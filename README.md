# Bioinformatics Toolkit

A comprehensive web-based toolkit for DNA/RNA sequence analysis built with Streamlit.

## Features

- **GC Content Calculator** - Calculate the proportion of guanine and cytosine bases
- **Reverse / Complement** - Generate reverse, complement, and reverse complement sequences
- **DNA Translation** - Translate DNA to amino acid sequences using the standard genetic code
- **Hamming Distance** - Calculate differences between equal-length sequences
- **Edit Distance** - Compute Levenshtein distance between sequences
- **Suffix Array** - Build suffix arrays and inverse suffix arrays
- **Overlap Graph** - Construct overlap graphs for sequence assembly
- **Boyer-Moore Search** - Find pattern occurrences using efficient string matching

## Installation

```bash
pip install streamlit numpy
```

## Usage

Run the application:

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

## Input Methods

Each tool supports two input methods:
- **Type Sequence** - Enter sequences directly in the text area
- **Upload FASTA** - Upload `.fa` or `.fasta` files

## Project Structure

```
├── app.py              # Streamlit web interface
├── bio_algorithms.py   # Core bioinformatics algorithms
├── GUI.py              # Legacy Tkinter interface
└── HAPPENN_dataset.fasta
```

## Requirements

- Python 3.8+
- Streamlit
- NumPy
