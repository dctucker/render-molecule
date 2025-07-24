# render-molecule

A simple Python API that renders molecules.

Currently supported formats:
* Input
  * SMILES
* Output
  * SVG
  * PNG

Compatible with MediaWiki external images by matching final part of the query string with a supported format (e.g. `&format=.png`).

## Installation
```
python3 -m venv .
bin/pip3 install -r requirements.txt
```

## Running
```
bin/python3 src/rest.py
```

## Usage
```
http://localhost:8555/render/?name=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&format=.png
```

Yields:

<img width="140" height="136" alt="image" src="https://github.com/user-attachments/assets/96f041fc-1aa2-4892-8782-3dbb6621d89c" />

URL encoding is recommended but not strictly required in all cases.
