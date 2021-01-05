from pathlib import Path

with open(Path(__file__).parent / "version.py") as fp: exec(fp.read())