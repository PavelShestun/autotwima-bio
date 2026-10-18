# autotwima-bio

Implementation of **Autotwima Automorphic-Twisted Magmas and $\Phi$-Groups** (verified in Lean 4) for bioinformatics coordinate transformations and DNA knot invariants.

## Requirements & Toolchain
- **uv** (Package manager & environment runner)
- **Ruff** (Linter & formatter)
- **ty** (Astral static type checker)
- **pytest** (Testing framework)

## Usage

```bash
uv sync                       # create .venv and install dependencies
uv run pytest                 # run tests
uv run ruff check .           # lint
uv run ruff format --check .  # check formatting
uv run ty check               # type check
uv run python data/download_data.py   # download phiX174 FASTA
uv run python scripts/benchmark.py    # run bioinformatics benchmark
```
