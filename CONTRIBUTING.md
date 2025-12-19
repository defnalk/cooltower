# Contributing to cooltower

Thank you for considering a contribution! This guide covers everything you need
to submit high-quality pull requests.

---

## Table of Contents

1. [Development setup](#development-setup)
2. [Code standards](#code-standards)
3. [Testing](#testing)
4. [Pull request process](#pull-request-process)
5. [Issue guidelines](#issue-guidelines)

---

## Development setup

```bash
git clone https://github.com/defnalk/cooltower.git
cd cooltower
make install        # installs in editable mode with all dev extras
make example        # verify the example runs cleanly
```

---

## Code standards

### Type hints
Every function signature and class attribute must be fully annotated.
`mypy --strict` must pass with zero errors.

```python
# Good
def humidity_ratio(T_db: float, T_wb: float, P: float = 101_325.0) -> float: ...

# Bad — missing return type
def humidity_ratio(T_db, T_wb, P=101_325.0): ...
```

### Docstrings
Use Google style on every public function, class, and module.
Include `Args`, `Returns`, and `Raises` sections.
Add a one-line `Example` where it aids understanding.

```python
def saturation_pressure(T_db: float) -> float:
    """Compute saturation vapour pressure via the Buck (1981) formula.

    Args:
        T_db: Dry-bulb temperature [°C].

    Returns:
        Saturation vapour pressure [Pa].

    Raises:
        ValueError: If T_db is outside [-100, 200] °C.

    Example:
        >>> round(saturation_pressure(25.0), 1)
        3169.9
    """
```

### `__all__`
Every module must define `__all__` listing all public names.

### Constants
No hardcoded numeric literals inside function bodies.
Add new physical constants to `cooltower/constants.py`.

### Logging
Use `logging.getLogger(__name__)` — no `print` statements.

### Error messages
Make error messages actionable: state the bad value, the constraint, and a hint.

```python
# Good
raise ValueError(
    f"Wet-bulb temperature T_wb = {T_wb} °C cannot exceed "
    f"dry-bulb temperature T_db = {T_db} °C."
)

# Bad
raise ValueError("Invalid temperature")
```

### Formatting & linting
```bash
make format      # auto-fix with ruff
make check       # lint + format check + mypy (must all pass before PR)
```

---

## Testing

- All new functionality must have corresponding unit tests in `tests/unit/`.
- New pipeline behaviour must have integration tests in `tests/integration/`.
- Coverage must remain ≥ 90 % (`pytest --cov-fail-under=90`).
- Use `@pytest.mark.parametrize` for boundary/edge cases.
- Use fixtures in `conftest.py` for shared state; do not duplicate setup.

```bash
make test            # full suite
make test-unit       # fast feedback during development
```

---

## Pull request process

1. **Branch** — branch from `main`: `git checkout -b feat/your-feature-name`
2. **Commit messages** — use the conventional commits format:
   - `feat: add Merkel number calculation`
   - `fix: guard against zero denominator in solve_air_flow_rate`
   - `docs: add wet-bulb iteration example`
   - `test: parametrize saturation_pressure boundary cases`
3. **Run checks locally** — `make check && make test` must pass cleanly.
4. **Open a PR** against `main` with:
   - A one-paragraph description of what changed and why.
   - A reference to any related issue (`Closes #42`).
5. **CI** — all GitHub Actions checks must be green before merge.
6. **Review** — at least one approving review is required.
7. **Changelog** — add an entry under `[Unreleased]` in `CHANGELOG.md`.

---

## Issue guidelines

### Bug reports
Include:
- Python version (`python --version`)
- `cooltower` version (`pip show cooltower`)
- Minimal reproducible example
- Expected vs. actual output / traceback

### Feature requests
Include:
- The engineering motivation (equation, reference, use-case)
- Proposed function signature and return type
- Any relevant references (textbook, standard, paper)

---

## Licensing

By contributing, you agree that your contributions will be licensed under
the MIT License that covers this project.
