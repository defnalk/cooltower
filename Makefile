.PHONY: install test lint format typecheck clean help

PYTHON := python3
PIP    := $(PYTHON) -m pip

help:          ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

install:       ## Install package in editable mode with dev extras
	$(PIP) install --upgrade pip
	$(PIP) install -e ".[dev,examples]"

test:          ## Run full test suite with coverage (fails under 90%)
	pytest --cov=cooltower --cov-report=term-missing --cov-fail-under=90 -v

test-unit:     ## Run unit tests only
	pytest tests/unit/ -v -m unit

test-integration: ## Run integration tests only
	pytest tests/integration/ -v -m integration

lint:          ## Run ruff linter
	ruff check src/ tests/

format:        ## Auto-format with ruff
	ruff format src/ tests/

format-check:  ## Check formatting without modifying files
	ruff format --check src/ tests/

typecheck:     ## Run mypy static type checker
	mypy src/cooltower

check: lint format-check typecheck ## Run all static checks

example:       ## Run the basic analysis example
	$(PYTHON) examples/basic_analysis.py

clean:         ## Remove build artifacts and caches
	rm -rf build/ dist/ .eggs/ *.egg-info src/*.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	rm -f .coverage coverage.xml
	rm -rf .mypy_cache .ruff_cache .pytest_cache htmlcov/
