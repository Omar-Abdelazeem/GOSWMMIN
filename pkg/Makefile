.PHONY: help install install-dev test test-cov lint format clean build upload-test upload docs

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

install:  ## Install package in current environment
	pip install -e .

install-dev:  ## Install package with development dependencies
	pip install -e .[dev,docs,examples]
	pre-commit install

test:  ## Run tests
	pytest tests/ -v

test-cov:  ## Run tests with coverage
	pytest tests/ -v --cov=goswmmin --cov-report=html --cov-report=term

lint:  ## Run linting checks
	flake8 goswmmin tests examples
	mypy goswmmin

format:  ## Format code with black and isort
	black goswmmin tests examples
	isort goswmmin tests examples

clean:  ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build:  ## Build package
	python -m build

upload-test:  ## Upload to test PyPI
	python -m twine upload --repository testpypi dist/*

upload:  ## Upload to PyPI
	python -m twine upload dist/*

docs:  ## Build documentation (placeholder)
	@echo "Documentation build not implemented yet"
