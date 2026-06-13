.PHONY: install-dev install-planar install-symmetry install-docs test lint format typecheck docs package pre-commit-install check

UV ?= uv

install-dev:
	$(UV) sync --extra dev

install-planar:
	$(UV) sync --extra dev --extra planar

install-symmetry:
	$(UV) sync --extra dev --extra symmetry

install-docs:
	$(UV) sync --extra dev --extra docs

test:
	$(UV) run pytest tests

lint:
	$(UV) run ruff check

format:
	$(UV) run ruff format

typecheck:
	$(UV) run mypy

docs:
	$(UV) run sphinx -W -b html docs/source docs/_build/html

package:
	$(UV) run python -m build

pre-commit-install:
	$(UV) run pre-commit install
	$(UV) run pre-commit install --hook-type pre-push

check: lint typecheck test
