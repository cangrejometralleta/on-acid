#!/bin/sh
# Builds the Package. Dependencies, lint, types, tests, then the artefact.
set -eu

cd -- "$(dirname -- "$0")"

PYTHON="${PYTHON:-python3}"
VENV="${VENV:-.venv}"

refuse() {
	printf '\342\235\214 %s\n' "$1" >&2
	printf '   %s\n' "$2" >&2
	exit 1
}

# Says where a venv keeps its executables; Windows Spells it Scripts.
venv_bin() {
	[ -d "$VENV/Scripts" ] && printf '%s/Scripts' "$VENV" || printf '%s/bin' "$VENV"
}

tool() {
	printf '%s/%s' "$(venv_bin)" "$1"
}

verify_interpreter() {
	command -v "$PYTHON" >/dev/null 2>&1 ||
		refuse "wanted a Python interpreter named $PYTHON." \
			"PYTHON=/path/to/python ./build.sh"

	"$PYTHON" -c 'import venv' >/dev/null 2>&1 ||
		refuse "wanted the venv module, which $PYTHON does not carry." \
			"install the python3-venv package for $PYTHON"
}

install_packages() {
	[ -x "$(tool python)" ] || "$PYTHON" -m venv "$VENV" ||
		refuse "wanted a virtual environment at $VENV." "$PYTHON -m venv $VENV"

	"$(tool python)" -m pip install --quiet --upgrade pip ruff mypy build ||
		refuse "wanted ruff, mypy and build in $VENV." \
			"$(tool python) -m pip install ruff mypy build"
}

# Lints rather than formats: the aligned tables are Deliberate,
# and a formatter that Reflows them Costs more than it Catches.
verify_lint() {
	"$(tool ruff)" check . ||
		refuse "wanted the lint to be clean." "$(tool ruff) check --fix ."
}

verify_types() {
	"$(tool mypy)" dissociation.py ||
		refuse "wanted the annotations to hold." "$(tool mypy) dissociation.py"
}

verify_tests() {
	"$PYTHON" -m unittest discover -p 'test_*.py' ||
		refuse "wanted every test to pass." "$PYTHON -m unittest discover -p 'test_*.py'"
}

verify_artefact() {
	"$(tool python)" -m build ||
		refuse "wanted a wheel and an sdist." "$(tool python) -m build"
}

verify_interpreter
install_packages
verify_lint
verify_types
verify_tests
verify_artefact
