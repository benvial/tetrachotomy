
PROJECT_NAME := "tetrachotomy"
USER_NAME := "benvial"


BRANCH := "$(git branch --show-current)"

PROJECT_DIR := "$(realpath $PWD)"

VERSION := """$(python3 -c "import toml; print(toml.load('pyproject.toml')['project']['version'])")"""



# Echo information
info:
    @echo {{PROJECT_NAME}} version {{VERSION}}, on branch {{BRANCH}}
    @echo directory {{PROJECT_DIR}}


# List recipes
list:
    just -l

# Make conda environment
conda-env:
    mamba env create -f environment.yml

# Install the python package locally in editable mode
install:
    pip install -e .

# Lint using flake8
lint:
	flake8 --exit-zero --ignore=E501,W503 {{PROJECT_NAME}} test/*.py examples/

# Check for duplicated code
dup:
	pylint --exit-zero -f colorized --disable=all --enable=similarities {{PROJECT_NAME}}

# Reformat code
style:
	@isort .
	@black .

# Push to github
gl:
    @git add -A
    @read -p "Enter commit message: " MSG; \
    git commit -a -m "$MSG"
    @git push origin {{BRANCH}}


# Show github repository
repo:
	xdg-open https://github.com/{{USER_NAME}}/{{PROJECT_NAME}}


# Cleanup
clean:
    # cd doc && make -s clean
    rm -rf build dist

# Clean, reformat and push to github
save: clean style gl


# Check we are on the main branch
checkmain:
	@if [ "{{ BRANCH }}" != "master" ]; then exit 1; fi

# Tag and push tags
tag: style checkmain
	@echo "Version v{{VERSION}}"
	# @git add -A
	# git commit -a -m "Publish v{{VERSION}}"
	# @git push origin {{BRANCH}}
	@git tag v{{VERSION}} || echo Ignoring tag since it already exists
	@git push --tags || echo Ignoring tag since it already exists on the remote


# Create python package
package: checkmain
	@rm -f dist/*
	@python3 -m build --sdist --wheel .

# Upload to pypi
pypi: package
	@twine upload dist/*


# Publish release on pypi
publish: tag pypi

# Make logo
logo:
	MPLBACKEND=agg python dev/logo.py 1523