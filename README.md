# gpc
Genotype Phenotype Correlation 


# Set up

Create a virtual environment. This is an optional step, otherwise the package will be installed into the active 
Python environment.

```shell
python3 -m venv venv
source venv/bin/activate
```

Install the package into the environment:

```shell
python3 -m pip install genophenocorr
```

# Run tests

Tests can be run after checking the source code from GitHub and installing the test dependencies:

```shell
cd genophenocorr
python3 -m pip install .[test]

pytest
```

# Build documentation

First, make sure you have the necessary dependencies:

```shell
cd genophenocorr
python3 -m pip install .[docs]
```
Setting up dependencies is a one-time action for a given Python environment.

Next, we can run the doc tests and build the HTML documentation:

```shell
cd docs

# Generate the API reference from Python docstrings
sphinx-apidoc --separate --module-first -d 2 -H "API reference" --follow-links -o apidocs ../src/genophenocorr

# Run the doc tests and build the documentation
make doctest html
```

The code above will run the documentation and fail if the docs are out of sync with the code.
Then, the docs will be built at `docs/build/html`
