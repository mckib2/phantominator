
# Remove any existing distribution archives
rm -rf dist

# Generate distribution archives
pip install --upgrade setuptools wheel
python setup.py sdist bdist_wheel

# Upload
pip install --upgrade twine
python -m twine upload dist/*
