rm -rf build/html
rm -rf api/*
cd ..
python setup.py build_sphinx -w
cd docs
