.PHONY: test build upload
test:
	py.test

verbose:
	py.test --verbose -r w

build:
	python3 setup.py bdist_wheel sdist

upload:
	twine upload dist/*

update:
	python3 setup.py sdist upload -r pypi
