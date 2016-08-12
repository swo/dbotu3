.PHONY: test
test:
	py.test

verbose:
	py.test --verbose -r w
