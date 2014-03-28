FILES=$(shell find familyanalyzer -type 'f')
TESTFILES=$(shell find test -type 'f')

.PHONY: test

all: install permissions test

install:
	@echo "Updating familyanalyzer"
	@pip install --upgrade .

clean:
	@rm -r build/ dist/ familyanalyzer.egg-info/

uninstall:
	@echo "uninstalling familyanalyzer"
	@pip uninstall -y familyanalyzer

permissions:
	@chmod 644 $(FILES) $(TESTFILES)
	@chmod 755 bin/familyanalyzer

test:
	@python -m unittest discover test/ 'test_*.py'
