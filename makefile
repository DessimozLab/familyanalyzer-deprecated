TOP := $(dir $(CURDIR)/$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))
FILES=$(shell find $(TOP) -type 'f')
TESTFILES=$(shell find test -type 'f')

.PHONY: test whereami

all: install permissions test

whereami:
	@echo $(TOP)

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
