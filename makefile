.PHONY: test
all: install test

install:
	@echo "Updating familyanalyzer"
	@pip install --upgrade .

clean:
	@echo "uninstalling familyanalyzer"
	@pip uninstall -y familyanalyzer

test:
	@python -m unittest discover test/ 'test_*.py'
