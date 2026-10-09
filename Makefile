PYTHON ?= python
RUN = cd pipeline && $(PYTHON) -W ignore

.PHONY: all analysis site figures report check clean

all: analysis site figures report

check:
	cd data && md5sum -c raw.md5

analysis:
	$(RUN) 01_inputs_and_gene_list.py
	$(RUN) 02_dependency_classes.py
	$(RUN) 03_single_gene.py
	$(RUN) 04_codependency.py
	$(RUN) 05_lineage_restricted.py
	$(RUN) 06_modules.py
	$(RUN) 07_literature.py
	$(RUN) 08_known_partners.py

site:
	$(RUN) 09_explorer_data.py
	$(RUN) 10_build_site.py

figures:
	$(RUN) 11_figures.py

report:
	$(RUN) 12_report.py

clean:
	rm -rf results/intermediate results/tables results/figures results/*.xlsx
