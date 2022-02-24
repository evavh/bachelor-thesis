run:
	python3 simulation/main.py
epaper:
	epaper -o Uni/Code/simulation simulation/*.py
	epaper -o Uni/Code/analysis analysis/*.py
test:
	testing/run_and_compare.sh
analyse:
	python3 analysis/main.py
get-data:
	cp simulation/output/* analysis/data/
