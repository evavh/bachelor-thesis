run:
	python3 simulation/main.py
epaper:
	send_to_ereader/send.sh
test:
	testing/run_and_compare.sh
analyse:
	python3 analysis/main.py
get-data:
	cp simulation/output/* analysis/data/
