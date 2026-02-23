.PHONY: all run lock clean

all: run

run:
	Rscript R/run_all.R

lock:
	Rscript -e "if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv', repos='https://cloud.r-project.org'); renv::snapshot(prompt = FALSE)"

clean:
	rm -f outputs/*.rds outputs/*.qs2 outputs/*.fst outputs/*.csv outputs/*.png outputs/*.txt outputs/*.yml
