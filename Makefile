CONFIG ?= config/config.yml

.PHONY: run step1 step2 step3 step4 step5 step6 renv

run:
	Rscript R/run_all.R --config $(CONFIG)

step1:
	Rscript R/01_load_and_snapshot.R --config $(CONFIG)
step2:
	Rscript R/02_remove_contaminants_and_recluster.R --config $(CONFIG)
step3:
	Rscript R/03_build_BASE_object.R --config $(CONFIG)
step4:
	Rscript R/04_QC_and_contamination_checks.R --config $(CONFIG)
step5:
	Rscript R/05_annotation_program_scoring.R --config $(CONFIG)
step6:
	Rscript R/06_DGE_cluster_markers_and_pseudobulk.R --config $(CONFIG)

renv:
	Rscript scripts/setup_renv.R
