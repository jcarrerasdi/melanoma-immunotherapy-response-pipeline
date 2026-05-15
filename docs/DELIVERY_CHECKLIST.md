# Delivery / GitHub checklist

## Before uploading to GitHub

- [ ] Confirm that the final thesis document has been updated with the final exact metrics.
- [ ] Confirm that `results/rf_cv_metrics_summary_final_100_no_leakage.csv` matches the internal validation values in the thesis.
- [ ] Confirm that `results/external_validation_metrics_by_dataset_final_100.csv` matches the external validation table in the thesis.
- [ ] Confirm that the gene-stability section uses the final 100-gene values: 292 unique genes, 29 genes in ≥4 folds, 13 genes in 5/5 folds.
- [ ] Remove or rename `docs/Memoria_original_uploaded.pdf` if it has not been updated.

## Recommended GitHub repository content

Keep:

- `README.md`
- `.gitignore`
- `scripts/`
- `data_raw/`
- `data_processed/metadata_*.csv`
- `results/`
- `figures/`
- `script_manifest.csv`
- `docs/` only if the final PDF is updated

Do not upload:

- `.RData`
- `.Rhistory`
- `.DS_Store`
- `__MACOSX/`
- unnecessary local cache files
- large regenerated `.rds` files unless explicitly required
