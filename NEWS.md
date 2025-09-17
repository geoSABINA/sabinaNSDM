# sabinaNSDM 1.1.0 (2025-09-18)

## Added
- **Spatial cross-validation (optional)** via `blockCV::cv_spatial()`.
  Available in `NSDM.Global()`, `NSDM.Regional()`, and `NSDM.Covariate()`
  Controlled by the new argument `spatialCV = list(k = <int>, size = <numeric>|NULL)`.
  This option builds a CV table for BIOMOD2, enforces class balance per TEST fold, and ignores `CV.nb.rep` and `CV.perc` when enabled. 
  It allows users to evaluate models with folds that respect spatial structure, reducing over-optimistic performance estimates.
- **Single-scale modeling option**: it is now possible to obtain results for a single model by providing data through the `regional` arguments in `NSDM.InputData` (leaving the global scale as `NULL`). 
  In this case, the workflow runs with `NSDM.InputData() → NSDM.FormattingData() → NSDM.SelectCovariates() → NSDM.Regional()`, allowing faster and simpler (non-nested) modeling, although the main focus of the package remains on nested workflows.

## Changed
- **Occurrence and absences (if provided) thinning** in `NSDM.FormattingData()` has been reimplemented internally, replacing the previous dependency on **ecospat**. The new approach is faster and more reproducible than before.

## Fixed
- **Spatial CV integration**: added defensive checks to stop with clear errors if a TEST fold is empty or lacks both classes.

## Documentation
- New `@param spatialCV` describing how to activate spatial CV, required fields, units for `size`, and its precedence over `CV.nb.rep` and `CV.perc`.

---

# sabinaNSDM 1.0.1 (2024-12-16)

## Added
- **Coefficient of variation (sd/mean) of ensemble probabilities**: introduced an additional output layer (`EMcv`) , stored alongside ensemble results to quantify uncertainty.

## Changed
- Cleaner console output and revised directory structure for saved absences/background.

## Fixed
- Correct behavior when no new scenarios are provided.

---

# sabinaNSDM 1.0.0 (2024-09-12)

## Added
- First public release of **sabinaNSDM**.
- Provides the complete NSDM workflow:
  - Data input and formatting.
  - Covariate selection.
  - Global and regional model fitting.
  - Covariate and multiply strategies.

