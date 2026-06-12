# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo does

`desihiz` generates merged catalogs for DESI high-redshift (hi-z) pilot observations targeting Lyman-alpha emitters (LAEs) and related populations. It combines:
- DESI spectra from custom healpix reductions at `$DESI_ROOT/users/raichoor/laelbg`
- FIBERMAP info, tractor photometry, photometric redshifts (COSMOS2020, CLAUDS)
- Visual Inspection (VI) results and exposure metadata

Code runs at NERSC (Perlmutter). The `DESI_ROOT` environment variable points to the DESI data root at NERSC.

## Installation

```bash
pip install -e .
# or the DESI standard way:
python setup.py install
```

The package lives in `py/desihiz/` and is installed as `desihizmerge`. All `bin/` scripts are installed as executables.

## Key commands

```bash
# Generate merged catalogs (run on a NERSC login/interactive node)
desi_hiz_merge --outfn $YOUR_OUTPUT_DIR/desi-odin.fits --img odin --numproc 32
desi_hiz_merge --outfn $YOUR_OUTPUT_DIR/desi-suprime.fits --img suprime --numproc 32

# Add extras (Redrock fits, continuum params, Zelda Lya fits, CNN classifications)
desi_hiz_extras --mergefn <merge.fits> --extra redrock
desi_hiz_extras --mergefn <merge.fits> --extra zelda

# Angular clustering catalog
desi_hiz_angclust --img odin --band N419

# Predicted n(z) for LAEs
desi_hiz_nz --outfn <out.fits> --img_selection hsc-wide_v20231206

# Custom coadds from DESI specprod
desi_custom_coadds --prognum <N> --lastnight <YYYYMMDD>

# Simulate coadded spectra
desi_simcoadd

# Run Redrock on GPU nodes (see doc/redrock_cmds.ascii for exact srun invocation)
```

## Code architecture

### Module layout (`py/desihiz/`)

**Core merge pipeline:**
- `hizmerge_io.py` — central I/O hub. Defines `allowed_imgs` (`odin`, `suprime`, `clauds`, `hscwide`, `ibis`, `merian`) and `allowed_img_cases` per survey. Contains essentially all shared logic: path resolution (`get_img_dir`, `get_specdirs`, `get_coaddfns`), photometry table construction (`get_phot_table`), spectra reading (`get_spec_table`), VI ingestion, extinction corrections, and `merge_cases` / `build_hs` which assemble the final FITS output.
- `hizmerge_odin.py`, `hizmerge_suprime.py`, `hizmerge_clauds.py`, `hizmerge_hscwide.py`, `hizmerge_ibis.py`, `hizmerge_merian.py` — survey-specific modules. Each provides `get_{survey}_{case}_infos()` functions that return target filenames, selection band metadata, and photometric catalog paths for each observation round (case).

**Extras pipeline** (post-merge value-adds):
- `extras_rr_cnn.py` — reads Redrock outputs (`read_zscan`) and CNN classifications
- `extras_phot_continuum.py` — spectrophotometric continuum fitting
- `extras_zelda.py` — Lyman-alpha line profile fitting using the `Lya_zelda` package

**Angular clustering:**
- `angclust_io.py` — builds per-band photometric target catalogs for angular clustering analysis (used externally by MJW)

**Simulation:**
- `simcoadd_utils.py` — simulates coadded spectra by injecting templates into real sky backgrounds; supports multiple noise methods and continuous (z, mag) grids

**Utilities:**
- `specphot_utils.py` — spectrophotometric helpers: filter curves, Gaussian smoothing (`get_smooth`), synthetic photometry
- `laelf_utils.py` / `laelf_data.py` — Lyman-alpha luminosity function predictions and flux-limit calculations (for HEPAP preparatory study)
- `hiz_nz.py` — predicted redshift distributions n(z)
- `plot_utils.py` — shared plotting helpers
- `suprime_analysis.py`, `suprime_djs.py`, `suprime_photoff_io.py`, `suprime_photspec_io.py` — Suprime-specific photometric analysis, DJS filter selection, photometric offsets
- `hsc_griz.py` — HSC g/r/i/z broad-band selection analysis

### Data flow for `desi_hiz_merge`

1. `get_img_infos(img, case)` → calls the survey-specific `get_{img}_{case}_infos()` to get target file paths and tile/fiber info
2. `get_coaddfns(img, case)` → resolves coadd FITS paths from specprod healpix directories
3. `create_coadd_merge()` → stacks spectra per healpix pixel
4. `get_spec_table()` → extracts Redrock results, FIBERMAP columns
5. `get_phot_table()` → cross-matches to tractor photometry, adds COSMOS2020/CLAUDS photo-z
6. `merge_cases()` + `build_hs()` → combines all cases, writes final FITS with HDU per data type

### Key environment variables

- `$DESI_ROOT` — root of DESI data at NERSC (e.g., `/global/cfs/cdirs/desi`)
- `$NERSC_HOST` — set automatically at NERSC; used for config logging
- `$RR_TEMPLATE_DIR` — Redrock template directory (see `doc/redrock_cmds.ascii`)

### Observation rounds ("cases")

Each imaging survey has one or more cases reflecting DESI tertiary program observation rounds:
- `cosmos_yr1`, `cosmos_yr2`, `cosmos_yr3` — COSMOS field, years 1–3
- `xmmlss_yr2`, `xmmlss_yr4` — XMM-LSS field
