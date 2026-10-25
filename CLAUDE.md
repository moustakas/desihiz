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
desi_hiz_merge --outfn $YOUR_OUTPUT_DIR/desi-protosteel.fits --img protosteel --numproc 32

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

# Run FastSpecFit on a merged catalog (reads/writes under $DESIHIZ_DIR/{img})
desi_hiz_fastspecfit --img protosteel --fastspec --mp 128
desi_hiz_fastspecfit --img protosteel --merge --coadd-type cumulative --mp 24
desi_hiz_fastspecfit --img odin --fastspec --vi-redshifts --mp 128

# Rerun Redrock with a custom template directory and gather the results
desi_hiz_redrock --img odin --run-redrock --template-dir $MY_TEMPLATE_DIR
desi_hiz_redrock --img odin --gather-redrock
```

## Code architecture

### Module layout (`py/desihiz/`)

**Core merge pipeline:**
- `hizmerge_io.py` — central I/O hub. Defines `allowed_imgs` (`odin`, `suprime`, `clauds`, `hscwide`, `ibis`, `merian`, `protosteel`) and `allowed_img_cases` per survey. Contains essentially all shared logic: path resolution (`get_img_dir`, `get_specdirs`, `get_coaddfns`), photometry table construction (`get_phot_table`), spectra reading (`get_spec_table`), VI ingestion, extinction corrections, and `merge_cases` / `build_hs` which assemble the final FITS output.
- `hizmerge_odin.py`, `hizmerge_suprime.py`, `hizmerge_clauds.py`, `hizmerge_hscwide.py`, `hizmerge_ibis.py`, `hizmerge_merian.py`, `hizmerge_protosteel.py` — survey-specific modules. Each provides `get_{survey}_{case}_infos()` functions that return target filenames, selection band metadata, and photometric catalog paths for each observation round (case).

**Extras pipeline** (post-merge value-adds):
- `extras_rr_cnn.py` — reads Redrock outputs (`read_zscan`) and CNN classifications
- `extras_phot_continuum.py` — spectrophotometric continuum fitting
- `extras_zelda.py` — Lyman-alpha line profile fitting using the `Lya_zelda` package

**FastSpecFit / Redrock post-processing** (operates on merged catalogs, img-parameterized):
- `fastspecfit_io.py` — runs FastSpecFit (fastspec/fastphot) on any img's healpix or (protosteel-only, today) cumulative tile coadds; merges per-coadd outputs; builds QA figures and an HTML QA browser. Catalogs and all outputs live under `$DESIHIZ_DIR/{img}/...` for every img (not just protosteel) — the merge catalog's *input* location and the derived-product *output* location are the same tree. The photometric extension to use (`PHOTINFO` vs `PHOTV2INFO`) is detected at runtime from the catalog's actual HDUs, not hardcoded per img. Cumulative/tile-coadd discovery loops over `hizmerge_io.get_img_cases(img)` + `get_specprod()`/`get_expids()`, so it is not protosteel-specific code — it simply finds nothing for a healpix-only specprod like `loa` (odin/suprime/clauds).
- `data/{odin,suprime,clauds,protosteel}-photinfo.yaml` — FastSpecFit photometric-parameter configs (bands, filters, flux columns) per img, checked into the repo and loaded via a package-relative path (works under `pip install -e .`).
- `redrock_io.py` — reruns Redrock on an img's healpix coadds (optionally with a custom `RR_TEMPLATE_DIR`, e.g. for testing new high-z templates) and gathers the results into a catalog row-matched to the merged catalog. Infrastructure only: building custom templates and validating against the VI subset is separate downstream work this unblocks.

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
- `$DESIHIZ_DIR` — root of the desihiz merged catalogs and FastSpecFit/Redrock derived products, for every img (defaults to `$DESI_ROOT/users/ioannis/desihiz` if unset). Note the underlying spectra FastSpecFit/Redrock actually read are NERSC-only regardless: odin/suprime/clauds coadds live under `$DESI_ROOT/users/raichoor/laelbg/loa` (~638GB), protosteel's under `$DESI_ROOT/spectro/redux/tertiary{51,52,55}` — so `desi_hiz_fastspecfit`/`desi_hiz_redrock` runs need to happen at NERSC even once `$DESIHIZ_DIR` itself is mirrored elsewhere.

### Observation rounds ("cases")

Each imaging survey has one or more cases reflecting DESI tertiary program observation rounds:
- `cosmos_yr1`, `cosmos_yr2`, `cosmos_yr3` — COSMOS field, years 1–3
- `xmmlss_yr2`, `xmmlss_yr4` — XMM-LSS field
- `ra130d5`, `ra140`, `cosmos` — protosteel (tertiary programs 0051, 0052, 0055; named by field center)

### protosteel notes

`protosteel` is the pilot for the DESI Steel weak-lensing redshift calibration sample: faint HSC-Wide galaxies selected on i-band magnitude only (22 < i_HSC < 23.5, no color cuts), `TERTIARY_TARGET = "STEEL"`.

Key differences from other surveys:
- Spectra live in the standard DESI spectro/redux tree (`$DESI_ROOT/spectro/redux/tertiary{51,52,55}/healpix/special/other/`) rather than the raichoor custom tree. Each program has its own specprod.
- Cases are named by HSC-Wide field center (`ra130d5`=RA~130.5°, `ra140`=RA~140°, `cosmos`=COSMOS field). Future fields (including one XMM-LSS and two more random HSC-Y3 pointings) will follow the same convention.
- `get_specdirs()` and `get_coaddfns()` return early for protosteel to use this non-standard path.
- `get_expids()` returns early for protosteel, reading `exposures-tertiary{NN}.fits` directly (FITS, not CSV).
- Photometry comes from per-case HSC parent catalogs (not a single all-sky file); see `_large_phot_fn` in `hizmerge_protosteel.py`. Row-selection via fitsio avoids loading the full catalog.
- COSMOS2020/CLAUDS photo-z are applied only to the `cosmos` case (the other two fields lie outside the COSMOS footprint). `get_clauds_fn()` returns `None` for non-cosmos cases.
- Fluxes are converted from nJy → nanomaggies on read and stored as `FLUX_{G,R,I,Z,Y}` / `FLUX_IVAR_{band}` and `FIBERFLUX_{band}` / `FIBERFLUX_IVAR_{band}`, matching the hscwide/suprime convention.
- No Galactic extinction columns (user computes these externally).
