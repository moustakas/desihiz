#!/usr/bin/env python


import os
import fitsio
import numpy as np
from astropy.table import Table
from desiutil.log import get_logger
from desihiz.hizmerge_io import (
    match_coord,
    get_init_infos,
    get_phot_fns,
)

log = get_logger()

# columns to read from the large HSC photometric catalog
_large_phot_cols = [
    "OBJECT_ID", "RA", "DEC", "I_CMODEL_MAG_CORR",
    "G_CMODEL_FLUX", "G_CMODEL_FLUXERR",
    "R_CMODEL_FLUX", "R_CMODEL_FLUXERR",
    "I_CMODEL_FLUX", "I_CMODEL_FLUXERR",
    "Z_CMODEL_FLUX", "Z_CMODEL_FLUXERR",
    "Y_CMODEL_FLUX", "Y_CMODEL_FLUXERR",
    "G_FIBER_FLUX", "G_FIBER_FLUXERR",
    "R_FIBER_FLUX", "R_FIBER_FLUXERR",
    "I_FIBER_FLUX", "I_FIBER_FLUXERR",
    "Z_FIBER_FLUX", "Z_FIBER_FLUXERR",
    "Y_FIBER_FLUX", "Y_FIBER_FLUXERR",
    "I_EXTENDEDNESS_VALUE",
    "I_CMODEL_FLAG",
    "I_MASK_BRIGHTSTAR_ANY",
    "M_I_BLENDEDNESS_ABS",
    "DNNZ_PHOTOZ_BEST", "DNNZ_PHOTOZ_RISK_BEST", "DNNZ_PHOTOZ_STD_BEST",
]

_NJY_PER_NANOMAGGY = 3631.0


_large_phot_basedir = os.path.join(
    os.getenv("DESI_ROOT", ""),
    "users", "ioannis", "desihiz", "protosteel", "phot",
)
#_large_phot_basedir = os.path.join(
#    os.getenv("DESI_ROOT", ""),
#    "users", "nweaverd", "nweaverd_desi", "Steel", "finalized_target_catalogs",
#)

_large_phot_fn = {
    "ra130d5": "hsc_icmodelmag_22-24_RA130d5_DEC000.fits.gz",
    "ra140":   "hsc_icmodelmag_22-24_RA140_DEC003.fits.gz",
    "cosmos":  "hsc_icmodelmag_22-24_COSMOS.fits.gz",
}


def get_protosteel_large_fn(case):
    """Full path to the per-case HSC photometric catalog for protosteel."""
    return os.path.join(_large_phot_basedir, _large_phot_fn[case])


def read_protosteel_large_phot(fn, objids):
    """
    Read photometric columns for specific OBJECT_IDs from the large HSC catalog.
    Uses fitsio row-selection to avoid loading the full catalog into memory.

    Args:
        fn: full path to the large HSC catalog (str)
        objids: OBJECT_ID values to retrieve (array of int)

    Returns:
        p: table with _large_phot_cols columns for the requested objects
    """
    all_objids = fitsio.read(fn, columns=["OBJECT_ID"])["OBJECT_ID"].ravel()
    rows = np.where(np.in1d(all_objids, objids))[0]
    p = Table(fitsio.read(fn, rows=rows, columns=_large_phot_cols))
    for key in p.colnames:
        p[key].name = p[key].name.upper()

    # convert nJy → nanomaggies, rename to standard FLUX/FIBERFLUX columns
    basename = os.path.basename(fn)
    for band in ["G", "R", "I", "Z", "Y"]:
        for prefix, raw_prefix in [("FLUX", "CMODEL"), ("FIBERFLUX", "FIBER")]:
            flux_raw = "{}_{}_{}" .format(band, raw_prefix, "FLUX")
            err_raw  = "{}_{}_{}".format(band, raw_prefix, "FLUXERR")
            flux_out = "{}_{}".format(prefix, band)
            ivar_out = "{}_IVAR_{}".format(prefix, band)

            p[flux_out] = p[flux_raw] / _NJY_PER_NANOMAGGY
            ivar = np.zeros(len(p), dtype=">f4")
            good = p[err_raw] > 0
            ivar[good] = (_NJY_PER_NANOMAGGY / p[err_raw][good]) ** 2
            p[ivar_out] = ivar

            del p[flux_raw], p[err_raw]

            log.info(
                "{}: convert {} and {} to {} (nanomaggies) and {}".format(
                    basename, flux_raw, err_raw, flux_out, ivar_out
                )
            )

    return p


def _get_protosteel_infos(prognum, fewcols_fn):
    """
    Shared logic for the per-case get_protosteel_cosmos_prNN_infos() functions.

    Args:
        prognum: zero-padded tertiary program number (e.g. "0051") (str)
        fewcols_fn: basename of the fewcols input catalog (str)

    Returns:
        mydict: same structure as get_init_infos(), keyed by "GRIZ"
    """
    fadir = os.path.join(
        os.getenv("DESI_ROOT"), "survey", "fiberassign", "special", "tertiary", prognum
    )

    # read fiberassign assign file
    fn = os.path.join(fadir, "tertiary-targets-{}-assign.fits".format(prognum))
    d = Table.read(fn)

    # select science targets
    sel = d["TERTIARY_TARGET"].astype(str) == "STEEL"
    d = d[sel]

    # read fewcols input catalog
    t = Table.read(os.path.join(fadir, "inputcats", fewcols_fn))
    for key in t.colnames:
        t[key].name = key.upper()

    # match fiberassign targets to fewcols by RA/Dec
    # coordinates should be identical (targets were selected from the fewcols catalog)
    iid, iit, _, _, _ = match_coord(
        d["RA"], d["DEC"],
        t["RA"], t["DEC"],
        search_radius=0.1,
        verbose=True,
    )
    if iid.size != len(d):
        log.warning(
            "tertiary{}: only {}/{} STEEL targets matched in fewcols catalog".format(
                prognum, iid.size, len(d)
            )
        )
    d, t = d[iid], t[iit]

    nrows = [len(d)]
    mydict = get_init_infos("protosteel", nrows)

    mydict["GRIZ"]["TARGETID"] = d["TARGETID"]
    mydict["GRIZ"]["TERTIARY_TARGET"] = d["TERTIARY_TARGET"]
    mydict["GRIZ"]["PHOT_RA"] = t["RA"]
    mydict["GRIZ"]["PHOT_DEC"] = t["DEC"]
    mydict["GRIZ"]["PHOT_SELECTION"] = np.full(
        len(d), "HSC_I", dtype=mydict["GRIZ"]["PHOT_SELECTION"].dtype
    )

    return mydict


def get_protosteel_ra130d5_infos():
    """
    Get minimal photometric infos for protosteel ra130d5 (tertiary 0051, RA~130.5 DEC~0).

    Returns:
        mydict: dictionary with {keys: arrays},
            with keys: TARGETID, TERTIARY_TARGET, PHOT_RA, PHOT_DEC, PHOT_SELECTION
    """
    return _get_protosteel_infos(
        "0051", "hsc_icmodelmag_22-24_fewcols_RA130d5_DEC000.fits.gz"
    )


def get_protosteel_ra140_infos():
    """
    Get minimal photometric infos for protosteel ra140 (tertiary 0052, RA~140 DEC~3).

    Returns:
        mydict: dictionary with {keys: arrays},
            with keys: TARGETID, TERTIARY_TARGET, PHOT_RA, PHOT_DEC, PHOT_SELECTION
    """
    return _get_protosteel_infos(
        "0052", "hsc_icmodelmag_22-24_fewcols_RA140_DEC003.fits.gz"
    )


def get_protosteel_cosmos_infos():
    """
    Get minimal photometric infos for protosteel cosmos (tertiary 0055, COSMOS field).

    Returns:
        mydict: dictionary with {keys: arrays},
            with keys: TARGETID, TERTIARY_TARGET, PHOT_RA, PHOT_DEC, PHOT_SELECTION
    """
    return _get_protosteel_infos(
        "0055", "hsc_icmodelmag_22-24_COSMOS_fewcols_withexisting.fits.gz"
    )


def get_protosteel_phot_infos(case, d, photdir=None):
    """
    Get the photometric information (OBJECT_ID, FILENAME) for a given protosteel case.

    Args:
        case: round of DESI observation (str)
        d: output of the get_spec_table() function
        photdir: unused for protosteel; paths are determined by get_phot_fns()

    Returns:
        objids: OBJECT_ID values from the fewcols catalog, row-matched to d (array of int)
        targfns: full paths to the large photometric catalog, row-matched to d (array of str)
    """
    objids = np.zeros(len(d), dtype=int)
    large_fn = get_protosteel_large_fn(case)
    targfns = np.zeros(len(d), dtype="S200")

    for band in ["GRIZ"]:

        ii_band = np.where(d[band])[0]
        fns = get_phot_fns("protosteel", case, band)
        log.info("{}\t{}\t{}\t{}".format(case, band, ii_band.size, fns))

        if fns is None:
            continue

        for fn in fns:

            sel = (d[band]) & (objids == 0)
            ii_band = np.where(sel)[0]
            log.info(
                "{}\t{}\t{}\t{}/{} targets not dealt with yet".format(
                    case, band, os.path.basename(fn), ii_band.size, d[band].sum()
                )
            )

            t = Table.read(fn)
            for key in t.colnames:
                t[key].name = key.upper()

            iid, iit, _, _, _ = match_coord(
                d["PHOT_RA"][ii_band],
                d["PHOT_DEC"][ii_band],
                t["RA"],
                t["DEC"],
                search_radius=1.0,
                verbose=True,
            )
            log.info(
                "{}\t{}\t{:04d}/{:04d}\t{}".format(
                    case, band, iid.size, ii_band.size, os.path.basename(fn)
                )
            )

            iid = ii_band[iid]
            objids[iid] = t["OBJECT_ID"][iit]
            targfns[iid] = large_fn

        assert ((d[band]) & (objids == 0)).sum() == 0

    return objids, targfns
