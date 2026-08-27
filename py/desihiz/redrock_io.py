#
# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-
"""
desihiz.redrock_io
===================

Rerun Redrock on desihiz healpix coadds with a custom template directory
(e.g. to test new high-z templates), and gather the results into a single
row-matched catalog.

This is infrastructure only: it reruns Redrock and collects the outputs.
Building the custom templates themselves, and validating the resulting
redshifts against the VI subset, is separate downstream work.

"""
import os
from glob import glob

import fitsio
import numpy as np
from astropy.table import Table, vstack

from desiutil.log import get_logger

from desihiz.hizmerge_io import allowed_imgs
from desihiz.fastspecfit_io import get_projectdir, read_catalog

log = get_logger()


def run_redrock(img, version=None, targetids=None, template_dir=None, rrdetails=False, overwrite=False):
    """Rerun Redrock on img's healpix coadds.

    Copies the targets of interest out of the original coadd (so the
    original DESI reduction is untouched) into
    $DESIHIZ_DIR/{img}/redux/{subdir}/coadd-{healpix}.fits, then runs
    Redrock against it, optionally with a custom RR_TEMPLATE_DIR.

    Parameters
    ----------
    img : str
        Element of allowed_imgs.
    version : str or None
        Catalog version; see fastspecfit_io.get_catalog().
    targetids : str or None
        Comma-separated TARGETIDs to process, or None for all.
    template_dir : str or None
        If set, exported as RR_TEMPLATE_DIR before running Redrock.
    rrdetails : bool
        If True, also write the per-target rrdetails HDF5 file.
    overwrite : bool
        If True, overwrite existing coadd/redrock files.
    """
    from redrock.external.desi import rrdesi
    from desispec.io import read_spectra, write_spectra

    assert img in allowed_imgs

    if template_dir is not None:
        os.environ["RR_TEMPLATE_DIR"] = template_dir

    baseoutdir = os.path.join(get_projectdir(), img, "redux")
    os.makedirs(baseoutdir, exist_ok=True)

    _, fspec = read_catalog(img, version=version, targetids=targetids)
    for orig_coaddfile in sorted(set(fspec["COADDFN"])):
        I = orig_coaddfile == fspec["COADDFN"]
        healpix = fspec["HEALPIX"][I][0]
        alltargetids = fspec["TARGETID"][I].data

        subdir = orig_coaddfile.split("/")[-2]
        outdir = os.path.join(baseoutdir, subdir)
        os.makedirs(outdir, exist_ok=True)

        coaddfile = os.path.join(outdir, "coadd-{}.fits".format(healpix))
        if not os.path.isfile(coaddfile) or overwrite:
            spec = read_spectra(orig_coaddfile, targetids=alltargetids)
            assert np.all(spec.target_ids() == alltargetids)
            log.info("Writing {} targets to {}".format(len(alltargetids), coaddfile))
            write_spectra(coaddfile, spec)
            del spec

        redrockfile = os.path.join(outdir, "redrock-{}.fits".format(healpix))
        rrdetailsfile = os.path.join(outdir, "rrdetails-{}.h5".format(healpix))
        if not os.path.isfile(redrockfile) or overwrite:
            cmd = "-i {} -o {} --gpu --max-gpuprocs=4 --mp 1 --zscan-galaxy=-0.005,4.0,3e-4".format(
                coaddfile, redrockfile
            )
            if rrdetails:
                cmd += " -d {}".format(rrdetailsfile)
            log.info("rrdesi {}".format(cmd))
            rrdesi(cmd.split())
        else:
            log.info("Skipping existing file {}".format(redrockfile))


def gather_redrock(img, version=None, overwrite=False):
    """Merge the rerun Redrock results into a catalog row-matched to img's catalog.

    Parameters
    ----------
    img : str
        Element of allowed_imgs.
    version : str or None
        Catalog version; see fastspecfit_io.get_catalog().
    overwrite : bool
        If True, overwrite an existing merged catalog.
    """
    from desitarget import geomask

    assert img in allowed_imgs

    projectdir = get_projectdir()
    redrockdir = os.path.join(projectdir, img, "redux")
    _, fspec = read_catalog(img, version=version)

    rrmergefile = os.path.join(projectdir, img, "redrock-{}.fits".format(img))
    if os.path.isfile(rrmergefile) and not overwrite:
        log.info("Output file {} exists; use --overwrite to overwrite.".format(rrmergefile))
        return

    rrfiles = glob(os.path.join(redrockdir, "*", "redrock-?????.fits"))
    if not rrfiles:
        log.warning("No redrock files found under {}; run --run-redrock first.".format(redrockdir))
        return

    zcat = vstack([Table(fitsio.read(rrfile, "REDSHIFTS")) for rrfile in rrfiles])
    zcat = zcat[geomask.match_to(zcat["TARGETID"], fspec["TARGETID"])]
    assert np.all(zcat["TARGETID"] == fspec["TARGETID"])

    zcat.write(rrmergefile, overwrite=True)
    log.info("Wrote {:,d} objects to {}".format(len(zcat), rrmergefile))
