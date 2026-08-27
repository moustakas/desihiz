#
# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-
"""
desihiz.fastspecfit_io
=======================

Run FastSpecFit on the merged desihiz catalogs (odin, suprime, clauds,
protosteel). Generalizes the original protosteel-only wrapper to any
img in ``allowed_imgs``, reusing the case/specprod bookkeeping already in
``hizmerge_io``.

Merged catalogs and all FastSpecFit outputs live under a single per-img
tree: ``$DESIHIZ_DIR/{img}/...``.

"""
import os
from glob import glob

import fitsio
import numpy as np
from astropy.table import Table

from desiutil.log import get_logger

from desihiz.hizmerge_io import allowed_imgs, get_img_cases, get_specprod, get_expids

log = get_logger()


def get_projectdir():
    """Base directory holding every img's merged catalogs and FastSpecFit outputs."""
    desihiz_dir = os.getenv("DESIHIZ_DIR")
    if desihiz_dir is not None:
        return desihiz_dir
    return os.path.join(os.getenv("DESI_ROOT"), "users", "ioannis", "desihiz")


def get_fphotofile(img):
    """Path to the packaged fastspecfit photometric-parameter yaml for img."""
    assert img in allowed_imgs
    return os.path.join(os.path.dirname(__file__), "data", "{}-photinfo.yaml".format(img))


def get_catalog(img, version=None):
    """Resolve the merged catalog path for img, defaulting to the latest version.

    Parameters
    ----------
    img : str
        Element of allowed_imgs.
    version : str or None
        Catalog version subdirectory (e.g. ``'v20260616'``). If None, the
        latest ``v*`` subdirectory under ``$DESIHIZ_DIR/{img}`` is used.

    Returns
    -------
    catalog : str
        Full path to the merged ``desi-{img}.fits`` catalog.
    """
    assert img in allowed_imgs
    projectdir = get_projectdir()

    if version is None:
        verdirs = sorted(glob(os.path.join(projectdir, img, "v*")))
        if not verdirs:
            raise FileNotFoundError(
                "No version directories found under {}".format(os.path.join(projectdir, img))
            )
        version = os.path.basename(verdirs[-1])

    return os.path.join(projectdir, img, version, "desi-{}.fits".format(img))


def get_photext(catalog):
    """Pick PHOTV2INFO if present in the catalog, else PHOTINFO."""
    extnames = [h.get_extname() for h in fitsio.FITS(catalog)]
    if "PHOTV2INFO" in extnames:
        return "PHOTV2INFO"
    return "PHOTINFO"


def _parse_fitsext(filepath):
    """Parse an optional FITS extension name from a 'path[EXTNAME]' string."""
    ext = None
    dirpath = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    if "[" in basename and "]" in basename:
        try:
            ext = basename[basename.find("[") + 1 : basename.find("]")]
            filepath = os.path.join(dirpath, basename[: basename.find("[")])
        except Exception:
            pass
    return filepath, ext


def _specprod_dir(coaddfile):
    """Extract the specprod root (ending in /healpix or /tiles) from a coadd path."""
    for tag in ("/healpix/", "/tiles/"):
        if tag in coaddfile:
            return coaddfile[: coaddfile.index(tag) + len(tag) - 1]
    raise ValueError("Cannot determine coadd type from path: {}".format(coaddfile))


def _outdir(img, prefix, coadd_type=None):
    """Return the top-level output directory for a given img and run prefix.

    Healpix paths need an explicit 'healpix' subdirectory because the path
    remainder after specprod_dir starts with the survey name, not 'healpix'.
    Cumulative paths already carry 'cumulative/' from the tiles directory
    structure, so no extra subdirectory is added.
    """
    base = os.path.join(get_projectdir(), img, prefix)
    if coadd_type == "healpix":
        return os.path.join(base, "healpix")
    return base


def _outfile(img, redrockfile, prefix, specprod_dir, coadd_type=None):
    """Derive the FastSpecFit output path from a redrock file path.

    Mirrors the convention in fastspecfit.mpi.plan:
      outfile = redrockfile.replace(specprod_dir, outdir).replace('redrock-', f'{prefix}-')
    """
    outdir = _outdir(img, prefix, coadd_type=coadd_type)
    fastfile = redrockfile.replace(specprod_dir, outdir).replace("redrock-", "{}-".format(prefix))
    if not prefix.startswith("fastphot"):
        fastfile += ".gz"
    return fastfile


def _get_cumulative_redrockfiles(img):
    """Return (redrockfile, tiles_dir) pairs across every tile-based case of img.

    Loops over img's cases, resolving each case's specprod via
    hizmerge_io.get_specprod() and its tiles via hizmerge_io.get_expids(),
    so this works for any img/case that has a standard tiles/cumulative
    tree -- today that is protosteel only, but nothing here is
    protosteel-specific. For a healpix-only img (e.g. odin/suprime/clauds,
    specprod "loa"), this simply finds no tiles directory and returns
    nothing.
    """
    from fastspecfit.mpi import findfiles

    redux_dir = os.path.join(os.getenv("DESI_ROOT"), "spectro", "redux")
    result = []
    seen_specprods = set()
    for case in get_img_cases(img):
        specprod = get_specprod(case)
        if specprod in seen_specprods:
            continue
        seen_specprods.add(specprod)

        tiles_dir = os.path.join(redux_dir, specprod, "tiles")
        if not os.path.isdir(tiles_dir):
            continue

        exp = get_expids(img, case)
        tiles = np.unique(exp["TILEID"])
        log.info("{}: {} unique tiles for case={}".format(specprod, len(tiles), case))
        redrockfiles = findfiles(tiles_dir, prefix="redrock", coadd_type="cumulative", tile=tiles)
        log.info("{}: {} cumulative redrock files found".format(specprod, len(redrockfiles)))
        for rr in redrockfiles:
            result.append((rr, tiles_dir))

    return result


def read_catalog(img, version=None, targetids=None, vi_redshifts=False, vi_quality_cut=2.0):
    """Read the img photometric and spectroscopic catalogs.

    Parameters
    ----------
    img : str
        Element of allowed_imgs.
    version : str or None
        Catalog version; see get_catalog().
    targetids : str or None
        Comma-separated list of TARGETIDs to keep, or None for all.
    vi_redshifts : bool
        If True, restrict to objects with valid visual-inspection redshifts.
    vi_quality_cut : float
        Minimum VI_QUALITY score when ``vi_redshifts`` is True.

    Returns
    -------
    fphoto, fspec : :class:`astropy.table.Table`
        Photometric and spectroscopic catalog tables, row-matched on TARGETID.
    """
    catalog = get_catalog(img, version=version)
    photext = get_photext(catalog)

    specpath, specext = _parse_fitsext("{}[SPECINFO]".format(catalog))
    photpath, photoext = _parse_fitsext("{}[{}]".format(catalog, photext))

    fphoto = Table(fitsio.read(photpath, ext=photoext))
    log.info("Read {:,d} objects from {}".format(len(fphoto), photpath))

    fspec = Table(fitsio.read(specpath, ext=specext))
    log.info("Read {:,d} objects from {}".format(len(fspec), specpath))
    assert np.all(fspec["TARGETID"] == fphoto["TARGETID"])

    if targetids is not None:
        I = np.isin(fspec["TARGETID"], np.array(targetids.split(","), dtype=np.int64))
        log.info("Trimming to {:,d} specified TARGETIDs.".format(I.sum()))
        fphoto = fphoto[I]
        fspec = fspec[I]

    if vi_redshifts:
        I = np.where(
            np.isfinite(fspec["VI_Z"])
            & (fspec["VI_Z"] > 1e-3)
            & (fspec["VI_QUALITY"] >= vi_quality_cut)
        )[0]
        log.info(
            "Trimming to {:,d} objects with VI_Z>0.001 and VI_QUALITY>={}.".format(
                len(I), vi_quality_cut
            )
        )
        fphoto = fphoto[I]
        fspec = fspec[I]

    return fphoto, fspec


def run_fastspec(
    img,
    version=None,
    targetids=None,
    fastphot=False,
    vi_redshifts=False,
    coadd_type="healpix",
    mp=1,
    ntargets=None,
    nolog=False,
    overwrite=False,
):
    """Run FastSpecFit on each img coadd file.

    For healpix coadds, iterates over unique coadd files in the catalog's
    COADDFN column. For cumulative tile coadds, discovers redrock files via
    _get_cumulative_redrockfiles(). Output files mirror the specprod
    directory structure under the project directory.

    Parameters
    ----------
    img : str
        Element of allowed_imgs.
    version : str or None
        Catalog version; see get_catalog().
    targetids : str or None
        Comma-separated TARGETIDs to process, or None for all.
    fastphot : bool
        If True, run fastphot (photometry only) instead of fastspec.
    vi_redshifts : bool
        If True, refit at visual-inspection redshifts.
    coadd_type : str
        Either ``'healpix'`` (default) or ``'cumulative'``.
    mp : int
        Number of multiprocessing workers.
    ntargets : int or None
        Limit the number of targets per coadd (useful for testing).
    nolog : bool
        If True, write output to stdout rather than a per-coadd log file.
    overwrite : bool
        If True, overwrite existing output files.
    """
    from desispec.parallel import stdouterr_redirected

    assert img in allowed_imgs

    prefix = "fastphot" if fastphot else "fastspec"
    if vi_redshifts:
        prefix += "-vi"

    catalog = get_catalog(img, version=version)
    fphotodir = "{}[{}]".format(catalog, get_photext(catalog))
    fphotofile = get_fphotofile(img)

    _, fspec = read_catalog(img, version=version, targetids=targetids, vi_redshifts=vi_redshifts)

    if coadd_type == "healpix":
        iteration = []
        for orig_coaddfile in sorted(set(fspec["COADDFN"])):
            I = fspec["COADDFN"] == orig_coaddfile
            redrockfile = os.path.join(
                os.path.dirname(orig_coaddfile),
                os.path.basename(orig_coaddfile).replace("coadd-", "redrock-"),
            )
            specprod_dir = _specprod_dir(orig_coaddfile)
            iteration.append((redrockfile, specprod_dir, I))
    else:
        iteration = []
        for redrockfile, specprod_dir in _get_cumulative_redrockfiles(img):
            I = None
            if vi_redshifts and os.path.isfile(redrockfile):
                # Scope the VI subset to just the TARGETIDs in this tile's
                # redrock file -- fspec carries no TILEID to do this via a
                # catalog-only cut, unlike the healpix branch above.
                rr_targetids = fitsio.read(redrockfile, "REDSHIFTS", columns=["TARGETID"])["TARGETID"]
                I = np.isin(fspec["TARGETID"], rr_targetids)
            iteration.append((redrockfile, specprod_dir, I))

    for redrockfile, specprod_dir, I in iteration:
        if not os.path.isfile(redrockfile):
            log.warning("Redrock file not found: {}".format(redrockfile))
            continue

        if vi_redshifts:
            vi_rows = fspec[I] if I is not None else fspec
            if len(vi_rows) == 0:
                log.info("No VI-redshift targets in {}; skipping.".format(redrockfile))
                continue

        fastfile = _outfile(img, redrockfile, prefix, specprod_dir, coadd_type=coadd_type)
        os.makedirs(os.path.dirname(fastfile), exist_ok=True)

        if os.path.isfile(fastfile) and not overwrite:
            log.info("Output file {} exists; use --overwrite to overwrite.".format(fastfile))
            continue

        cmdargs = (
            "{} -o {} --mp {}"
            " --fphotodir={} --fphotofile={}"
            " --ignore-quasarnet"
        ).format(redrockfile, fastfile, mp, fphotodir, fphotofile)
        if ntargets:
            cmdargs += " --ntargets {}".format(ntargets)
        if vi_redshifts:
            targetids_str = ",".join(vi_rows["TARGETID"].astype(str))
            input_redshifts = ",".join(vi_rows["VI_Z"].astype(str))
            cmdargs += " --targetids {} --input-redshifts {}".format(targetids_str, input_redshifts)

        if fastphot:
            from fastspecfit.fastspecfit import fastphot as fast
            log.info("fastphot {}".format(cmdargs))
        else:
            from fastspecfit.fastspecfit import fastspec as fast
            log.info("fastspec {}".format(cmdargs))

        if nolog:
            fast(args=cmdargs.split())
        else:
            logfile = fastfile.replace(".fits.gz", ".log").replace(".fits", ".log")
            with stdouterr_redirected(to=logfile, overwrite=True):
                fast(args=cmdargs.split())


def merge_fastspec(
    img, version=None, mp=1, vi_redshifts=False, fastphot=False, coadd_type="healpix", overwrite=False
):
    """Merge all per-coadd FastSpecFit outputs into a single catalog.

    For healpix, derives the expected output file list from the catalog's
    COADDFN column and warns about any missing files. For cumulative tiles,
    enumerates files via the same redrock-file list used by run_fastspec.

    Parameters
    ----------
    img : str
        Element of allowed_imgs.
    version : str or None
        Catalog version; see get_catalog().
    mp : int
        Number of multiprocessing workers passed to the merge step.
    vi_redshifts : bool
        If True, merge the VI-redshift refit outputs.
    fastphot : bool
        If True, merge fastphot outputs instead of fastspec.
    coadd_type : str
        Either ``'healpix'`` (default) or ``'cumulative'``.
    overwrite : bool
        If True, overwrite an existing merged catalog.
    """
    from fastspecfit.mpi import _domerge

    assert img in allowed_imgs

    prefix = "fastphot" if fastphot else "fastspec"
    if vi_redshifts:
        prefix += "-vi"

    mergefile = os.path.join(get_projectdir(), img, "{}-{}-{}.fits".format(prefix, img, coadd_type))

    if os.path.isfile(mergefile) and not overwrite:
        log.info("Output file {} exists; use --overwrite to overwrite.".format(mergefile))
        return

    fastfiles = []
    if coadd_type == "healpix":
        _, fspec = read_catalog(img, version=version, vi_redshifts=vi_redshifts)
        for orig_coaddfile in sorted(set(fspec["COADDFN"])):
            redrockfile = os.path.join(
                os.path.dirname(orig_coaddfile),
                os.path.basename(orig_coaddfile).replace("coadd-", "redrock-"),
            )
            specprod_dir = _specprod_dir(orig_coaddfile)
            fastfile = _outfile(img, redrockfile, prefix, specprod_dir, coadd_type=coadd_type)
            if os.path.isfile(fastfile):
                fastfiles.append(fastfile)
            else:
                log.warning("FastSpec file not found, skipping: {}".format(fastfile))
    else:
        for redrockfile, tiles_dir in _get_cumulative_redrockfiles(img):
            fastfile = _outfile(img, redrockfile, prefix, tiles_dir, coadd_type=coadd_type)
            if os.path.isfile(fastfile):
                fastfiles.append(fastfile)
            else:
                log.warning("FastSpec file not found, skipping: {}".format(fastfile))

    if not fastfiles:
        log.warning("No {} files found; run --fastspec first.".format(prefix))
        return

    _domerge(fastfiles, fastphot=fastphot, mergefile=mergefile, mp=mp)


def fastspec_qa(
    img,
    version=None,
    fastphot=False,
    vi_redshifts=False,
    coadd_type="healpix",
    mp=1,
    ntargets=None,
    overwrite=False,
):
    """Build FastSpecFit QA figures for all img coadd outputs.

    Iterates over the same (redrockfile, fastfile) pairs as run_fastspec and
    runs fastspecfit-qa on each output file that exists.

    Parameters
    ----------
    img : str
        Element of allowed_imgs.
    version : str or None
        Catalog version; see get_catalog().
    fastphot : bool
        If True, generate QA for fastphot outputs instead of fastspec.
    vi_redshifts : bool
        If True, generate QA for the VI-redshift refit outputs.
    coadd_type : str
        Either ``'healpix'`` (default) or ``'cumulative'``.
    mp : int
        Number of multiprocessing workers.
    ntargets : int or None
        Limit the number of targets per coadd (useful for testing).
    overwrite : bool
        If True, overwrite existing QA figures.
    """
    catalog = get_catalog(img, version=version)
    fphotodir = "{}[{}]".format(catalog, get_photext(catalog))
    fphotofile = get_fphotofile(img)

    prefix = "fastphot" if fastphot else "fastspec"
    if vi_redshifts:
        prefix += "-vi"

    qadir = os.path.join(get_projectdir(), img, "qa-{}".format(prefix))

    if coadd_type == "healpix":
        _, fspec = read_catalog(img, version=version, vi_redshifts=vi_redshifts)
        pairs = []
        for orig_coaddfile in sorted(set(fspec["COADDFN"])):
            redrockfile = os.path.join(
                os.path.dirname(orig_coaddfile),
                os.path.basename(orig_coaddfile).replace("coadd-", "redrock-"),
            )
            specprod_dir = _specprod_dir(orig_coaddfile)
            fastfile = _outfile(img, redrockfile, prefix, specprod_dir, coadd_type=coadd_type)
            outdir = os.path.dirname(fastfile).replace(_outdir(img, prefix), qadir)
            pairs.append((redrockfile, fastfile, outdir))
    else:
        pairs = []
        for rr, td in _get_cumulative_redrockfiles(img):
            fastfile = _outfile(img, rr, prefix, td, coadd_type=coadd_type)
            outdir = os.path.dirname(fastfile).replace(_outdir(img, prefix), qadir)
            pairs.append((rr, fastfile, outdir))

    for redrockfile, fastfile, outdir in pairs:
        if not os.path.isfile(redrockfile):
            log.warning("Redrock file not found: {}".format(redrockfile))
            continue
        if not os.path.isfile(fastfile):
            log.warning("FastSpec file not found, skipping QA: {}".format(fastfile))
            continue

        from fastspecfit.qa import fastqa
        cmdargs = "{} -o {} --redrockfiles {} --mp {} --fphotodir {} --fphotofile {}".format(
            fastfile, outdir, redrockfile, mp, fphotodir, fphotofile
        )
        if ntargets:
            cmdargs += " --ntargets {}".format(ntargets)
        if overwrite:
            cmdargs += " --overwrite"
        log.info("fastqa {}".format(cmdargs))
        fastqa(args=cmdargs.split())


def build_qa_html(img, fastphot=False, vi_redshifts=False, coadd_type="healpix"):
    """Build an HTML navigation page for the QA PNG figures.

    Walks the QA directory, builds a JSON manifest of all PNG files, and
    writes a single self-contained ``index.html`` with client-side navigation
    and search (no server required).

    Parameters
    ----------
    img : str
        Element of allowed_imgs.
    fastphot : bool
        If True, build the page for fastphot QA outputs.
    vi_redshifts : bool
        If True, build the page for VI-redshift refit QA outputs.
    coadd_type : str
        Either ``'healpix'`` (default) or ``'cumulative'``.
    """
    import json
    from glob import glob as _glob

    prefix = "fastphot" if fastphot else "fastspec"
    if vi_redshifts:
        prefix += "-vi"

    qadir = os.path.join(get_projectdir(), img, "qa-{}".format(prefix))
    root = os.path.join(qadir, coadd_type)

    if not os.path.isdir(root):
        log.warning("QA directory not found: {}".format(root))
        return

    records = []
    if coadd_type == "healpix":
        for sv in sorted(os.scandir(root), key=lambda e: e.name):
            if not sv.is_dir():
                continue
            for prog in sorted(os.scandir(sv.path), key=lambda e: e.name):
                if not prog.is_dir():
                    continue
                for grp in sorted(os.scandir(prog.path), key=lambda e: e.name):
                    if not grp.is_dir():
                        continue
                    for hpx in sorted(os.scandir(grp.path), key=lambda e: e.name):
                        if not hpx.is_dir():
                            continue
                        for png in sorted(_glob(os.path.join(hpx.path, "*.png"))):
                            targetid = os.path.basename(png).replace(".png", "").split("-")[-1]
                            records.append(
                                {
                                    "survey": sv.name,
                                    "program": prog.name,
                                    "group": grp.name,
                                    "healpix": hpx.name,
                                    "targetid": targetid,
                                    "path": os.path.relpath(png, start=root),
                                }
                            )
    elif coadd_type == "cumulative":
        for tile in sorted(os.scandir(root), key=lambda e: e.name):
            if not tile.is_dir():
                continue
            for night in sorted(os.scandir(tile.path), key=lambda e: e.name):
                if not night.is_dir():
                    continue
                for png in sorted(_glob(os.path.join(night.path, "*.png"))):
                    targetid = os.path.basename(png).replace(".png", "").split("-")[-1]
                    records.append(
                        {
                            "survey": "",
                            "program": "",
                            "group": tile.name,
                            "healpix": night.name,
                            "targetid": targetid,
                            "path": os.path.relpath(png, start=root),
                        }
                    )

    log.info("Found {} QA PNGs under {}".format(len(records), root))
    if not records:
        log.warning("No QA files found; run --qa first.")
        return

    title = "{} QA · {}".format(img, prefix)
    manifest_json = json.dumps(records, separators=(",", ":"))

    if coadd_type == "healpix":
        labels_json = json.dumps(
            {
                "l1": "group",
                "l1p": "groups",
                "l1_title": "Healpix groups",
                "l2": "pixel",
                "l2p": "pixels",
                "l2_title": "Healpix pixels",
                "search_hint": "healpix or targetid",
            }
        )
    else:
        labels_json = json.dumps(
            {
                "l1": "tile",
                "l1p": "tiles",
                "l1_title": "Tiles",
                "l2": "thrunight",
                "l2p": "thrunights",
                "l2_title": "Thrunights",
                "search_hint": "tile or targetid",
            }
        )

    # Use placeholder substitution to avoid f-string brace escaping in CSS/JS.
    template = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
body{font-family:sans-serif;max-width:1000px;margin:0 auto;padding:1em;color:#333}
h1{font-size:1.15em;border-bottom:1px solid #ccc;padding-bottom:.3em}
#search{width:100%;box-sizing:border-box;padding:.5em;font-size:1em;margin-bottom:.75em}
#nav{margin-bottom:.6em;font-size:.9em;color:#555}
.row{padding:.4em .5em;border-bottom:1px solid #eee;cursor:pointer}
.row:hover{background:#f5f5f5}
.grid{display:flex;flex-wrap:wrap;gap:.5em;margin-top:.5em}
.thumb{width:180px;text-align:center}
.thumb img{width:180px;border:1px solid #ddd;display:block}
.tid{font-size:.65em;color:#666;word-break:break-all;margin-top:.2em}
a{color:#0066cc;text-decoration:none}
a:hover{text-decoration:underline}
h3{margin:.8em 0 .3em;font-size:.9em;color:#555;text-transform:uppercase;letter-spacing:.05em}
</style>
</head>
<body>
<h1>__TITLE__</h1>
<input type="text" id="search" placeholder="search __SEARCH_HINT__…" oninput="doSearch()">
<div id="nav"></div>
<div id="content"></div>
<script>
const DATA=__DATA__;
const L=__LABELS__;
const byGroup={};
DATA.forEach(r=>{
  if(!byGroup[r.group])byGroup[r.group]={};
  if(!byGroup[r.group][r.healpix])byGroup[r.group][r.healpix]=[];
  byGroup[r.group][r.healpix].push(r);
});
function breadcrumb(parts){return parts.filter(Boolean).join(' / ');}
function showRoot(){
  document.getElementById('search').value='';
  document.getElementById('nav').innerHTML='';
  const groups=Object.keys(byGroup).sort((a,b)=>+a-+b);
  let h='';
  groups.forEach(g=>{
    const pix=Object.keys(byGroup[g]);
    const n=pix.reduce((s,p)=>s+byGroup[g][p].length,0);
    h+=`<div class="row" onclick="showGroup('${g}')">${g} · ${pix.length} ${pix.length>1?L.l2p:L.l2} · ${n} spectra</div>`;
  });
  document.getElementById('content').innerHTML=h||'<p>No QA files found.</p>';
}
function showGroup(g){
  document.getElementById('search').value='';
  const r0=Object.values(byGroup[g])[0][0];
  document.getElementById('nav').innerHTML=
    `<a href="#" onclick="showRoot()">‹ ${L.l1p}</a>   ${breadcrumb([r0.survey,r0.program,g])}`;
  const pix=Object.keys(byGroup[g]).sort((a,b)=>+a-+b);
  let h='';
  pix.forEach(p=>{
    const n=byGroup[g][p].length;
    h+=`<div class="row" onclick="showPixel('${g}','${p}')">${p} · ${n} spectr${n>1?'a':'um'}</div>`;
  });
  document.getElementById('content').innerHTML=h;
}
function showPixel(g,p){
  document.getElementById('search').value='';
  const recs=byGroup[g][p],r0=recs[0];
  document.getElementById('nav').innerHTML=
    `<a href="#" onclick="showRoot()">‹ ${L.l1p}</a>  `+
    `<a href="#" onclick="showGroup('${g}')">‹ ${g}</a>  `+
    `${breadcrumb([r0.survey,r0.program,g,p])} · ${recs.length} spectra`;
  let h='<div class="grid">';
  recs.forEach(r=>{
    h+=`<div class="thumb"><a href="${r.path}" target="_blank"><img src="${r.path}" loading="lazy"><div class="tid">${r.targetid}</div></a></div>`;
  });
  document.getElementById('content').innerHTML=h+'</div>';
}
function doSearch(){
  const q=document.getElementById('search').value.trim();
  if(!q){showRoot();return;}
  document.getElementById('nav').innerHTML='';
  let h='';
  const mg=Object.keys(byGroup).filter(g=>g.includes(q)).sort((a,b)=>+a-+b);
  if(mg.length){
    h+=`<h3>${L.l1_title}</h3>`;
    mg.forEach(g=>{
      const pix=Object.keys(byGroup[g]);
      const n=pix.reduce((s,p)=>s+byGroup[g][p].length,0);
      h+=`<div class="row" onclick="showGroup('${g}')">${g} · ${pix.length} ${L.l2p} · ${n} spectra</div>`;
    });
  }
  const mp=[];
  Object.keys(byGroup).forEach(g=>Object.keys(byGroup[g]).forEach(p=>{if(p.includes(q))mp.push([g,p]);}));
  if(mp.length){
    h+=`<h3>${L.l2_title}</h3>`;
    mp.sort((a,b)=>+a[1]-+b[1]).forEach(([g,p])=>{
      const n=byGroup[g][p].length;
      h+=`<div class="row" onclick="showPixel('${g}','${p}')">${p} (${L.l1} ${g}) · ${n} spectra</div>`;
    });
  }
  const mt=DATA.filter(r=>r.targetid.includes(q));
  if(mt.length){
    h+=`<h3>Targets (${mt.length}${mt.length>50?' — showing first 50':''})</h3>`;
    mt.slice(0,50).forEach(r=>{
      h+=`<div class="row"><a href="${r.path}" target="_blank">${r.targetid}</a>   ${L.l1} ${r.group}</div>`;
    });
  }
  document.getElementById('content').innerHTML=h||`<p>No results for “${q}”.</p>`;
}
showRoot();
</script>
</body>
</html>"""

    html = (
        template.replace("__TITLE__", title)
        .replace("__SEARCH_HINT__", json.loads(labels_json)["search_hint"])
        .replace("__DATA__", manifest_json)
        .replace("__LABELS__", labels_json)
    )
    htmlfile = os.path.join(root, "index.html")
    with open(htmlfile, "w") as f:
        f.write(html)
    log.info("Wrote {}".format(htmlfile))
