# Notes: custom high-z Redrock templates

Status: idea/assessment stage, nothing implemented yet. Written 2026-08-26 after
reviewing `~/code/fastspecfit-projects/laelbg-templates/laelbg-templates`
(stale prototype -- old fastspecfit I/O, not reusable as-is) and discussing
strategy. Picks up from the `fastspecfit_io.py`/`redrock_io.py` refactor
already merged into this repo.

## The problem

Pipeline Redrock redshifts for the odin/suprime/clauds high-z (LAE) samples
are unreliable at z~2-4: the default GALAXY template doesn't give Redrock a
genuinely competitive z~2-4 solution to fit against, so lower-z interlopers
([OII]/Lya confusion etc.) can win the chi2 competition even when the true
redshift is high. This is a *template* problem, not a S/N or statistics
problem -- `hizmerge_io.py` already carries `RR_Z0`/`RR_Z1` dual-redshift
bookkeeping as a post-hoc workaround for exactly this.

Goal: build a Redrock-compatible high-z template (or archetype set) so the
redshift search itself can distinguish the correct high-z solution, using
the VI-confirmed odin/suprime/clauds samples as ground truth.

## Finding: `--constrain-age` already does part of what we want

FastSpecFit has a `--constrain-age` flag (confirmed in the checked-out
source, `fastspecfit/continuum.py:106`, `_younger_than_universe()`): it
restricts the fitted SPS template ages to be younger than the age of the
universe at the object's redshift.

- The **old** `fastspecfit-desihiz` script used it.
- The **protosteel** script (basis for the new `fastspecfit_io.py`) never
  had it, so it was dropped in the refactor.
- For protosteel (low-z-ish weak-lensing calibration sample) this doesn't
  matter much either way.
- For odin/suprime/clauds at z=2-4 it matters a lot: without it, FastSpecFit
  can assign unphysically old stellar populations (e.g. 10 Gyr at z=3),
  visibly distorting the derived continuum shape.

**TODO**: re-add `--constrain-age` as a passthrough CLI flag in
`desi_hiz_fastspecfit` before generating any fits meant to seed template
construction.

Separate, still-open question: does FastSpecFit 3.6.1 expose a way to
restrict the underlying SPS *template grid* itself (metallicity, burst mode,
etc.), as opposed to this runtime per-fit age cut? Not yet checked.

## Reframe: the low-z side probably doesn't need rebuilding

The `laelbg-templates` prototype's `--fuji-templates` path builds a fresh
low-z (z<1.6) BGS/ELG/LRG template library from the official DESI SV3
fastspecfit catalog. This is likely unnecessary now:

- Redrock's *existing default* GALAXY template was already trained on a
  large, diverse low-z DESI sample -- the low-z regime already works.
- The millions of fastspecfit fits from the recently-completed Loa run
  don't change that calculus; statistics were never the bottleneck there.

**Recommendation**: don't rebuild the low-z basis from scratch. Instead,
add a new high-z-specific template/subtype *alongside* the existing default
GALAXY template, and let Redrock's native multi-template chi2 scan do the
disambiguation. Revisit only if the augmented approach turns out to need a
refreshed low-z draw for some reason (e.g. joint-fit basis choice below).

## Basis-construction method: PCA vs. NMF vs. archetypes

Framing: is this "compress a smooth continuum of galaxy diversity" (favors
PCA/NMF) or "discriminate between two structurally different discrete
populations" (favors archetypes)? Low-z "normal" galaxies and high-z
Lya-emitting galaxies are not points on one smooth manifold -- they're
qualitatively different SEDs. That points toward archetypes.

- **PCA**: what Redrock's native linear-template mode was designed around.
  Cheapest, most battle-tested path. Eigenspectra can go negative (fine
  mathematically -- Redrock's fit-time linear coefficients aren't
  constrained in sign regardless of how the basis was built).
- **NMF**: non-negative, "parts-based," more interpretable components, but
  no compression-optimality guarantee over PCA (needs more components for
  equivalent reconstruction fidelity). Its main selling point doesn't
  translate into better redshift recovery, since Redrock's fit-time
  coefficients aren't sign-constrained either way.
- **Archetypes** (the CPLEX idea): matches Redrock's native archetype-mode
  fitting, which cross-correlates against discrete real/model spectra
  directly rather than a compressed linear basis. Good match for a
  genuinely bimodal population, and for a sample (VI-confirmed high-z) that
  is statistics-limited rather than statistics-rich -- no pressure to
  compress lossily when you want to preserve distinguishing UV/Lya
  features rather than average them away.

**Lean**: use archetype mode for the high-z addition, built directly from
VI-confirmed, age-constrained FastSpecFit rest-frame models. Leave the
existing default low-z template untouched initially.

`NonnegMFPy` (the NMF package the old prototype imports) is not currently a
dependency anywhere in this repo or environment -- confirmed not
importable. Only add it if NMF ends up being the chosen path.

## Proposed validation

Rerun Redrock via `desi_hiz_redrock --template-dir` twice -- once with the
stock template set, once with stock + new high-z archetypes -- and compare
recovered z against `VI_Z` for the VI-confirmed subset. Track:
- purity / catastrophic-failure rate
- completeness of correct high-z recovery vs. low-z aliasing

## Open questions / things to check before committing to a path

1. Does DESI/Anand's group already have archetype-mode infrastructure for
   large-scale-structure redshift refinement that could be reused instead
   of building this from scratch? (Ask before searching -- don't assume.)
2. What is FastSpecFit 3.6.1's actual API (if any) for restricting the
   underlying SPS template grid itself, beyond `--constrain-age`?
3. Science choices inherited from the old prototype that still need a
   decision if any low-z work does happen: age cut (z<1.6), `DELTACHI2>40`
   purity cut, per-class exemplar count (`ntemp=100`), NMF `ncomp=10` --
   none of these should be carried over silently.

## Code housekeeping (not urgent, noted for later)

- `restframe_templates()` in the old prototype is stale: it was written
  against an older FastSpecFit I/O layer (`read_fastspecfit`,
  `cache_templates()` usage) that has since changed. Needs a rewrite
  against current FastSpecFit I/O, not a straight port, whenever this
  becomes active work.
- The prototype script also has a live `pdb.set_trace()` and a dead
  `if False:` block left in from exploration -- drop both.
- Eventual home for any of this, following the `fastspecfit_io.py` /
  `redrock_io.py` pattern already in this repo: a new
  `py/desihiz/templates_io.py` + `bin/desi_hiz_templates`.
