"""Modelled read depth for `krsim vcf-host-depth --depth_model sWGA|gDNA`.

Each host's depth along the genome is simulated from a model fitted to MalariaGEN Pf8 samples:
sWGA AR(8) and gDNA ARMA(2,2), both with Student t kicks, on log window depth (500 bp windows).
The models describe depth relative to a sample's own typical (median) depth, so:

    window depth = host's typical depth x the model's relative value for that window

Per host:
  1. its typical depth is drawn once from a log-normal (centre depth_mean, SD of log depth_spread;
     depth_spread = 0 gives every host depth_mean)
  2. a relative depth track along each chromosome: kicks e ~ mu + scale * t(nu), filtered by the
     AR/MA model (y_t = sum ar_i y_t-i + e_t + sum ma_i e_t-i), then exp(y + offset) - pseudo
     (undoing log(depth / median + pseudo) minus its mean, as in the fits), at least 0; windows within
     end_width bp of a chromosome end multiplied by end_factor (sWGA: few reads at chromosome ends);
     then scaled so the host's median window is exactly 1, as a real sample's median counts its ends
  3. capped at `cap` times the typical depth
The read depth at each site is then Poisson(window depth).

Defaults: src/krsim/data/depth_model_defaults.json, made from the Pf8 fits by
real_data/pf8/model_comparison/make_depth_model_defaults.py.
"""

import json
import os
from importlib import resources

import numpy as np
from scipy.signal import lfilter

from krsim.output.pf8_profiles import CHROMS, CHROM_LENGTHS, WINDOW_SIZE

MODEL_TYPES = ("sWGA", "gDNA")
NUMBER_KEYS = ("mu", "scale", "nu", "offset", "cap", "end_width", "end_factor", "depth_mean", "depth_spread")
LIST_KEYS = ("ar", "ma")
SETTINGS = LIST_KEYS + NUMBER_KEYS
BURN_IN = 2000  # windows simulated and dropped before each chromosome starts, so the track starts settled


def load_defaults():
    """The defaults file: {"sWGA": {...}, "gDNA": {...}, "pseudo": ..., ...}."""
    with resources.files("krsim.data").joinpath("depth_model_defaults.json").open() as f:
        return json.load(f)


def model_settings(depth_model, depth_mean=None, depth_spread=None, params_file=None):
    """
    The settings for one run: the defaults for this sample type, then any values in params_file
    (a JSON object with any of SETTINGS), then --depth_mean / --depth_spread. Checked before use.
    """
    if depth_model not in MODEL_TYPES:
        raise ValueError(f"--depth_model must be one of {', '.join(MODEL_TYPES)}, got {depth_model}")
    defaults = load_defaults()
    s = {k: defaults[depth_model][k] for k in SETTINGS}
    s["model"], s["pseudo"] = defaults[depth_model]["model"], defaults["pseudo"]

    if params_file is not None:
        try:
            with open(params_file) as f:
                given = json.load(f)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{params_file} is not valid JSON ({exc})") from exc
        if not isinstance(given, dict):
            raise ValueError(f"{params_file} must hold one JSON object, e.g. {{\"nu\": 3.0, \"cap\": 1000}}")
        unknown = sorted(set(given) - set(SETTINGS))
        if unknown:
            raise ValueError(f"{params_file}: unknown setting{'s' if len(unknown) > 1 else ''} "
                             f"{', '.join(unknown)} (can set: {', '.join(SETTINGS)})")
        s.update(given)
        s["model"] = f"{s['model']} (changed by {os.path.basename(params_file)})"
    if depth_mean is not None:
        s["depth_mean"] = depth_mean
    if depth_spread is not None:
        s["depth_spread"] = depth_spread

    # check the values
    for k in LIST_KEYS:
        if not isinstance(s[k], list) or not all(isinstance(v, (int, float)) and not isinstance(v, bool)
                                                 for v in s[k]):
            raise ValueError(f"{k} must be a list of numbers, e.g. [1.1, -0.3]")
        s[k] = [float(v) for v in s[k]]
    for k in NUMBER_KEYS:
        if not isinstance(s[k], (int, float)) or isinstance(s[k], bool) or not np.isfinite(s[k]):
            raise ValueError(f"{k} must be a number, got {s[k]!r}")
        s[k] = float(s[k])
    for k in ("scale", "nu", "cap", "depth_mean"):
        if s[k] <= 0:
            raise ValueError(f"{k} must be above 0, got {s[k]}")
    for k in ("end_width", "end_factor", "depth_spread"):
        if s[k] < 0:
            raise ValueError(f"{k} must be 0 or more, got {s[k]}")
    if s["ar"]:
        # a usable AR part keeps depth from drifting off: every root of 1 - ar1 z - ar2 z^2 ... lies outside |z| = 1
        roots = np.roots(np.r_[-np.array(s["ar"])[::-1], 1.0])
        if np.any(np.abs(roots) <= 1.0):
            raise ValueError(f"ar = {s['ar']} would make depth drift without bound (the AR part is not "
                             "stationary); use smaller coefficients")
    return s


def host_depths(rng, s, n_hosts):
    """Each host's typical (median) depth, drawn once from the log-normal."""
    if s["depth_spread"] == 0:
        return np.full(n_hosts, s["depth_mean"])
    return s["depth_mean"] * np.exp(rng.normal(0.0, s["depth_spread"], n_hosts))


def relative_track(rng, s):
    """{3D7 chromosome: relative depth per 500 bp window}: ends lowered, median 1 over the genome, capped."""
    ar_poly = np.r_[1.0, -np.array(s["ar"])]
    ma_poly = np.r_[1.0, np.array(s["ma"])]
    track = {}
    for chrom in CHROMS:
        n = -(-CHROM_LENGTHS[chrom] // WINDOW_SIZE)
        kicks = s["mu"] + s["scale"] * rng.standard_t(s["nu"], n + BURN_IN)
        y = lfilter(ma_poly, ar_poly, kicks)[BURN_IN:]
        # cap y before exp so a rare huge run cannot overflow; the real cap comes below
        t = np.maximum(np.exp(np.minimum(y + s["offset"], 50.0)) - s["pseudo"], 0.0)
        # ends first, so the host's median (below) counts the low ends, as a real sample's median does
        end = min(int(round(s["end_width"] / WINDOW_SIZE)), len(t))
        if end > 0:
            t[:end] *= s["end_factor"]
            t[len(t) - end:] *= s["end_factor"]
        track[chrom] = t

    median = np.median(np.concatenate(list(track.values())))
    if not median > 0:
        median = np.mean(np.concatenate(list(track.values()))) or 1.0
    return {chrom: np.minimum(t / median, s["cap"]) for chrom, t in track.items()}


def describe(s):
    """One line for the VCF header: every value used."""
    def num(v):
        return f"{v:.6g}"
    parts = [f'model="{s["model"]}"', f'ar="{",".join(num(v) for v in s["ar"])}"',
             f'ma="{",".join(num(v) for v in s["ma"])}"']
    parts += [f"{k}={num(s[k])}" for k in NUMBER_KEYS] + [f"pseudo={num(s['pseudo'])}"]
    return ",".join(parts)
