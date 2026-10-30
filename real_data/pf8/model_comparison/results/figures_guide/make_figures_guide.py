"""One PDF with every depth-model figure and its explanation: figures_guide.pdf (open in Preview)."""
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
TEXT, MUTED = "#0b0b0b", "#52514e"

PAGES = [
    ("Kicks: the random surprise at each window", "1_kicks_diagram.png", [
        ("What's on it", [
            "Panel 1: each window's new value = a prediction from the past (blue circle) + a random kick (orange arrow).",
            "Panel 2: kicks are drawn from a Gaussian with mean 0 and spread σ. 68% fall within ±1σ, 95% within ±2σ, 99.7% within ±3σ.",
            "Panel 3: what one kick of size 1 does to the windows after it, under three models.",
        ]),
        ("Takeaways", [
            "The kick is new randomness, independent of the past. The prediction is what depends on the past. σ is fitted per sample.",
            "AR(1) remembers 1 previous depth: a kick fades gradually and never fully stops.",
            "MA(3) remembers 3 previous kicks: a kick pushes the next 3 windows, then stops.",
            "ARMA(1,3) does both: a bump for 3 windows, then a fade.",
        ]),
    ]),
    ("AR: how far back does the model need to look?", "2_ar_lookback.png", [
        ("What's on it", [
            "x-axis: how far back the model looks (kb), with the AR order underneath (AR(1) = 0.5 kb, AR(8) = 4 kb).",
            "y-axis: held-out score, i.e. how far behind each sample's best model (lower = better, log scale).",
            "Lines = median per sample type; faint dots = individual samples.",
        ]),
        ("Takeaways", [
            "sWGA (blue) falls steeply, from about 2,900 at AR(1) to about 15 at AR(8), then flattens. "
            "Looking back about 3-4 kb captures what matters; sWGA depth moves in stretches a few kb long.",
            "gDNA (orange) flattens at about 45-70 and never reaches 0. No AR does as well as gDNA's best model, ARMA(2,2), however far back it looks.",
        ]),
    ]),
    ("ARMA and Gaussian-process models", "3_arma_gp.png", [
        ("What's on it", [
            "Median held-out score per model, sWGA (blue) and gDNA (orange) side by side, with the number on each bar.",
            "Left: the three ARMA models. Right: the five GP kernels.",
        ]),
        ("Takeaways", [
            "ARMA, gDNA: ARMA(2,2) scores 0 (best for gDNA); ARMA(1,1) and (2,1) score about 125, so both look-back parts are needed.",
            "ARMA, sWGA: ARMA(2,2) (80) is behind the mid-length ARs (about 15).",
            "GP: all far behind. Squared exponential, the smoothest kernel, is the worst for both types: depth is rough, not smooth.",
            "Check: GP exponential (2,868) matches AR(1) (2,867), because they are the same model.",
        ]),
    ]),
    ("Every model in order, per sample type", "4_ranking.png", [
        ("What's on it", [
            "Models sorted from worst (bottom) to best (top) by median held-out score. Shorter bar = better.",
            "The coloured bar is the best model; '(n wins)' = best for n samples.",
        ]),
        ("Takeaways", [
            "gDNA: ARMA(2,2) is far ahead and wins 9 of 11 samples.",
            "sWGA: the top is a near-tie of AR(7)-AR(20); the wins are spread across several models, so there is no single sWGA winner.",
            "Both: AR(1), GP exponential and GP squared exponential are at the bottom.",
        ]),
    ]),
    ("Which model wins for each sample?", "5_ranks.png", [
        ("What's on it", [
            "Rows = samples (sWGA above the black line, gDNA below). Columns = models.",
            "Number = that model's rank for that sample (1 = best, 20 = worst); darker = better.",
        ]),
        ("Takeaways", [
            "gDNA: a solid dark first column, with ARMA(2,2) ranked 1 in 9 samples.",
            "The exceptions are the two TRAC2 samples (QC0637-C, RCN12993), whose best models are AR(5)-AR(6), like sWGA.",
            "sWGA: the dark cells are scattered across AR(5)-AR(20), with a different best for different samples.",
            "The right-hand columns (smooth GPs, AR(1)) are pale for every sample.",
        ]),
    ]),
    ("Every sample: real vs AR(1) vs the best model", "6_examples.png", [
        ("What's on it", [
            "One panel per sample (sWGA on top, gDNA below), ordered from AR(1) closest to the best model to furthest.",
            "Lines: how similar windows are by distance over the first 5 kb: real (black), AR(1) (grey dashed), the type's best model (colour).",
            "Each panel has its own y-axis, so small differences show. The title gives how far AR(1) is behind the sample's best model.",
        ]),
        ("Takeaways", [
            "sWGA: in most samples AR(1) misses badly. The real curve bends (a different rate of fading at short and longer distances), "
            "which a single exponential cannot follow; AR(15) tracks it almost exactly. Only RCN02790 is close for AR(1).",
            "gDNA: ARMA(2,2) follows the real curve in every sample. AR(1) looks close by eye in typical gDNA (gaps 258-646), "
            "but is clearly off for the two TRAC2 samples (2,898 and 3,499).",
        ]),
    ]),
    ("gDNA in detail: what ARMA(2,2) says", "7_gdna_arma22.png", [
        ("What's on it", [
            "Left: the 9 typical gDNA samples (TRAC2 left out), sorted by depth. ARMA(2,2) split into a fast part (blue), "
            "a slow part (green) and a 'noise' part (grey, negative). Text: how long the fast and slow ups and downs last.",
            "Right: real similarity-by-distance curves (black) against the fitted ARMA(2,2) (orange).",
        ]),
        ("Takeaways", [
            "Almost all the variation is the fast part, lasting about 0.7 kb, in every sample.",
            "The slow part is small (about 5%) and lasts 15-35 kb.",
            "The 'noise' share is negative (about -20%): there is no independent per-window noise. Instead, neighbouring windows "
            "are extra alike (probably reads straddling window boundaries).",
            "The same shape from 10x to 177x: it does not depend on sequencing depth. The fit (orange) follows the real curves (black).",
        ]),
    ]),
    ("All scores: AIC, BIC and held-out", "8_scores.png", [
        ("What's on it", [
            "Every model's score in every sample (dots) and the median (bar), sWGA and gDNA side by side, under AIC, BIC and held-out.",
        ]),
        ("Takeaways", [
            "BIC and held-out mostly agree. AIC's penalty is too small for this much data: for sWGA it often picks the biggest model, AR(20).",
            "Held-out is the main score; BIC is the backup.",
        ]),
    ]),
]


def page(pdf, title, image, sections):
    fig = plt.figure(figsize=(11.69, 8.27))   # A4 landscape
    fig.text(0.04, 0.955, title, fontsize=17, color=TEXT, weight="bold", va="top")
    img = mpimg.imread(os.path.join(FIG, image))
    h, w = img.shape[:2]
    box_w, box_h = 0.92, 0.56                  # space for the image, as a fraction of the page
    page_ratio = 8.27 / 11.69
    scale = min(box_w, box_h * page_ratio * w / h) / box_w
    iw = box_w * scale
    ih = iw * h / w / page_ratio
    ax = fig.add_axes([0.04 + (box_w - iw) / 2, 0.91 - ih, iw, ih])
    ax.imshow(img)
    ax.axis("off")
    y = 0.91 - ih - 0.035
    for heading, bullets in sections:
        fig.text(0.04, y, heading, fontsize=11, color=TEXT, weight="bold", va="top")
        y -= 0.03
        for b in bullets:
            lines = textwrap.wrap(b, 165)
            fig.text(0.05, y, "•  " + "\n    ".join(lines), fontsize=9, color=MUTED, va="top", linespacing=1.35)
            y -= 0.024 * len(lines) + 0.008
        y -= 0.01
    pdf.savefig(fig)
    plt.close(fig)


with PdfPages(os.path.join(HERE, "figures_guide.pdf")) as pdf:
    for title, image, sections in PAGES:
        page(pdf, title, image, sections)
print("wrote", os.path.join(HERE, "figures_guide.pdf"))
