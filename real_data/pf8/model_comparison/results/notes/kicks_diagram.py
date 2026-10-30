import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import lfilter

TEXT, MUTED, GRID = "#0b0b0b", "#8a8984", "#e6e5e1"
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#b9b8b3"


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9, labelcolor=TEXT)


fig = plt.figure(figsize=(14, 10.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1], hspace=0.42, wspace=0.22)

# --- 1: new value = prediction + kick, window by window (AR(1), a = 0.7) ---
ax = fig.add_subplot(gs[0, :])
a = 0.7
kicks = np.array([0.0, 0.5, -0.3, 0.8, -0.6, 0.2, -0.4, 0.6, 0.1, -0.5, 0.3])
y = np.zeros(len(kicks))
pred = np.zeros(len(kicks))
y[0] = 0.4
for t in range(1, len(kicks)):
    pred[t] = a * y[t - 1]
    y[t] = pred[t] + kicks[t]
w = np.arange(len(kicks))
ax.axhline(0, color=GRID, lw=1)
ax.plot(w, y, color=TEXT, lw=1.2, zorder=1)
for t in range(1, len(kicks)):
    ax.scatter(t, pred[t], s=45, facecolor="white", edgecolor=BLUE, lw=1.8, zorder=3)
    ax.annotate("", xy=(t, y[t]), xytext=(t, pred[t]),
                arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2, shrinkA=4, shrinkB=4), zorder=2)
    ax.annotate("", xy=(t - 0.1, pred[t]), xytext=(t - 0.9, y[t - 1]),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1, ls="--", shrinkA=5, shrinkB=5, alpha=0.6))
ax.scatter(w, y, s=55, color=TEXT, zorder=4)
t = 3
ax.annotate("prediction\n(0.7 × previous value)", xy=(t, pred[t]), xytext=(t + 0.35, pred[t] - 0.45),
            fontsize=9, color=BLUE, arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.8))
ax.annotate(f"kick = +{kicks[t]:.1f}\n(new random surprise)", xy=(t, (pred[t] + y[t]) / 2), xytext=(t + 0.35, 1.25),
            fontsize=9, color=ORANGE, arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8))
ax.annotate("new value\n= prediction + kick", xy=(t, y[t]), xytext=(t - 1.6, 1.25),
            fontsize=9, color=TEXT, arrowprops=dict(arrowstyle="-", color=TEXT, lw=0.8))
ax.set_xticks(w, [f"window {i}" for i in w], fontsize=8)
ax.set_ylabel("log depth", fontsize=10)
ax.set_ylim(-0.75, 1.55)
ax.set_title("1. Every window: new value = prediction from the past (blue ○) + a random kick (orange arrow)",
             loc="left", fontsize=11, color=TEXT)
style(ax)

# --- 2: the kicks are drawn from a Gaussian ---
ax = fig.add_subplot(gs[1, 0])
x = np.linspace(-4, 4, 400)
pdf = np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi)
ax.plot(x, pdf, color=ORANGE, lw=2)
for k, alpha, label in ((3, 0.12, "99.7%"), (2, 0.2, "95%"), (1, 0.35, "68%")):
    m = np.abs(x) <= k
    ax.fill_between(x[m], pdf[m], color=ORANGE, alpha=alpha, lw=0)
for k, label, yy in ((1, "68% of kicks\nwithin ±1σ", 0.2), (2, "95% within ±2σ", 0.06), (3, "99.7% within ±3σ", 0.012)):
    ax.annotate(label, xy=(k, yy), xytext=(k + 0.15, yy + 0.02), fontsize=8, color=TEXT)
ax.axvline(0, color=MUTED, lw=0.8, ls=":")
ax.set_xticks(range(-3, 4), ["−3σ", "−2σ", "−1σ", "0", "+1σ", "+2σ", "+3σ"])
ax.set_yticks([])
ax.set_xlabel("size of the kick (on the log-depth scale)", fontsize=10)
ax.set_title("2. Kicks are drawn from a Gaussian: mean 0, spread σ\n(same σ at every window; σ is fitted for each sample)",
             loc="left", fontsize=11, color=TEXT)
style(ax)

# --- 3: what one kick does over the following windows ---
ax = fig.add_subplot(gs[1, 1])
n = 10
imp = np.zeros(n); imp[0] = 1
ar1 = lfilter([1], [1, -0.5], imp)
ma3 = lfilter([1, 0.6, 0.3, 0.1], [1], imp)
arma13 = lfilter([1, 0.6, 0.3, 0.1], [1, -0.5], imp)
lags = np.arange(n)
ax.plot(lags, ar1, "o-", color=BLUE, lw=1.8, label="AR(1): remembers 1 previous depth → fades, never stops")
ax.plot(lags, ma3, "s-", color=AQUA, lw=1.8, label="MA(3): remembers 3 previous kicks → stops after 3")
ax.plot(lags, arma13, "^-", color=ORANGE, lw=1.8, label="ARMA(1,3): both → a bump, then a fade")
ax.axhline(0, color=GRID, lw=1)
ax.set_xticks(lags, ["kick"] + [f"+{i}" for i in lags[1:]])
ax.set_xlabel("windows after the kick", fontsize=10)
ax.set_ylabel("effect on depth", fontsize=10)
ax.set_title("3. What one kick of 1 does to the windows after it\n(example weights: a = 0.5; b = 0.6, 0.3, 0.1)",
             loc="left", fontsize=11, color=TEXT)
ax.legend(frameon=False, fontsize=8.5, loc="upper right")
style(ax)

fig.suptitle("Kicks: the random surprise added at each window", x=0.01, ha="left", fontsize=14, color=TEXT)
fig.savefig("/Users/bty208/Desktop/RA_POSTDOC_IMPERIAL/pf8_samples/model/notes/kicks_diagram.png",
            dpi=160, facecolor="white", bbox_inches="tight")
