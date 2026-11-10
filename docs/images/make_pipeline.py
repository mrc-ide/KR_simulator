"""
Draws docs/images/pipeline.png, the krsim pipeline figure on the website's home page.

    python docs/images/make_pipeline.py

Grey boxes are krsim commands, pink pills the files each writes (inside its dashed outline), and arrows
show which files each command reads. Coordinates are in pixels of a 2,100 x 1,640 canvas, y downwards.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

W, H = 2100, 1640
INK, RED, GREY, PINK = "#0f2333", "#e03131", "#e8e8e8", "#f0cdeb"
FONT = {"family": "DejaVu Sans Mono", "weight": "bold"}

fig = plt.figure(figsize=(W / 100, H / 100), dpi=140)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")


def group(x0, y0, x1, y1):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=INK, lw=2.2, ls=(0, (5, 4)), zorder=1))


def command(x0, y0, x1, y1, text):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=GREY, ec=INK, lw=3, zorder=3))
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, text, color=RED, fontsize=30, ha="center", va="center", zorder=4, **FONT)


def pill(x0, y0, x1, y1, text):
    h = y1 - y0
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, h, boxstyle=f"round,pad=0,rounding_size={h / 2}",
                                fc=PINK, ec=INK, lw=3, zorder=3))
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, text, color="black", fontsize=27, ha="center", va="center", zorder=4,
            **FONT)


def arrow(x, y0, y1):
    ax.annotate("", xy=(x, y1), xytext=(x, y0), zorder=5,
                arrowprops=dict(arrowstyle="-|>,head_width=0.45,head_length=0.9", lw=7, color="black",
                                shrinkA=0, shrinkB=0))


# simulate
group(490, 50, 2050, 263)
command(462, 3, 951, 97, "krsim simulate")
pill(531, 133, 1062, 227, "ts_{seed}.trees")
pill(1126, 133, 1941, 227, "sim_ind_node_map_{seed}.csv")
# mutate
group(490, 445, 1197, 662)
command(462, 405, 951, 497, "krsim mutate")
pill(541, 535, 1142, 628, "ts_mut_{seed}.trees")
# coal-mut
group(33, 853, 738, 1203)
command(3, 807, 520, 900, "krsim coal-mut")
pill(83, 937, 690, 1030, "coal_T_{seed}.csv")
pill(83, 1060, 690, 1155, "mut_T_{seed}.csv")
# vcf-strain
group(873, 921, 1699, 1167)
command(840, 880, 1331, 973, "krsim vcf-strain")
pill(916, 1022, 1655, 1117, "output_{seed}_strain.vcf")
# vcf-host-depth
group(1080, 1365, 2050, 1620)
command(1050, 1322, 1715, 1416, "krsim vcf-host-depth")
pill(1110, 1470, 2010, 1565, "output_{seed}_host_depth.vcf")

arrow(1008, 227, 533)   # ts -> mutate
arrow(595, 628, 853)    # ts_mut -> coal-mut
arrow(1085, 628, 878)   # ts_mut -> vcf-strain
arrow(1380, 1117, 1320)  # strain VCF -> vcf-host-depth
arrow(1848, 227, 1363)  # host-strain map -> vcf-host-depth

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline.png")
fig.savefig(out, dpi=140, facecolor="white")
print(f"wrote {out}")
