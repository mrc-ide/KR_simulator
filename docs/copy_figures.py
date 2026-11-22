"""Copy the figures the website uses from figures/ (their only home in git) into docs/images/.

Quarto runs this before every render (pre-render in _quarto.yml), locally and on GitHub. docs/images/ is not in
git. To show a new figure on a page, add it to WEBSITE and use images/<name> in the page.
"""

import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURES = os.path.join(HERE, "..", "figures")
IMAGES = os.path.join(HERE, "images")

WEBSITE = [
    "simulator/KR_model.png",
    "simulator/KR_tskit.png",
    "simulator/pipeline.png",
    "simulator/scaling_pipeline.png",
    "pf8_data/depth_along_genome_PF1066-Cx.png",
    "pf8_data/depth_spread_sWGA_samples.png",
]

os.makedirs(IMAGES, exist_ok=True)
for f in WEBSITE:
    shutil.copy2(os.path.join(FIGURES, f), os.path.join(IMAGES, os.path.basename(f)))
