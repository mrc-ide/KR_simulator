# Tutorial

[`tutorial.ipynb`](tutorial.ipynb) runs every `krsim` command on a small simulation (3 hosts, 10 generations, seed 123) and shows what each one writes. It is saved with its outputs, so you can read it here on GitHub without running anything.

To run it yourself, install krsim from the repo root (`pip install -e .`, see the main [README](../README.md#installation)), then open the notebook from this folder:

```bash
pip install notebook      # Jupyter isn't installed with krsim
cd tutorial
jupyter notebook tutorial.ipynb
```

Running it writes its files to `tutorial/output/` (not tracked by git), including the two interactive HTML visualizations, which the notebook links to.
