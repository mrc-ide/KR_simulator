# Tutorial

[`tutorial.ipynb`](tutorial.ipynb) runs every `krsim` command on a small simulation (3 hosts, 10 generations, seed 123) and shows what each one writes. It is saved with its outputs, so you can read it here on GitHub without running anything.

To run it yourself, install krsim as in the main [README](../README.md#installation) (the environment includes Jupyter and matplotlib), then open the notebook with the `krsim` environment active:

```bash
jupyter notebook tutorial.ipynb
```

If you installed from a clone for development, use `pip install -e ".[tutorial]"` from the repo root to get the same extras.

Running it writes its files to `output/` next to the notebook (`tutorial/output/` in a clone, not tracked by git), including the two interactive HTML visualizations. Copies of those two pages from the same run are kept in [`examples/`](examples/) so the notebook's links work without running anything. GitHub shows them as source code, so download them and open them in a browser.
