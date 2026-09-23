"""True pairwise IBD fraction between strains, from the coalescence time table."""

from krsim.sim.params import GENOME_LENGTH


def strain_ibd_fractions(coal, time_boundary):
    """
    For each strain pair in `coal` (the coal_T table), the fraction of the
    genome whose coalescence time is <= `time_boundary`. Returns a DataFrame
    with columns sample_1, sample_2, frac_IBD.
    """
    # segments whose coalescence time falls at or below the boundary count as IBD
    ibd = coal[coal["coal_time"] <= time_boundary].copy()
    ibd["segment_length"] = ibd["right"] - ibd["left"]

    ibd_length = ibd.groupby(["sample_1", "sample_2"])["segment_length"].sum()

    pairs = coal[["sample_1", "sample_2"]].drop_duplicates().reset_index(drop=True)
    pairs["frac_IBD"] = pairs.apply(
        lambda row: ibd_length.get((row["sample_1"], row["sample_2"]), 0.0) / GENOME_LENGTH,
        axis=1,
    )
    return pairs
