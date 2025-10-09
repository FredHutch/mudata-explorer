from typing import Tuple
import pandas as pd


def parse_df(abund: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # The samples are in the columns, so transpose the table
    abund = abund.T

    # Build a taxonomy table by splitting up the taxonomic strings

    tax = pd.DataFrame([
        dict(zip(
            ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
            cname.split(";")
        ))
        for cname in abund.columns
    ], index=abund.columns)

    # Rename the taxa, both in the table and the taxonomy
    name_map = {
        col: col.split(";")[-1]
        for col in abund.columns
    }

    # Clean up those names a bit more
    name_map = {
        kw: val.split("__")[-1] if len(val) > 4 and val[1:].startswith("__") else val
        for kw, val in name_map.items()
    }

    abund.rename(columns=name_map, inplace=True)
    tax.rename(index=name_map, inplace=True)

    # Make sure that the number of dimensions adds up
    assert abund.shape[1] == tax.shape[0], "Number of taxa must match"

    return abund, tax
