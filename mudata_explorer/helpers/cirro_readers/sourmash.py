from anndata import AnnData
from cirro import DataPortalDataset
from cirro.sdk.file import DataPortalFile
from mudata import MuData
import pandas as pd
import streamlit as st
from typing import Dict, Tuple, List, Optional
from mudata_explorer.helpers.cirro_readers import util
from mudata_explorer.parsers import sourmash
from mudata_explorer.parsers.microbiome import parse_adata


def read(
    dataset: DataPortalDataset
) -> Optional[MuData]:
    """Read datasets produced by sourmash."""

    adata = _read_sourmash_as_anndata(dataset)
    if adata is None:
        return

    return parse_adata(adata, groupby_var=False)


def _read_sourmash_as_anndata(dataset: DataPortalDataset) -> Optional[AnnData]:
    # Get the list of files from the dataset
    files = util.list_files(dataset)

    # Get the options for lineage summaries which are provided
    lineage_summary_prefix = "data/tax_metagenome/all.tax.metagenome."
    lineage_summary_suffix = ".lineage_summary.tsv"
    options = {
        f.name[len(lineage_summary_prefix):-len(lineage_summary_suffix)].title(): f
        for f in files
        if f.name.startswith(lineage_summary_prefix) and f.name.endswith(lineage_summary_suffix)
    }

    # Read in the relative abundance table,
    # which also contains the taxonomic assignment for each feature
    with st.container(border=1):
        rel_abund = _read_rel_abund(options)
    if rel_abund is None:
        return
    abund, tax = rel_abund

    # Read the samplesheet saved with the dataset
    sample_meta = util.sample_metadata(dataset)

    # The sample metadata must have the same rows as the abundance
    sample_meta = sample_meta.reindex(abund.index)

    # Make an AnnData object
    adata = AnnData(
        X=abund,
        var=tax,
        obs=sample_meta
    )
    return adata


def _read_rel_abund(
    options: Dict[str, DataPortalFile]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read the relative abundance table from a dataset."""

    # The user must select which taxonomic level to use
    tax_level = st.selectbox(
        "Taxonomic Level",
        list(options.keys()),
        index=0
    )

    # Read the selected file
    abund = options[tax_level].read_csv(sep="\t", index_col=0)

    # Parse the table using the syntax of Sourmash
    return sourmash.parse_df(abund)
