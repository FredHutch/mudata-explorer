from typing import Optional
from anndata import AnnData
from biom import load_table
import pandas as pd
from mudata_explorer.helpers import cirro
from mudata_explorer.parsers import curatedMetagenomicData
from mudata_explorer.parsers import microbiome, metaphlan
from mudata_explorer.apps.helpers import load_mdata
import streamlit as st
from tempfile import NamedTemporaryFile


def _read_table(
    label: str,
    exts=["csv", "xlsx", "tsv"],
    index_col=0,
    **kwargs
) -> Optional[pd.DataFrame]:
    file = st.file_uploader(label, type=exts)
    if file is None:
        return

    # Read the table
    if file.name.endswith(".csv"):
        return pd.read_csv(
            file,
            index_col=index_col,
            skiprows=_find_skiprows(file, sep=","),
            **kwargs
        )
    elif file.name.endswith(".tsv"):
        return pd.read_csv(
            file,
            sep='\t',
            index_col=index_col,
            skiprows=_find_skiprows(file, sep="\t"),
            **kwargs
        )
    else:
        return pd.read_excel(file, index_col=index_col, **kwargs)


def _find_skiprows(file, sep=","):
    """
    Find the number of rows to skip in a CSV file.
    This is useful for reading files with headers that are not in the first row.

    1. Find the number of fields in every row
    2. Find the modal number of fields
    3. Find the first row with the modal number of fields
    """
    lines = file.readlines()
    file.seek(0)

    # Find the number of fields in each row
    num_fields = [len(line.decode().split(sep)) for line in lines]

    # Find the modal number of fields
    modal_num_fields = max(set(num_fields), key=num_fields.count)

    # Find the first row with the modal number of fields
    for i, num in enumerate(num_fields):
        if num == modal_num_fields:
            return i
    return 0


def _load_data_csv():
    st.write("""Load Data From Abundance Tables (CSV/TSV/XLSX)

- Required: Table of organism abundances in wide format (samples x organisms)
- Optional: Metadata for samples - the first column should be the sample names
- Optional: Metadata for organisms - the first column should be the organism names
""")

    abund = _read_table("Upload the abundance table")
    if abund is None:
        return
    
    # If the samples are in the columns, transpose the table
    if st.selectbox("Each sample is a different:", ["Row", "Column"]) == "Column":
        abund = abund.T

    # If the values need to be converted to proportions
    if st.checkbox("Convert values to proportions"):
        abund = abund.apply(lambda x: x / x.sum(), axis=1)

    # Optional: Upload sample metadata
    obs = _read_table("Upload sample metadata (optional)")

    # Let the user add a prefix to the sample names
    prefix = st.text_input("Prefix for sample names (optional)")
    if prefix and len(prefix) > 0:
        # Add the prefix to the sample names
        abund = abund.rename(index=lambda x: f"{prefix}{x}")
        if obs is not None:
            obs = obs.rename(index=lambda x: f"{prefix}{x}")

    # Optional: Upload organism metadata
    var = _read_table("Upload organism metadata (optional)")

    if var is not None:
        groupby_var = st.checkbox(
            "Collapse organisms by taxonomic rank (from the organism metadata)"
        )
    else:
        groupby_var = False

    # Align the axes of the abundance table and the metadata
    if obs is not None:
        obs = obs.reindex(abund.index)

    # Parse the data into an AnnData object
    adata = AnnData(abund, obs=obs, var=var)

    mdata = microbiome.parse_adata(adata, groupby_var=groupby_var)
    load_mdata(mdata)


def _load_data_biom():
    st.write('Load Data From BIOM File')

    biom = st.file_uploader("Upload a BIOM file", type=["biom"])
    if biom is None:
        return

    with NamedTemporaryFile() as f:
        f.write(biom.getvalue())
        f.seek(0)
        table = load_table(f.name)
    adata = table.to_anndata()

    mdata = microbiome.parse_adata(adata, groupby_var=True)
    load_mdata(mdata)


def _load_data_cmd():
    st.markdown(
        '''Load Curated Metagenomic Data
        
Analyze abundance data produced by the
[Curated Metagenomic Data](https://waldronlab.io/curatedMetagenomicData/)
project.''')

    file = st.file_uploader("Curated Metagenomic Data")
    if file is None:
        return

    # Parse the curatedMetagenomicData format as AnnData
    adata = curatedMetagenomicData.parse_df(
        pd.read_csv(sep="\t")
        .rename(columns=dict(sample="sample_id"))
    )

    # Run the microbiome analysis
    mdata = microbiome.parse_adata(adata, groupby_var=True)
    load_mdata(mdata)


def _load_data_metaphlan():
    st.write("""Load Data From MetaPhlAn Merged Abundance Table""")

    abund = _read_table("Upload the abundance table", comment="#")
    if abund is None:
        return

    # The user must select which taxonomic level to use
    tax_level = st.selectbox(
        "Taxonomic Level",
        ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"],
        index=6
    )

    # Parse the table using the syntax of MetaPhlAn
    abund, tax = metaphlan.parse_df(abund, tax_level)

    # Show the top organisms
    st.write("Top Organisms:")
    st.dataframe(
        abund
            .mean()
            .sort_values(ascending=False).head(10)
            .rename_axis("Organism")
            .rename("Mean Abundance"),
        use_container_width=True
    )

    # Optional: Upload sample metadata
    obs = _read_table("Upload sample metadata (optional)")

    # Parse the data into an AnnData object
    adata = AnnData(abund, obs=obs, var=tax)

    mdata = microbiome.parse_adata(adata, groupby_var=False)
    load_mdata(mdata)


def run():
    # Let people load data from Cirro
    with st.container(border=1):
        st.write("#### Load From Cirro")
        cirro.load_from_cirro(
            filter_process_ids=[
                "ingest_biom",
                "process-nf-core-ampliseq-2-4-0",
                "curated_metagenomic_data",
                "process-hutch-metaphlan-paired-1-0",
                "hutch-gig-map-align-pangenome-1-0",
                "hutch-sourmash-tax"
            ],
            show_link=False
        )

    with st.container(border=1):
        st.write("#### Upload Local Files")

        source = st.selectbox(
            "File Format",
            options=[
                "Abundance Tables (CSV/TSV/XLSX)",
                "BIOM File",
                "Curated Metagenomic Data",
                "MetaPhlAn Merged Abundance Table"
            ],
            index=None
        )
        if source == "Abundance Tables (CSV/TSV/XLSX)":
            _load_data_csv()
        elif source == "BIOM File":
            _load_data_biom()
        elif source == "Curated Metagenomic Data":
            _load_data_cmd()
        elif source == "MetaPhlAn Merged Abundance Table":
            _load_data_metaphlan()
