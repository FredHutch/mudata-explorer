from typing import Optional, Union, List
from cirro import CirroApi, DataPortal, DataPortalProject
from cirro import DataPortalDataset
from cirro.auth.device_code import DeviceCodeAuth
from cirro.config import AppConfig, list_tenants
from io import StringIO
from mudata import MuData
from mudata_explorer.app.mdata import get_mdata, set_mdata
from mudata_explorer.app.hash import get_dat_hash, set_mdata_hash
from mudata_explorer.helpers.io import hydrate_uns
from mudata_explorer.helpers.cirro_readers import util, mudata, ampliseq, biom, metaphlan, gig_map
from mudata_explorer.helpers.cirro_readers import sourmash
from mudata_explorer.helpers.cirro_readers import differential_abundance
from mudata_explorer.helpers.cirro_readers import curatedMetagenomicData
import streamlit as st
from streamlit.delta_generator import DeltaGenerator
try:
    from streamlit.runtime.scriptrunner_utils import script_run_context
except ImportError:
    from streamlit.runtime.scriptrunner import script_run_context
from tempfile import TemporaryDirectory
from threading import Thread
from time import sleep


def setup_cirro_client():

    if st.query_params.get("domain") is None:
        tenant_dict = {
            tenant['displayName']: tenant['domain']
            for tenant in list_tenants()
        }

        # Let the user select a tenant
        tenant = st.selectbox(
            "Organization",
            ["< select for login >"] + list(tenant_dict.keys())
        )

        domain = tenant_dict.get(tenant)

    else:
        domain = st.query_params["domain"]

    if domain:
        _cirro_login(domain, st.empty())


def _cirro_login(domain: str, container: DeltaGenerator):

    # Connect to Cirro - capturing the login URL
    auth_io = StringIO()
    cirro_login_thread = Thread(
        target=_cirro_login_sub,
        args=(auth_io, domain)
    )
    script_run_context.add_script_run_ctx(cirro_login_thread)

    cirro_login_thread.start()

    login_string = auth_io.getvalue()

    while len(login_string) == 0 and cirro_login_thread.is_alive():
        sleep(1)
        login_string = auth_io.getvalue()

    container.write(login_string)
    cirro_login_thread.join()
    container.empty()
    st.rerun()


def _cirro_login_sub(auth_io: StringIO, base_url: str):

    app_config = AppConfig(base_url=base_url)

    st.session_state['Cirro-auth_info'] = DeviceCodeAuth(
        region=app_config.region,
        client_id=app_config.client_id,
        auth_endpoint=app_config.auth_endpoint,
        enable_cache=False,
        auth_io=auth_io
    )

    st.session_state['Cirro-client'] = CirroApi(
        auth_info=st.session_state['Cirro-auth_info'],
        base_url=base_url
    )
    st.session_state['Cirro'] = DataPortal(
        client=st.session_state['Cirro-client']
    )


def load_from_cirro(
    filter_process_ids: Optional[List[str]] = None,
    show_link=True,
    id="main"
):

    # If the Cirro client has not been set up,
    # prompt the user to set it up
    if not st.session_state.get("Cirro"):
        setup_cirro_client()
        return

    # Ask the user to select a project
    project = _select_project("load")
    if project is None:
        return

    # Ask the user to select a dataset
    dataset = _select_dataset(project, filter_process_ids=filter_process_ids)
    if dataset is None:
        return
        
    # If the selection was hardcoded, just load the dataset
    # But if not, ask the user if they want to load it
    if show_link and dataset.id != st.query_params.get("dataset_id"):

        cols = st.columns([1, 5])

        # Show a link to the dataset
        url = _load_dataset_url(project.id, dataset.id)
        cols[1].code(url, language=None)

        # Ask the user if they want to load the dataset
        if not cols[0].button("Load", use_container_width=True):
            return

    # Read the MuData object from the contents of the dataset
    mdata = _read_dataset(dataset)

    # If no data was read, stop
    if mdata is None:
        return

    # Get the hash of the data
    _, hash, _ = get_dat_hash(mdata)

    hydrate_uns(mdata)
    set_mdata(mdata, full=True, id=id)
    set_mdata_hash(hash, id=id)
    st.switch_page("pages/views.py")


def _load_dataset_url(
    project_id: str,
    dataset_id: str
) -> str:
    """Get the URL to load the dataset from Cirro."""
    base_url = st.session_state.get("Cirro")._client.configuration.base_url
    app_url = "https://mudata-explorer.streamlit.app"
    query_args = "&".join([
        f"{kw}={val}"
        for kw, val in dict(
            domain=base_url,
            project_id=project_id,
            dataset_id=dataset_id
        ).items()
    ])
    return f"{app_url}/load?{query_args}"


def save_to_cirro(id="main"):

    st.write("#### Save to Cirro")

    # If the Cirro client has not been set up,
    # prompt the user to set it up
    if not st.session_state.get("Cirro"):
        setup_cirro_client()
        return

    # Ask the user to select a project
    project = _select_project("save")
    if project is None:
        return

    # Ask the user to provide a name and description for the dataset
    name = st.text_input("Name")
    description = st.text_area("Description")

    if not name:
        st.error("Name is required")
        return

    if not description:
        description = ""

    if not st.button("Save"):
        return

    # Get the active dataset
    mdata = get_mdata(id=id, full=True)

    # If there is no active dataset
    if mdata is None:
        st.error("No dataset to save")
        return

    # Get the binary blob, hash, and file size
    blob, hash, size = get_dat_hash(mdata)

    assert blob is not None

    # Set the file name
    basename = name.replace(' ', '-').lower().replace('/', '-')
    while "--" in basename:
        basename = basename.replace("--", "-")
    fn = f"{basename}-{hash}.h5mu"

    # Write out the dataset to a temporary file
    # and upload it to Cirro
    with TemporaryDirectory() as tmp:

        # Write the MuData object to the file
        with open(f"{tmp}/{fn}", "wb")as handle:
            handle.write(blob)

        # Upload the file to Cirro
        try:
            dataset = project.upload_dataset(
                name=name,
                description=description,
                process="mudata-h5mu",
                upload_folder=tmp
            )
        except Exception as e:
            st.exception(e)
            return

    st.write(
        f"Wrote {name} to {project.name} ({size})"
    )
    # Show a link to the dataset
    url = _load_dataset_url(project.id, dataset.id)
    st.markdown(f"[Permalink]({url})")


def _select_project(key: str) -> DataPortalProject:
    """Get the list of projects available from Cirro."""

    portal: DataPortal = st.session_state.get("Cirro")
    if not portal:
        util.clear_cirro_client()

    # Try to get the list of projects
    try:
        projects = portal.list_projects()
    # If there is an error
    except Exception as e:
        # Report it to the user
        st.exception(e)
        util.clear_cirro_client()

    # If a project is supplied in the query string, use it
    if "project_id" in st.query_params:
        for p in projects:
            if p.id == st.query_params["project_id"]:
                return p
        raise ValueError(f"Project {st.query_params['project_id']} not found")
    
    # Sort the projects by name
    projects.sort(key=lambda p: p.name)

    # Give the user the option to select a project
    project = st.selectbox(
        "Project",
        [p.name for p in projects],
        key=f"{key}_project",
        index=None,
        placeholder="< select a project >"
    )

    if st.session_state.get(f"{key}_project") is None:
        return

    # Return the project object
    for p in projects:
        if p.name == st.session_state.get(f"{key}_project"):
            return p
    raise ValueError(f"Project '{project}' not found")


def _select_dataset(
    project: DataPortalProject,
    filter_process_ids: Optional[List[str]] = None
) -> DataPortalDataset:

    # Try to get the list of datasets in the project
    try:
        datasets = project.list_datasets()
    except Exception as e:
        st.exception(e)
        util.clear_cirro_client()

    # Filter down to the datasets which have parsing configured
    # If filter_process_ids is provided, only those will be allowed
    datasets = [
        d for d in datasets
        if _read_dataset(
            d,
            check_only=True,
            filter_process_ids=filter_process_ids
        )
    ]

    # Sort by date
    datasets.sort(key=lambda d: d.created_at, reverse=True)

    # If a dataset is supplied in the query string, use it
    if "dataset_id" in st.query_params:
        for p in datasets:
            if p.id == st.query_params["dataset_id"]:
                return p
        raise ValueError(f"Dataset {st.query_params['dataset_id']} not found")

    # Give the user the option to select a dataset
    if datasets:
        st.write(":arrow_down_small: Select a dataset")
        selection = st.dataframe(
            [
                dict(
                    Name=ds.name,
                    Description=ds.description,
                    Created=ds.created_at,
                    User=ds.created_by,
                    id=ds.id
                )
                for ds in datasets
            ],
            column_config=dict(id=None),
            hide_index=True,
            key="select-dataset",
            selection_mode="single-row",
            on_select="rerun"
        )
    else:
        st.write("No datasets available")
        return

    if len(selection.get("selection", {}).get("rows", [])) == 0:
        return

    # Return the dataset object
    return datasets[selection["selection"]["rows"][0]]


def _read_dataset(
    dataset: DataPortalDataset,
    check_only: bool = False,
    filter_process_ids: Optional[List[str]] = None
) -> Union[bool, Optional[MuData]]:
    """
    Read a dataset from Cirro.
    Report any errors that arise.
    If check_only is True, only check if the dataset's type is valid
    """

    # If a list of process IDs is supplied, only allow those
    if filter_process_ids is not None:
        if dataset.process_id not in filter_process_ids:
            if check_only:
                return False
            else:
                return None

    # MuData datasets
    if dataset.process_id == "mudata-h5mu":
        if check_only:
            return True
        else:
            return mudata.read(dataset)

    # nf-core/ampliseq datasets
    elif dataset.process_id == "process-nf-core-ampliseq-2-4-0":
        if check_only:
            return True
        else:
            return ampliseq.read(dataset)

    # biom datasets
    elif dataset.process_id == "ingest_biom":
        if check_only:
            return True
        else:
            return biom.read(dataset)

    # curatedMetagenomicData datasets
    elif dataset.process_id == "curated_metagenomic_data":
        if check_only:
            return True
        else:
            return curatedMetagenomicData.read(dataset)

    # CirroBio/nf-differential-abundance datasets
    elif dataset.process_id in [
        f"differential-abundance-{suffix}"
        for suffix in [
            "metaphlan",
            "ampliseq",
            "CMD",
            "gig-map-align-reads",
            "table",
            "sourmash"
        ]
    ]:
        if check_only:
            return True
        else:
            return differential_abundance.read(dataset)

    # Metaphlan taxonomic classification datasets
    elif dataset.process_id == "process-hutch-metaphlan-paired-1-0":
        if check_only:
            return True
        else:
            return metaphlan.read(dataset)

    # Sourmash taxonomic classification datasets
    elif dataset.process_id == "hutch-sourmash-tax":
        if check_only:
            return True
        else:
            return sourmash.read(dataset)

    # gig-map alignment datasets
    elif dataset.process_id == "hutch-gig-map-align-pangenome-1-0":
        if check_only:
            return True
        else:
            return gig_map.read(dataset)

    elif check_only:
        return False
