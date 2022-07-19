import streamlit as st
import pandas as pd
import scanpy
from SinCel_funcs import *

st.set_page_config(layout="wide")
#Navigation bar
#st.markdown('<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-TX8t27EcRE3e/ihU7zmQxVncDAy5uIKz4rEkgIXeMed4M0jlfIDPvg6uqKI2xXr2" crossorigin="anonymous">', unsafe_allow_html=True)
#
#st.markdown("""
#<nav class="navbar fixed-top navbar-expand-lg navbar-dark" style="background-color: #21d88e;">
#  <a class="navbar-brand" href="">ScRNA-pipl</a>
#  <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarSupportedContent" aria-controls="navbarSupportedContent" aria-expanded="false" aria-label="Toggle navigation">
#    <span class="navbar-toggler-icon"></span>
#  </button>
#  <div class="collapse navbar-collapse" id="navbarNav">
#    <ul class="navbar-nav">
#      <li class="nav-item active">
#        <a class="nav-link disabled" href="#">Home <span class="sr-only">(current)</span></a>
#      </li>
#    </ul>
#  </div>
#</nav>
#""", unsafe_allow_html=True)

print('Running scRNA Pipeline APP')

TITLE = st.container()
#MARKDOWN = st.container()
#FILE_UPLOADER = st.container()

##st.markdown(
##    '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@4.5.3/dist/css/bootstrap.min.css" integrity="sha384-TX8t27EcRE3e/ihU7zmQxVncDAy5uIKz4rEkgIXeMed4M0jlfIDPvg6uqKI2xXr2" crossorigin="anonymous">',
##    unsafe_allow_html=True,
##)
##
##query_params = st.experimental_get_query_params()
##tabs = ["Home", "About", "Contact"]
##if "tab" in query_params:
##    active_tab = query_params["tab"][0]
##else:
##    active_tab = "Home"
##
##if active_tab not in tabs:
##    st.experimental_set_query_params(tab="Home")
##    active_tab = "Home"
##
##li_items = "".join(
##    f"""
##    <li class="nav-item">
##        <a class="nav-link{' active' if t==active_tab else ''}" href="/?tab={t}">{t}</a>
##    </li>
##    """
##    for t in tabs
##)
##tabs_html = f"""
##    <ul class="nav nav-tabs">
##    {li_items}
##    </ul>
##"""
##st.markdown(tabs_html, unsafe_allow_html=True)
##st.markdown("<br>", unsafe_allow_html=True)

with TITLE:
#    st.title('scStudio')
    st.markdown('<p style="font-family:sans-serif; color:#04410D; font-weight:bold; font-size: 30px;">scStudio</p>', \
        unsafe_allow_html = True)
    st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-size: 18px;">Catalogue of tools for analysing single-cell sequencing data to identifying biomarkers with added exploratory data analysis tools and visualization tools.This tool also supports between the cell types and within the cell conditions</p>', \
                unsafe_allow_html = True)

#@st.cache(suppress_st_warning=True, persist=True)
def file_uploader():
    col1, col2, col3 = st.columns([3, 1, 1])
    with open('app/STYLE.css') as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
    with col1:
        uploaded_file = st.file_uploader("Choose .h5ad file")
    if uploaded_file is not None:
        dataf = scanpy.read_h5ad(uploaded_file)
#        dataf = dataf.raw.to_adata()
        ob, va = dataf.shape
        all_celltypes, all_samples = dataf.obs['cell_type'].nunique(), dataf.obs['patient_id'].nunique()
        with col2:
            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cells<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{ob}</strong></p>''',  unsafe_allow_html=True)
        with col2:
            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Genes<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{va}</strong>''',  unsafe_allow_html=True)
#        with col3:
#            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cell Types</p><p style="font-family:Arial; font-size: 20px;"><strong>{all_celltypes}</strong></p>''',  unsafe_allow_html=True)
#        with col3:
#            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Total samples</p><p style="font-family:Arial; font-size: 20px;"><strong>{all_samples}</strong></p>''',  unsafe_allow_html=True)
        with col3:
            st.image('cell-icon-12.jpg', output_format = 'JPEG', use_column_width = 'always')
        return dataf, True
    else:
        st.markdown('<p style="font-family:sans-serif; color:#03673D; font-size: 15px;">Upload your file</p>', \
            unsafe_allow_html = True)
        ob, va, all_celltypes, all_samples = 0,0,0,0
        with col2:
            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cells<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{ob}</strong></p>''',  unsafe_allow_html=True)
        with col2:
            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Genes<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{va}</strong>''',  unsafe_allow_html=True)
#        with col3:
#            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cell Types</p><p style="font-family:Arial; font-size: 20px;"><strong>{all_celltypes}</strong></p>''',  unsafe_allow_html=True)
#        with col3:
#            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Total samples</p><p style="font-family:Arial; font-size: 20px;"><strong>{all_samples}</strong></p>''',  unsafe_allow_html=True)
        with col3:
            st.image('cell-icon-12.jpg', output_format = 'JPEG', use_column_width = 'always')       
        return False,False

    

menu = st.sidebar.radio(
    "",
    ("Summary", "Vizualization", "Differential Gene Expression"),
)

fileup_container = st.container()

with open('app/STYLE.css') as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)


adata, upfile = file_uploader()

if upfile:
    if menu == 'Summary':
        meta_info(adata)
#    elif menu == 'Quality Control':
#        qc_info(adata)
    elif menu == 'Vizualization':
        viz_info(adata)
    elif menu == 'Differential Gene Expression':
        deg_info(adata)



st.sidebar.markdown('---')
st.sidebar.write('contact raghukiran812@gmail.com for suggestions')

#9603986566
#31dfa4

