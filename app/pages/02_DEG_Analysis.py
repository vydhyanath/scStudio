#from ..Summary import on_select
import streamlit as st
import pandas as pd
import scanpy
import numpy as np
import holoviews as hv
import os
from altair.vegalite.v4.schema.channels import Tooltip
from collections import Counter
import altair as alt
from holoviews import dim


st.set_page_config(layout="wide", page_icon = "/home/raghukiran/scStudio/icon1.png")
st.header('Differential Gene Expression Data analysis')

def get_count(adata,include_only = 'category'):
    obs_dtypes = adata.obs.dtypes.to_frame(name='type')
    if include_only == 'category':
        cols = [col for col in obs_dtypes.index if str(obs_dtypes.loc[col].type) == 'category' \
            if adata.obs[col].nunique() != adata.shape[0]]
    else:
        cols = [col for col in obs_dtypes.index if str(obs_dtypes.loc[col].type) == 'category' or \
            str(obs_dtypes.loc[col].type).rstrip('3264') == 'int' or str(obs_dtypes.loc[col].type).rstrip('3264') == 'float' \
            if adata.obs[col].nunique() != adata.shape[0]]

    obs_dic = dict(zip(cols, \
        [Counter(adata.obs[col]) for col in cols]))
    return obs_dic

def on_group():
    if st.session_state.cell_group:
        st.session_state.update_cell_group = st.session_state.cell_group

def onchange_A():
    if st.session_state.setA:
        st.session_state.update_setA = st.session_state.setA

def onchange_B():
    if st.session_state.setB:
        st.session_state.update_setB = st.session_state.setB

def update_pval():
    if st.session_state.adjpval:
        st.session_state.update_adjpval = st.session_state.adjpval

def update_method():
    if st.session_state.method:
        st.session_state.update_method = st.session_state.method

def update_click():
    if st.session_state.click:
        st.session_state.update_click = st.session_state.click
#@st.experimental_memo(suppress_st_warning=True)
def perform_cell_deg(adata):
    param_container = st.container()
    col1, col2, col3 = st.columns([2,1,1])
    adata_raw = adata.raw.to_adata()
    adata_raw.uns['log1p']['base'] = None
    methods = ['t-test', 'wilcoxon', 't-test_overestim_var']
    with col1:
        st.selectbox('Group by', get_count(adata_raw).keys(), key = 'cell_group', on_change = on_group)
        if 'update_cell_group' not in st.session_state:
            st.session_state.update_cell_group = st.session_state.cell_group
    setA_values = adata_raw.obs[st.session_state.update_cell_group].unique()
    setB_values = adata_raw.obs[st.session_state.update_cell_group].unique().tolist() + ['rest']
    with col2:
        st.selectbox('Set A', setA_values, key = 'setA', on_change = onchange_A)
        if 'update_setA' not in st.session_state:
            st.session_state.update_setA = st.session_state.setA
    with col3:
        st.selectbox('Set B', setB_values, key = 'setB', on_change = onchange_B)
        if 'update_setB' not in st.session_state:
            st.session_state.update_setB = st.session_state.setB
    
    acol1, acol2, _,_ = st.columns([1.5,0.25, 0.25, 0.25])
    with acol1:
        st.radio('Method', methods, key = 'method', on_change = update_method, horizontal = True)
        if 'update_method' not in st.session_state:
            st.session_state.update_method = st.session_state.method

    with acol2:
        st.number_input('AdjPval cutoff', value = 0.05,min_value = 0.0, max_value = 1.0, key = 'adjpval', on_change = update_pval)
        if 'update_adjpval' not in st.session_state:
            st.session_state.update_adjpval = st.session_state.adjpval
        print(st.session_state.update_adjpval)
    _, _, _,bcol4 = st.columns([1.5,0.25, 0.25, 0.25])
    with bcol4:
        run = st.button('Run')
#        st.button('Run', key = 'click', on_click = update_click)
#        if 'update_click' not in st.session_state:
#            st.session_state.update_click = st.session_state.click

#    if st.session_state.update_click:
    if run:
        st.spinner(f"Running !")
        scanpy.tl.rank_genes_groups(adata_raw, st.session_state.update_cell_group, method = st.session_state.update_method, reference = st.session_state.update_setB, \
                    group = [st.session_state.update_setA])
        st.success(f"successfully ran {st.session_state.update_method}")
        return scanpy.get.rank_genes_groups_df(adata_raw,st.session_state.update_setA, key = 'rank_genes_groups', pval_cutoff=st.session_state.update_adjpval)


def between_cell_deg(adata):
    st.markdown('<p style="font-family:Arial; color:#CA0327; font-weight:bold; font-size:20px">Between Cell Groups</p>', unsafe_allow_html=True)
    st.bokeh_chart(hv.render(hv.Table(perform_cell_deg(adata)).opts(width = 800, height = 500), backend = 'bokeh'))
#    st.info('Work in Progress')

def pseudo_bulk_deg(adata):
    st.markdown('<p style="font-family:Arial; color:#CA0327; font-weight:bold; font-size:20px">Between Sample Groups</p>', unsafe_allow_html=True)
    st.info('Work in Progress')

def deg_info(adata):
    st.sidebar.markdown('<p style="font-family:sans-serif; font-weight:bold ;color:#ff9633 ;font-size: 18px;">Differential Gene Expression</p>',\
        unsafe_allow_html = True)
    deg_sec = st.sidebar.selectbox(
        "Identify differential genes between the cell types or with in the cell conditions",
        ("Cell Types", "Pseudo-bulk analysis")
    )
    if deg_sec == "Cell Types":
        between_cell_deg(adata)
    
    elif deg_sec == "Pseudo-bulk analysis":
        pseudo_bulk_deg(adata)

#adata, upfile = file_uploader()

#if upfile:
#    deg_info(adata)

col1, col2 = st.columns([1, 1])

with open('app/STYLE.css') as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

if 'OBS' not in st.session_state:
    st.session_state['OBS'] = 0
if 'VAR' not in st.session_state:
    st.session_state['VAR'] = 0
if 'uploaded' not in st.session_state:
    st.session_state['uploaded'] = False

with col1:
    st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cells<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{st.session_state['OBS']}</strong></p>''',  unsafe_allow_html=True)
with col2:
    st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Genes<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{st.session_state['VAR']}</strong>''',  unsafe_allow_html=True)

if st.session_state["uploaded"]:    
    deg_info(st.session_state["adata_obj"])
else:
    st.write("upload the file")