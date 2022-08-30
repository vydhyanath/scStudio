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

#from SinCel_funcs import *
hv.extension('bokeh', logo=False)
st.set_page_config(layout="wide", page_icon = "/home/raghukiran/scStudio/icon1.png")

print('Running scRNA Pipeline APP')

def plot_bar(datafr, value_list):
    col1,col2 = st.columns([2, 1])
    with col2:
        st.selectbox('Variable to plot',value_list, key = 'col', on_change = on_select)
        if 'update_column' not in st.session_state:
            st.session_state.update_column = st.session_state.col
    with col1:
        bar_plot = hv.Bars(data = datafr[st.session_state.update_column].value_counts()).opts(height = 500, \
            width = 800, invert_axes = True, tools = ['hover'], xticks = 10)
        return st.bokeh_chart(hv.render(bar_plot, backend = 'bokeh'))
        

        


TITLE = st.container()


with TITLE:
#    st.title('scStudio')
    st.markdown('<p style="font-family:sans-serif; color:#04410D; font-weight:bold; font-size: 45px;">scStudio</p>', \
        unsafe_allow_html = True)
    st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-size: 18px;">Catalogue of tools for analysing single-cell sequencing data to identifying biomarkers with added exploratory data analysis tools and visualization tools.This tool also supports between the cell types and within the cell conditions</p>', \
                unsafe_allow_html = True)

def on_upload():
    if st.session_state.upload:
        st.session_state.up = st.session_state.upload

def file_uploader():
#    col1, col2, col3 = st.columns([4, 1.5, 0.7])
    with open('app/STYLE.css') as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
#    with col1:
    uploaded_file = st.file_uploader("Choose .h5ad file", key = 'upload', on_change = on_upload, type = ['h5ad'])
    if 'up' not in st.session_state:
        st.session_state.up = st.session_state.upload
    if st.session_state.up is not None:
#        st.write(f"uploaded file : {uploaded_file.name}")
        dataf = scanpy.read_h5ad(st.session_state.up)
        ob, va = dataf.shape
        
#        with col2:
#            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cells<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{ob}</strong></p>''',  unsafe_allow_html=True)
#        with col2:
#            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Genes<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{va}</strong>''',  unsafe_allow_html=True)
#        with col3:
#            st.image('cell-icon-12.jpg', output_format = 'JPEG', use_column_width = 'always')

        return dataf, True, ob, va
    else:
        st.markdown('<p style="font-family:sans-serif; color:#03673D; font-size: 15px;">Upload your file</p>', \
            unsafe_allow_html = True)
        ob, va, all_celltypes, all_samples = 0,0,0,0
#        with col2:
#            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cells<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{ob}</strong></p>''',  unsafe_allow_html=True)
#        with col2:
#            st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Genes<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{va}</strong>''',  unsafe_allow_html=True)

#        with col3:
#            st.image('cell-icon-12.jpg', output_format = 'JPEG')       
        return False,False, ob, va




with open('app/STYLE.css') as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)


adata, upfile, OBS, VAR = file_uploader()

def on_select():
    if st.session_state.col:
        st.session_state.update_column = st.session_state.col

def meta_info(adata):
    st.header('General Meta data')
    adata_stats = lambda x : adata.obs[x].nunique() \
        if adata.obs[x].dtype == 'category' else adata.obs[x].sum() 
    TABLE = st.container()


    tab1, tab2 = st.tabs(['Plot','Table'])

    meta_df = pd.DataFrame({col : adata_stats(col) for col in adata.obs.columns}.items(), \
            columns = ['Fields', 'Unique Values'])
    table_plot = hv.Table(meta_df)


    obs_dtypes = adata.obs.dtypes.to_frame(name='type')

    def get_count():
        cols = [col for col in obs_dtypes.index if str(obs_dtypes.loc[col].type) == 'category' \
            if adata.obs[col].nunique() != adata.shape[0]]
        obs_dic = dict(zip(cols, \
            [Counter(adata.obs[col]) for col in cols]))
        return obs_dic

    print(get_count())        
    obs_dict = get_count()


    acol1,acol2 = st.columns([2, 1])
    with tab1:
        plot_bar(adata.obs, obs_dict.keys())
#        with acol1:
#            bar_plot = hv.Bars(data = adata.obs[st.session_state.update_column].value_counts()).opts(height = 500, \
#                width = 800, invert_axes = True, tools = ['hover'], xticks = 10)
#            st.bokeh_chart(hv.render(bar_plot, backend = 'bokeh'))
#        with acol2:
#            st.selectbox('Variable to plot',obs_dict.keys(), key = 'col', on_change = on_select)
#            if 'update_column' not in st.session_state:
#                st.session_state.update_column = st.session_state.col
    with tab2:
        st.bokeh_chart(hv.render(table_plot, backend = 'bokeh'))

col1, col2 = st.columns([1, 1])

with col1:
    st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cells<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{OBS}</strong></p>''',  unsafe_allow_html=True)
with col2:
    st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Genes<br></p><p style="font-family:Arial; font-size: 20px;"><strong>{VAR}</strong>''',  unsafe_allow_html=True)



st.session_state["uploaded"] = upfile
st.session_state["adata_obj"] = adata
st.session_state["OBS"] = OBS
st.session_state["VAR"] = VAR
#st.session_state["file_name"] = file_name

#st.sidebar.markdown(f'''><p style="font-family:Arial; font-size: 15px;">file name<br></p><p style="font-family:Arial; font-size: 16px;"><strong>{file_name.rstrip('.h5ad')}</strong></p>''',  \
#    unsafe_allow_html=True)


if upfile:
    meta_info(adata)

st.sidebar.markdown('---')
st.sidebar.write('contact raghukiran812@gmail.com for suggestions')

