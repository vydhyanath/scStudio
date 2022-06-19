import os
from altair.vegalite.v4.schema.channels import Tooltip
from matplotlib.pyplot import title
import pandas as pd
import numpy as np
import scanpy
import streamlit as st
from collections import Counter
import altair as alt
import holoviews as hv
from holoviews import dim

hv.extension('bokeh', logo=False)

#def plot_bar(data, var = None):
#    if type(data) == Counter:
#        st_bar = st.bar_chart(pd.DataFrame(data.values(), index = data.keys(), columns = [var]))
#    elif type(data) == dict:
#        st_bar = st.bar_chart(pd.DataFrame(data[var].values(), index = data[var].keys(), columns = [var]))
#    else:
#        st_bar = st.bar_chat(data)
#    return st_bar

#def plot_bar(data, opt = None):
#    datafrm = pd.DataFrame(data[opt].items(), columns=[opt, 'frequency'])
#    st_bar = alt.Chart(datafrm).mark_bar().encode(
#        x = opt, y = 'frequency', tooltip = [alt.Tooltip(opt, title = 'observation'), \
#        alt.Tooltip('frequency', title = 'frequency')]
#    ).properties(height = 600)
#    return st_bar


#def plot_density(data):




def on_select():
    if st.session_state.col:
        st.session_state.update_column = st.session_state.col

#def on_select():
#    if st.session_state.nrows:
#        st.session_state.show = st.session_state.nrows

def between_cell_deg(adata):
    st.markdown('<p style="font-family:Arial; font-size:25px">Differential Gene Expression between Cell Groups</p>', unsafe_allow_html=True)

def pseudo_bulk_deg(adata):
    st.markdown('<p style="font-family:Arial; font-size:25px">Differential Gene Expression within Cells Groups</p>', unsafe_allow_html=True)

#@st.cache(suppress_st_warning=True, persist=True)
def meta_info(adata):
    st.header('General Meta data')
    adata_stats = lambda x : adata.obs[x].nunique() \
        if adata.obs[x].dtype == 'category' else adata.obs[x].sum().round() 
    TABLE = st.container()


    col1, col2, col3= st.columns([2.5, 0.5, 0.5])

    with TABLE:
        meta_df = pd.DataFrame({col : adata_stats(col) for col in adata.obs.columns}.items(), \
            columns = ['Fields', 'Unique Values'])
        table_plot = hv.Table(meta_df)
        st.bokeh_chart(hv.render(table_plot, backend = 'bokeh'))
#***********************important******************************        
#        sel, _ = st.columns([0.45,3.55])
#        with sel:
#            sel_val = st.selectbox('show rows', (5, 10, 50, 100), \
#                on_change = on_select, key = 'nrows')
#            if 'show' not in st.session_state:
#                st.session_state.show = st.session_state.nrows
#        with open('apps/STYLE.css') as f:
#            st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
#
#
#        st.write(adata.obs.head(st.session_state.show))
        obs_dtypes = adata.obs.dtypes.to_frame(name='type')
        obs_dict = dict(zip(obs_dtypes.index, \
            [Counter(adata.obs[col]) for col in obs_dtypes.index \
                if str(obs_dtypes.loc[col].type) == 'category']))

        

    mp_bar = st.container() 
    with mp_bar:
        st.markdown('<p style="font-family:sans-serif; font-size: 18px;">general statistics of single cells observation</p>', \
                unsafe_allow_html = True)

    acol1, _,acol3 = st.columns([1.5, 0.1, 0.35])
    with acol3:
        opt = st.selectbox('Variable to plot',obs_dict.keys(), key = 'col', on_change = on_select)
        if 'update_column' not in st.session_state:
            st.session_state.update_column = st.session_state.col
    with acol1:
#        st.altair_chart(plot_bar(obs_dict, st.session_state.update_column), use_container_width= True)
        bar_plot = hv.Bars(data = adata.obs[st.session_state.update_column].value_counts()).opts(height = 500, \
            width = 700, invert_axes = True, tools = ['hover'], xticks = 10)
        st.bokeh_chart(hv.render(bar_plot, backend = 'bokeh'))












#@st.cache(suppress_st_warning=True, persist=True)
def  qc_info(adata):
    st.header('Quality control of single cell Sequencing data')
    st.sidebar.markdown('<p style="font-family:sans-serif; font-weight:bold ;color:#FDE0D1 ;font-size: 18px;">Quality Control</p>',\
        unsafe_allow_html = True)
    st.sidebar.markdown('<p style="font-family:sans-serif ;font-size: 12px;">Quality Statistics for Single Cell Data</p>',\
        unsafe_allow_html = True)
    
        
    hist_counts = hv.Histogram((np.histogram(adata.obs['log1p_total_counts'], bins=100))).opts(width = 400, \
        height = 300, xlabel = 'counts depth').relabel('Histograms of count depth per cell')
    hist_genes = hv.Histogram((np.histogram(adata.obs['log1p_n_genes_by_counts'], bins=100))).opts(width = 400, \
        height = 300, xlabel = 'genes count').relabel('Histogram of the number of genes detected per cell')
    elb_plt = hv.Points(list(enumerate(adata.obs['n_counts'].sort_values(ascending = False).values))).\
        opts(logy = True, width = 400, height = 300, xlabel = 'Barcode Rank', ylabel = 'log counts depth').\
            relabel('Count depth distribution from high to low count depths')
    genes_counts = hv.Points((adata.obs['n_counts'].values, \
        adata.obs['n_genes'].values)).opts(xlabel = 'count depth', \
            ylabel = 'genes per cell', width = 400, height = 300).relabel('Number of genes versus the count depth')

    qc_layout = hv.Layout(hist_counts + hist_genes + elb_plt + \
        genes_counts).opts(shared_axes = False).cols(2)

    st.bokeh_chart(hv.render(qc_layout, backend = 'bokeh'))

def plot_violine(anndata, top_genes = 10):
    genes = anndata.to_df().mean(0).sort_values(ascending = False).head(top_genes).index
    melt_df  = anndata.to_df().loc[:,genes].melt()
    violine = hv.Violin(melt_df, ['index'], \
            'value', label = f'Top {top_genes} Highly expressed genes').opts(height = 500, \
                width = 600, violin_color=dim('index').str(), cmap='Set1', xrotation=25, xlabel = 'Genes')
    return violine

def plot_pca(anndata, color_by):
    df = pd.DataFrame(anndata.obsm['X_pca'][:, :2], columns = ['pca1', 'pca2'])
    df[color_by] = anndata.obs[color_by].to_list()
    if df[color_by].dtype == 'object':
        points = hv.Points(df).opts(color = color_by, width = 900, height = 500, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], cmap = 'set1', xaxis = None,\
                yaxis = None)
    else:
        points = hv.Points(df).opts(color = color_by, width = 900, height = 500, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], xaxis = None, yaxis = None)
    return points

def plot_tsne(anndata, color_by):
    df = pd.DataFrame(anndata.obsm['X_tsne'], columns = ['tsne1', 'tsne2'])
    df[color_by] = anndata.obs[color_by].to_list()
    if df[color_by].dtype == 'object':
        points = hv.Points(df).opts(color = color_by, width = 900, height = 500, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], cmap = 'set1', xaxis = None,\
                yaxis = None)
    else:
        points = hv.Points(df).opts(color = color_by, width = 900, height = 500, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], xaxis = None,\
                yaxis = None)
    return points

def plot_umap(anndata, color_by):
    df = pd.DataFrame(anndata.obsm['X_umap'], columns = ['umap1', 'umap2'])
    df[color_by] = anndata.obs[color_by].to_list()
    if df[color_by].dtype == 'object':
        points = hv.Points(df).opts(color = color_by, width = 900, height = 500, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], cmap = 'set1', xaxis = None,\
                yaxis = None)
    else:
        points = hv.Points(df).opts(color = color_by, width = 900, height = 500, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], xaxis = None,\
                yaxis = None)
    return points




#@st.cache(suppress_st_warning=True, persist=True)
def viz_info(adata):
    st.header('Single Cell data Vizualization')
    gex = st.container()
    cor = st.container()
    
    st.sidebar.markdown('<p style="font-family:sans-serif; font-weight:bold ;color:#FDE0D1 ;font-size: 18px;">Vizualization</p>',\
        unsafe_allow_html = True)
    st.sidebar.markdown('Gene Expression Plots')
    violine_plot = st.sidebar.button('Violine Plot')
    if True:
        st.markdown('Violine Plot of highly expressed genes')
        violine_genes = plot_violine(adata)
        st.bokeh_chart(hv.render(violine_genes, backend = 'bokeh'))

    clust_map = st.sidebar.button('Cluster Map')
    if clust_map:
        st.markdown('Hierarchical Clustering of up regulated and down regulated genes')
        st.info('Work in progress')

    dot_matrix = st.sidebar.button('Dot Matrix')
    if dot_matrix:
        st.markdown('Dot of upregulated and down regulated genes')
        st.info('Work in progress')

    st.sidebar.markdown('Single Cell Clusters')

    pca_plot = st.sidebar.button('PCA')
    if pca_plot:
        with st.spinner('Running...'):
            scanpy.pp.pca(adata)
        st.success('Done!')
        st.markdown("PCA PLOT")
        pca = plot_pca(adata, 'cell_type')
        st.bokeh_chart(hv.render(pca, backend = 'bokeh'))

    tsne_plot = st.sidebar.button('TSNE')
    if tsne_plot:
        with st.spinner('Running...'):        
            scanpy.tl.tsne(adata)
        st.success('Done!')
        tsne =  plot_tsne(adata, 'n_counts')
        st.markdown("TSNE PLOT")
        st.bokeh_chart(hv.render(tsne, backend = 'bokeh'))


    umap_plot = st.sidebar.button('UMAP')
    if umap_plot:
        with st.spinner('Running...'):
            scanpy.pp.neighbors(adata)
            scanpy.tl.umap(adata, min_dist=0.5, spread=1.0, random_state=1, n_components=2)
        st.success('Done!')
        umap =  plot_umap(adata, 'n_counts')
        st.markdown("UMAP PLOT")
        st.bokeh_chart(hv.render(umap, backend = 'bokeh'))

    heatmap = st.sidebar.button('Distance Matrix')










        
#@st.cache(suppress_st_warning=True, persist=True)
def deg_info(adata):
    st.header('Differential Gene Expression Data analysis')
    st.sidebar.markdown('<p style="font-family:sans-serif; font-weight:bold ;color:#FDE0D1 ;font-size: 18px;">Differential Gene Expression</p>',\
        unsafe_allow_html = True)
    deg_sec = st.sidebar.selectbox(
        "Identify differential genes between the cell types or with in the cell conditions",
        ("Cell Types", "Pseudo-bulk analysis")
    )
    if deg_sec == "Cell Types":
        between_cell_deg(adata)
    
    elif deg_sec == "Pseudo-bulk analysis":
        pseudo_bulk_deg(adata)







#def Null_adata():
#    _, col2, col3 = st.columns([3, 1, 1])
#
#    with open('apps/STYLE.css') as f:
#        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
#    with col2:
#        st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cells<br></p><p style="font-family:Arial; font-size: 20px;"><strong>0</strong></p>''',  unsafe_allow_html=True)
#    with col3:
#        st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Genes<br></p><p style="font-family:Arial; font-size: 20px;"><strong>0</strong>''',  unsafe_allow_html=True)
#    _, col4 = st.columns([3,2])
#    with col4:
#        st.markdown(f'''><p style="font-family:Arial; font-size: 15px;">Cell Types</p><p style="font-family:Arial; font-size: 20px;"><strong>0</strong></p>''',  unsafe_allow_html=True)