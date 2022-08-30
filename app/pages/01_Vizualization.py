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
st.header('Single Cell data Vizualization')

#def plot_violine(anndata, top_genes):
#    genes = anndata.to_df().mean(0).sort_values(ascending = False).head(top_genes).index
#    melt_df  = anndata.to_df().loc[:,genes].melt(var_name = 'index')
#    violine = hv.Violin(melt_df, ['index'], \
#            'value', label = f'Top {top_genes} Highly expressed genes').opts(height = 500, \
#                width = 800, violin_color=dim('index').str(), cmap='Set1', xrotation=25, xlabel = 'Genes')
#    return violine, melt_df

def on_change_val_vio():
    if st.session_state.val_vio:
        st.session_state.update_val_vio = st.session_state.val_vio

#@st.cache(suppress_st_warning=True, allow_output_mutation = True)
def violine_plotter(anndata):
    col1, col2 = st.columns([3,1])
    with col1:
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-size: 25px;">Violine Plot of highly expressed genes</p>',\
        unsafe_allow_html = True)
    grouper = st.container()
    with grouper:
        with col2:
            st.text_input("show",10, key = 'val_vio', on_change = on_change_val_vio)
        if 'update_val_vio' not in st.session_state:
            st.session_state.update_val_vio = st.session_state.val_vio    
        genes = anndata.to_df().mean(0).sort_values(ascending = False).head(int(st.session_state.update_val_vio)).index

        melt_df  = anndata.to_df().loc[:,genes].melt(var_name = 'index')
        with col1:
            violine = hv.Violin(melt_df, ['index'], \
                    'value', label = f'Top {st.session_state.update_val_vio} Highly expressed genes').opts(height = 500, \
                        width = 1400, violin_color=dim('index').str(), cmap='Set1', xrotation=25, xlabel = 'Genes')
    return violine, melt_df

def plot_heatmap(anndata, groupby,top_n_genes = None, markers = None):
    anndata = anndata.raw.to_adata()
    if anndata.obs[groupby].nunique() >  1:
        scanpy.tl.dendrogram(anndata, groupby)
        if markers != None and top_n_genes == None:
            anndata = anndata[:, markers]
        if top_n_genes != None and markers == None:
            genes = anndata.to_df().mean(0).sort_values(ascending = False).head(top_n_genes).index
            anndata = anndata[:, genes]
        data_frame = anndata.to_df()
        data_frame['label'] = anndata.obs[groupby]
        mean_exp = data_frame.groupby('label').agg('mean')
        mean_exp = mean_exp.loc[anndata.uns[f'dendrogram_{groupby}']['categories_ordered'],:]
        mean_exp = mean_exp.melt(ignore_index = False,var_name = 'index',value_name =  'mean_exp').reset_index()
        print(mean_exp)
        mean_exp = mean_exp.loc[:, ['index', 'label', 'mean_exp']]
        heatmap = hv.HeatMap(mean_exp, label = f'Top {top_n_genes} Highly expressed genes').opts(cmap = 'autumn', \
            height = 500, width = 800, colorbar = True, colorbar_position = 'left', tools = ['hover'],\
                line_width = 0, line_color = None, invert_yaxis=True, xrotation=45)
        return heatmap, mean_exp
    else:
        return None, None

#@st.experimental_memo(suppress_st_warning=True)
def heatmap_plotter(anndata, values, markers = None):
    anndata = anndata.raw.to_adata()
    clust_col1,clust_col2 = st.columns([2,0.5])
    with clust_col1:
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-size: 25px;">Hierarchical Clustering of Top expressed Genes</p>',\
            unsafe_allow_html = True)
    with clust_col2:
        st.text_input("show",10, key = 'val_clust', on_change = on_change_val_clust)
    if 'update_val_clust' not in st.session_state:
        st.session_state.update_val_clust = st.session_state.val_clust
    with clust_col1:
        st.selectbox('Group by', values, key = 'clust', on_change = group_on_change_clust)
    if 'update_clust' not in st.session_state:
        st.session_state.update_clust  = st.session_state.clust
    if anndata.obs[st.session_state.update_clust].nunique() >  1:
        scanpy.tl.dendrogram(anndata, st.session_state.update_clust)
        if markers != None and st.session_state.update_val_clust == None:
            anndata = anndata[:, markers]
        if st.session_state.update_val_clust != None and markers == None:
            genes = anndata.to_df().mean(0).sort_values(ascending = False).head(int(st.session_state.update_val_clust)).index
            anndata = anndata[:, genes]
        data_frame = anndata.to_df()
        data_frame['label'] = anndata.obs[st.session_state.update_clust]
        mean_exp = data_frame.groupby('label').agg('mean')
        mean_exp = mean_exp.loc[anndata.uns[f'dendrogram_{st.session_state.update_clust}']['categories_ordered'],:]
        mean_exp = mean_exp.melt(ignore_index = False,var_name = 'index',value_name =  'mean_exp').reset_index()
        print(mean_exp)
        mean_exp = mean_exp.loc[:, ['index', 'label', 'mean_exp']]
        heatmap = hv.HeatMap(mean_exp, label = f'Top {st.session_state.update_val_clust} Highly expressed genes').opts(cmap = 'autumn', \
            height = 500, width = 1400, colorbar = True, colorbar_position = 'left', tools = ['hover'],\
                line_width = 0, line_color = None, invert_yaxis=True, xrotation=45) 
        return heatmap, mean_exp
    else:
        return None, None   

def Plot_dotmatrix(anndata, groupby, top_n_genes = 10,markers = None):
    anndata = anndata.raw.to_adata()
    if anndata.obs[groupby].nunique() >  1:
        scanpy.tl.dendrogram(anndata, groupby)
        if markers != None and top_n_genes == None:
            anndata = anndata[:, markers]
        if top_n_genes != None and markers == None:
            genes = anndata.to_df().mean(0).sort_values(ascending = False).head(top_n_genes).index
            anndata = anndata[:, genes]

        data_frame = anndata.to_df()
        data_frame['label'] = anndata.obs[groupby]
        mean_exp = data_frame.groupby('label').agg('mean')
        exp_counts = anndata.to_df().applymap(lambda x : 1 if x > 0 else 0)
        exp_counts['label'] = anndata.obs[groupby]
        cell_counts = exp_counts.groupby('label').agg('sum')
        mean_exp = mean_exp.loc[anndata.uns[f'dendrogram_{groupby}']['categories_ordered'],:]
        cell_counts = cell_counts.loc[anndata.uns[f'dendrogram_{groupby}']['categories_ordered'],:]
        mean_exp = mean_exp.melt(ignore_index = False, value_name =  'mean_exp', var_name = 'index').reset_index()
        cell_counts = cell_counts.melt(ignore_index = False, value_name =  'Cell Fraction', var_name = 'index').reset_index()
        comp_df = pd.merge(left = cell_counts, right = mean_exp, left_on = ['label', 'index'], \
            right_on = ['label', 'index'], how = 'inner')
        dot_plot = hv.Points(comp_df, label = f'Top {top_n_genes} Highly expressed genes').opts(height = 500, width = 800, line_width = 2,\
                        invert_axes = True, color = 'mean_exp', cmap = 'autumn', xrotation=45,\
                        colorbar = True, size = dim('Cell Fraction') * 0.05, tools = ['hover'], invert_yaxis=True)
        return dot_plot, comp_df
    else:
        return None, None

#@st.experimental_memo(suppress_st_warning=True)
def dotmatrix_plotter(anndata, values, markers = None):
    anndata = anndata.raw.to_adata()
    dot_col1,dot_col2 = st.columns([2,0.5])
    with dot_col1:
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-size: 25px;">Dot plot of Top expressed Genes</p>',\
        unsafe_allow_html = True)
    with dot_col2:
        st.text_input("show",10, key = 'val_dot', on_change = on_change_val_dot)
    if 'update_val_dot' not in st.session_state:
        st.session_state.update_val_dot = st.session_state.val_dot
    with dot_col1:
        st.selectbox('Group by', values, key = 'dot', on_change = group_on_change_dot)
    if 'update_dot' not in st.session_state:
        st.session_state.update_dot  = st.session_state.dot
    if anndata.obs[st.session_state.update_dot].nunique() >  1:
        scanpy.tl.dendrogram(anndata, st.session_state.update_dot)
        if markers != None and st.session_state.update_val_dot == None:
            anndata = anndata[:, markers]
        if st.session_state.update_val_dot != None and markers == None:
            genes = anndata.to_df().mean(0).sort_values(ascending = False).head(int(st.session_state.update_val_dot)).index
            anndata = anndata[:, genes]

        data_frame = anndata.to_df()
        data_frame['label'] = anndata.obs[st.session_state.update_dot]
        mean_exp = data_frame.groupby('label').agg('mean')
        exp_counts = anndata.to_df().applymap(lambda x : 1 if x > 0 else 0)
        exp_counts['label'] = anndata.obs[st.session_state.update_dot]
        cell_counts = exp_counts.groupby('label').agg('sum')
        mean_exp = mean_exp.loc[anndata.uns[f'dendrogram_{st.session_state.update_dot}']['categories_ordered'],:]
        cell_counts = cell_counts.loc[anndata.uns[f'dendrogram_{st.session_state.update_dot}']['categories_ordered'],:]
        mean_exp = mean_exp.melt(ignore_index = False, value_name =  'mean_exp', var_name = 'index').reset_index()
        cell_counts = cell_counts.melt(ignore_index = False, value_name =  'Cell Fraction', var_name = 'index').reset_index()
        comp_df = pd.merge(left = cell_counts, right = mean_exp, left_on = ['label', 'index'], \
            right_on = ['label', 'index'], how = 'inner')
        dot_plot = hv.Points(comp_df, label = f'Top {st.session_state.update_val_dot} Highly expressed genes').opts(height = 500, width = 1400, line_width = 2,\
                        invert_axes = True, color = 'mean_exp', cmap = 'autumn', xrotation=45,\
                        colorbar = True, size = dim('Cell Fraction') * 0.05, tools = ['hover'], invert_yaxis=True)
        return dot_plot, comp_df
    else:
        return None, None

@st.cache(suppress_st_warning=True, allow_output_mutation = True)
def plot_pca(anndata, color_by):
    df = pd.DataFrame(anndata.obsm['X_pca'][:, :2], columns = ['pca1', 'pca2'])
    df[color_by] = anndata.obs[color_by].to_list()
    if df[color_by].dtype == 'object':
        points = hv.Points(df).opts(color = color_by, width = 1400, height = 800, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], cmap = 'set1', xaxis = None,\
                yaxis = None)
    else:
        points = hv.Points(df).opts(color = color_by, width = 1400, height = 800, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], xaxis = None, yaxis = None, colorbar = True)
    return points
@st.cache(suppress_st_warning=True, allow_output_mutation = True)
def plot_tsne(anndata, color_by):
    df = pd.DataFrame(anndata.obsm['X_tsne'], columns = ['tsne1', 'tsne2'])
    df[color_by] = anndata.obs[color_by].to_list()
    if df[color_by].dtype == 'object':
        points = hv.Points(df).opts(color = color_by, width = 1400, height = 800, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], cmap = 'set1', xaxis = None,\
                yaxis = None)
    else:
        points = hv.Points(df).opts(color = color_by, width = 1400, height = 800, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], xaxis = None,\
                yaxis = None, colorbar = True)
    return points
@st.cache(suppress_st_warning=True, allow_output_mutation = True)
def plot_umap(anndata, color_by):

    df = pd.DataFrame(anndata.obsm['X_umap'], columns = ['umap1', 'umap2'])
    df[color_by] = anndata.obs[color_by].to_list()
    if df[color_by].dtype == 'object':
        points = hv.Points(df).opts(color = color_by, width = 1400, height = 800, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], cmap = 'set1', xaxis = None,\
                yaxis = None)
    else:
        points = hv.Points(df).opts(color = color_by, width = 1400, height = 800, \
            legend_position = 'right', legend_offset = (0,300), tools = ['hover'], xaxis = None,\
                yaxis = None, colorbar = True)
    return points

def group_on_change_clust():
    if st.session_state.clust:
        st.session_state.update_clust = st.session_state.clust
def group_on_change_dot():
    if st.session_state.dot:
        st.session_state.update_dot = st.session_state.dot


def on_change_val_clust():
    if st.session_state.val_clust:
        st.session_state.update_val_clust = st.session_state.val_clust

def on_change_val_dot():
    if st.session_state.val_dot:
        st.session_state.update_val_dot = st.session_state.val_dot

def group_on_change_pca():
    if st.session_state.pca:
        st.session_state.update_pca = st.session_state.pca

def group_on_change_tsne():
    if st.session_state.tsne:
        st.session_state.update_tsne = st.session_state.tsne

def group_on_change_umap():
    if st.session_state.umap:
        st.session_state.update_umap = st.session_state.umap



def viz_info(adata):
    gex = st.container()
    cor = st.container()
    st.sidebar.markdown('<p style="font-family:sans-serif; font-weight:bold ;color:#ff9633;font-size: 18px;">Vizualization</p>',\
        unsafe_allow_html = True)
    st.sidebar.markdown('Gene Expression Plots')
    violine_plot = st.sidebar.checkbox('Show Violine Plot')
    if violine_plot:
        violine_genes, violine_df = violine_plotter(adata)
        vio1, vio2 = st.tabs(['Plot','Table'])
        with vio1:
            st.bokeh_chart(hv.render(violine_genes, backend = 'bokeh'))
        with vio2:
            st.bokeh_chart(hv.render(hv.Table(violine_df).opts(width = 800, height = 500), backend = 'bokeh'))


    obs_dtypes = adata.obs.dtypes.to_frame(name='type')
    def get_count(include_only = 'category'):
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

    
    bx_val = get_count().keys() 
    all_val = get_count('all').keys()

    clust_map = st.sidebar.checkbox('Show Cluster Map')
    if clust_map:
        heatmap, heatmap_df = heatmap_plotter(adata, bx_val)
        if heatmap == None and heatmap_df == None:
            st.warning('There are 1 unique value in selected group')
        else:
            clust1, clust2 = st.tabs(['Plot', 'Table'])
            with clust1:
                st.bokeh_chart(hv.render(heatmap, backend = 'bokeh'))  
            with clust2:
                st.bokeh_chart(hv.render(hv.Table(heatmap_df).opts(width = 800, height = 500), backend = 'bokeh')) 

    dot_matrix = st.sidebar.checkbox('Show Dot Matrix')
    if dot_matrix:
        dotplot, dotplot_df = dotmatrix_plotter(adata,bx_val)
        if dotplot == None and dotplot_df == None:
            st.warning('There are 1 unique value in selected group')
        else:
            dot1, dot2 = st.tabs(['Plot', 'Table'])
            with dot1:
                st.bokeh_chart(hv.render(dotplot, backend = 'bokeh'))
            with dot2:
                st.bokeh_chart(hv.render(hv.Table(dotplot_df).opts(width = 800, height = 500), backend = 'bokeh'))    

    st.sidebar.markdown('Single Cell Clusters')

    show_cluster = st.sidebar.radio("Show", ('PCA','TSNE','UMAP'))

    if show_cluster == 'PCA':

        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-weight:bold; font-size: 18px;">PCA PLOT</p>',unsafe_allow_html = True)
        st.selectbox('Group by', all_val, key = 'pca', on_change = group_on_change_pca)
        if 'update_pca' not in st.session_state:
            st.session_state.update_pca = st.session_state.pca
        pca = plot_pca(adata, st.session_state.update_pca)
        st.bokeh_chart(hv.render(pca, backend = 'bokeh'))


    elif show_cluster == 'TSNE':

        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-weight:bold; font-size: 18px;">TSNE PLOT</p>',unsafe_allow_html = True)
        st.selectbox('Group by', all_val, key = 'tsne', on_change = group_on_change_tsne)
        if 'update_tsne' not in st.session_state:
            st.session_state.update_tsne = st.session_state.tsne

        tsne = plot_tsne(adata, st.session_state.update_tsne)
        st.bokeh_chart(hv.render(tsne, backend = 'bokeh'))

    else:

        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-weight:bold; font-size: 18px;">UMAP PLOT</p>',unsafe_allow_html = True)
        st.selectbox('Group by', all_val, key = 'umap', on_change = group_on_change_umap)
        if 'update_umap' not in st.session_state:
            st.session_state.update_umap = st.session_state.umap
        umap = plot_umap(adata, st.session_state.update_umap)
        st.bokeh_chart(hv.render(umap, backend = 'bokeh'))


    heatmap = st.sidebar.button('Distance Matrix')

#adata, upfile = file_uploader()

#if upfile:
#    viz_info(adata)

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
    viz_info(st.session_state["adata_obj"])
else:
    st.write("upload the file")