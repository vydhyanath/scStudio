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
    st.markdown('<p style="font-family:Arial; color:#CA0327; font-weight:bold; font-size:20px">Between Cell Groups</p>', unsafe_allow_html=True)
    st.info('Work in Progress')
def pseudo_bulk_deg(adata):
    st.markdown('<p style="font-family:Arial; color:#CA0327; font-weight:bold; font-size:20px">Between Sample Groups</p>', unsafe_allow_html=True)
    st.info('Work in Progress')

#@st.cache(suppress_st_warning=True, persist=True)
def meta_info(adata):
    st.header('General Meta data')
    adata_stats = lambda x : adata.obs[x].nunique() \
        if adata.obs[x].dtype == 'category' else adata.obs[x].sum() 
    TABLE = st.container()


    tab1, tab2 = st.tabs(['Plot','Table'])
#    col1, col2, col3= st.columns([2.5, 0.5, 0.5])

    meta_df = pd.DataFrame({col : adata_stats(col) for col in adata.obs.columns}.items(), \
            columns = ['Fields', 'Unique Values'])
    table_plot = hv.Table(meta_df)

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
#    print(obs_dtypes)

    def get_count():
        cols = [col for col in obs_dtypes.index if str(obs_dtypes.loc[col].type) == 'category' \
            if adata.obs[col].nunique() != adata.shape[0]]
        obs_dic = dict(zip(cols, \
            [Counter(adata.obs[col]) for col in cols]))
        return obs_dic

    print(get_count())        
    obs_dict = get_count()
#    obs_dict = dict(zip(obs_dtypes.index, \
#        [Counter(adata.obs[col]) for col in obs_dtypes.index \
#            if str(obs_dtypes.loc[col].type) == 'category']))
#    for key,value in obs_dict.items():
#        if value == len(adata.obs_names):
#    del obs_dict['barcodes']
#    print(obs_dict.keys())
        

#    mp_bar = st.container() 
#    with mp_bar:
#        st.markdown('<p style="font-family:sans-serif; font-size: 18px;">general statistics of single cells observation</p>', \
#                unsafe_allow_html = True)

    acol1, _,acol3 = st.columns([1.5, 0.1, 0.35])
    with acol3:
        with tab1:
            st.selectbox('Variable to plot',obs_dict.keys(), key = 'col', on_change = on_select)
            if 'update_column' not in st.session_state:
                st.session_state.update_column = st.session_state.col
    with acol1:
#        st.altair_chart(plot_bar(obs_dict, st.session_state.update_column), use_container_width= True)
        bar_plot = hv.Bars(data = adata.obs[st.session_state.update_column].value_counts()).opts(height = 500, \
            width = 700, invert_axes = True, tools = ['hover'], xticks = 10)
        with tab1:
            st.bokeh_chart(hv.render(bar_plot, backend = 'bokeh'))
        with tab2:
            st.bokeh_chart(hv.render(table_plot, backend = 'bokeh'))



#def  qc_info(adata):
#    st.header('Quality control of single cell Sequencing data')
#    st.sidebar.markdown('<p style="font-family:sans-serif; font-weight:bold ;color:#04410D ;font-size: 18px;">Quality Control</p>',\
#        unsafe_allow_html = True)
#    st.sidebar.markdown('<p style="font-family:sans-serif ;font-size: 12px;">Quality Statistics for Single Cell Data</p>',\
#        unsafe_allow_html = True)
#    
#        
#    hist_counts = hv.Histogram((np.histogram(adata.obs['log1p_total_counts'], bins=100))).opts(width = 400, \
#        height = 300, xlabel = 'counts depth').relabel('Histograms of count depth per cell')
#    hist_genes = hv.Histogram((np.histogram(adata.obs['log1p_n_genes_by_counts'], bins=100))).opts(width = 400, \
#        height = 300, xlabel = 'genes count').relabel('Histogram of the number of genes detected per cell')
#    elb_plt = hv.Points(list(enumerate(adata.obs['n_counts'].sort_values(ascending = False).values))).\
#        opts(logy = True, width = 400, height = 300, xlabel = 'Barcode Rank', ylabel = 'log counts depth').\
#            relabel('Count depth distribution from high to low count depths')
#    genes_counts = hv.Points((adata.obs['n_counts'].values, \
#        adata.obs['n_genes'].values)).opts(xlabel = 'count depth', \
#            ylabel = 'genes per cell', width = 400, height = 300).relabel('Number of genes versus the count depth')
#
#    qc_layout = hv.Layout(hist_counts + hist_genes + elb_plt + \
#        genes_counts).opts(shared_axes = False).cols(2)
#
#    st.bokeh_chart(hv.render(qc_layout, backend = 'bokeh'))

def plot_violine(anndata, top_genes):
    genes = anndata.to_df().mean(0).sort_values(ascending = False).head(top_genes).index
    melt_df  = anndata.to_df().loc[:,genes].melt(var_name = 'index')
    violine = hv.Violin(melt_df, ['index'], \
            'value', label = f'Top {top_genes} Highly expressed genes').opts(height = 500, \
                width = 800, violin_color=dim('index').str(), cmap='Set1', xrotation=25, xlabel = 'Genes')
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

def group_on_change_clust():
    if st.session_state.clust:
        st.session_state.update_clust = st.session_state.clust
def group_on_change_dot():
    if st.session_state.dot:
        st.session_state.update_dot = st.session_state.dot

def on_change_val_vio():
    if st.session_state.val_vio:
        st.session_state.update_val_vio = st.session_state.val_vio

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
#    @st.experimental_singleton(suppress_st_warning=True)
#    def preprocess():
#        with st.spinner('Running...PCA'):
#            scanpy.pp.pca(adata)
#            st.success('Done! PCA')  
#        with st.spinner('Running...TSNE'):
#            scanpy.tl.tsne(adata)
#            st.success('Done! TSNE')
#        with st.spinner('Running...UMAP'):
#            scanpy.pp.neighbors(adata)
#            scanpy.tl.umap(adata, min_dist=0.5, spread=1.0, random_state=1, n_components=2)   
#            st.success('Done! UMAP')
#        return adata

#    adata = preprocess()   

    st.header('Single Cell data Vizualization')
    gex = st.container()
    cor = st.container()
#    print(adata)
    st.sidebar.markdown('<p style="font-family:sans-serif; font-weight:bold ;color:#04410D;font-size: 18px;">Vizualization</p>',\
        unsafe_allow_html = True)
    st.sidebar.markdown('Gene Expression Plots')
    violine_plot = st.sidebar.checkbox('Show Violine Plot')
    if violine_plot:
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-size: 18px;">Violine Plot of highly expressed genes</p>',\
        unsafe_allow_html = True)
        st.text_input("show",10, key = 'val_vio', on_change = on_change_val_vio)
        if 'update_val_vio' not in st.session_state:
            st.session_state.update_val_vio = st.session_state.val_vio
        violine_genes, violine_df = plot_violine(adata, int(st.session_state.update_val_vio))
        #print(violine_df)
        vio1, vio2 = st.tabs(['Plot','Table'])
        with vio1:
            st.bokeh_chart(hv.render(violine_genes, backend = 'bokeh'))
        with vio2:
            st.bokeh_chart(hv.render(hv.Table(violine_df).opts(width = 800, height = 500), backend = 'bokeh'))

    obs_dtypes = adata.obs.dtypes.to_frame(name='type')
    def get_count():
        cols = [col for col in obs_dtypes.index if str(obs_dtypes.loc[col].type) == 'category' \
            if adata.obs[col].nunique() != adata.shape[0]]
        obs_dic = dict(zip(cols, \
            [Counter(adata.obs[col]) for col in cols]))
        return obs_dic

    
    bx_val = get_count().keys() 

    clust_map = st.sidebar.checkbox('Show Cluster Map')
    if clust_map:
        clust_col1,clust_col2 = st.columns([2,0.5])
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-size: 18px;">Hierarchical Clustering of Top expressed Genes</p>',\
        unsafe_allow_html = True)
        with clust_col2:
            st.text_input("show",10, key = 'val_clust', on_change = on_change_val_clust)
        if 'update_val_clust' not in st.session_state:
            st.session_state.update_val_clust = st.session_state.val_clust
        with clust_col1:
            st.selectbox('Group by', bx_val, key = 'clust', on_change = group_on_change_clust)
        if 'update_clust' not in st.session_state:
            st.session_state.update_clust  = st.session_state.clust
        
        heatmap, heatmap_df = plot_heatmap(adata, st.session_state.update_clust, \
            top_n_genes=int(st.session_state.update_val_clust))
        if heatmap == None and heatmap_df == None:
            st.warning('There are 1 unique value in selected group')
        else:
            clust1, clust2 = st.tabs(['Plot', 'Table'])
            with clust1:
                st.bokeh_chart(hv.render(heatmap, backend = 'bokeh'))  
            with clust2:
                st.bokeh_chart(hv.render(hv.Table(heatmap_df).opts(width = 800, height = 500), backend = 'bokeh')) 
 #       st.info('Work in progress')

    dot_matrix = st.sidebar.checkbox('Show Dot Matrix')
    if dot_matrix:
        dot_col1,dot_col2 = st.columns([2,0.5])
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-size: 18px;">Dot plot of Top expressed Genes</p>',\
        unsafe_allow_html = True)
        with dot_col2:
            st.text_input("show",10, key = 'val_dot', on_change = on_change_val_dot)
        if 'update_val_dot' not in st.session_state:
            st.session_state.update_val_dot = st.session_state.val_dot
        with dot_col1:
            st.selectbox('Group by', bx_val, key = 'dot', on_change = group_on_change_dot)
        if 'update_dot' not in st.session_state:
            st.session_state.update_dot  = st.session_state.dot
        dotplot, dotplot_df = Plot_dotmatrix(adata, st.session_state.update_dot, \
            top_n_genes=int(st.session_state.update_val_dot))
        if dotplot == None and dotplot_df == None:
            st.warning('There are 1 unique value in selected group')
        else:
            dot1, dot2 = st.tabs(['Plot', 'Table'])
            with dot1:
                st.bokeh_chart(hv.render(dotplot, backend = 'bokeh'))
            with dot2:
                st.bokeh_chart(hv.render(hv.Table(dotplot_df).opts(width = 800, height = 500), backend = 'bokeh'))    

    st.sidebar.markdown('Single Cell Clusters')

#    pca_plot = st.sidebar.button('PCA')
    show_cluster = st.sidebar.radio("Show", ('PCA','TSNE','UMAP'))

    if show_cluster == 'PCA':
#        with st.spinner('Running...'):
#            scanpy.pp.pca(adata)
#        st.success('Done!')
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-weight:bold; font-size: 18px;">PCA PLOT</p>',unsafe_allow_html = True)
        st.selectbox('Group by', bx_val, key = 'pca', on_change = group_on_change_pca)
        if 'update_pca' not in st.session_state:
            st.session_state.update_pca = st.session_state.pca
        pca = plot_pca(adata, st.session_state.update_pca)
        st.bokeh_chart(hv.render(pca, backend = 'bokeh'))

#    tsne_plot = st.sidebar.button('TSNE')
    elif show_cluster == 'TSNE':
#        with st.spinner('Running...'):        
#            scanpy.tl.tsne(adata)
#        st.success('Done!')
#        tsne =  plot_tsne(adata, 'n_counts')
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-weight:bold; font-size: 18px;">TSNE PLOT</p>',unsafe_allow_html = True)
        st.selectbox('Group by', bx_val, key = 'tsne', on_change = group_on_change_tsne)
        if 'update_tsne' not in st.session_state:
            st.session_state.update_tsne = st.session_state.tsne
#        preprocess_TSNE()
        tsne = plot_tsne(adata, st.session_state.update_tsne)
        st.bokeh_chart(hv.render(tsne, backend = 'bokeh'))

#    umap_plot = st.sidebar.button('UMAP')
    else:
#        with st.spinner('Running...'):
#            scanpy.pp.neighbors(adata)
#            scanpy.tl.umap(adata, min_dist=0.5, spread=1.0, random_state=1, n_components=2)
#        st.success('Done!')
#        umap =  plot_umap(adata, 'n_counts')
        st.markdown('<p style="font-family:sans-serif; color:#CA0327; font-weight:bold; font-size: 18px;">UMAP PLOT</p>',unsafe_allow_html = True)
        st.selectbox('Group by', bx_val, key = 'umap', on_change = group_on_change_umap)
        if 'update_umap' not in st.session_state:
            st.session_state.update_umap = st.session_state.umap
#        preprocess_UMAP()
        umap = plot_umap(adata, st.session_state.update_umap)
        st.bokeh_chart(hv.render(umap, backend = 'bokeh'))


    heatmap = st.sidebar.button('Distance Matrix')










        
#@st.cache(suppress_st_warning=True, persist=True)
def deg_info(adata):
    st.header('Differential Gene Expression Data analysis')
    st.sidebar.markdown('<p style="font-family:sans-serif; font-weight:bold ;color:#04410D ;font-size: 18px;">Differential Gene Expression</p>',\
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