import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import Descriptors, PandasTools
import cirpy 
from cirpy import Molecule
import plotly.graph_objs as go
import dash_bootstrap_components as dbc 
from jupyter_dash import JupyterDash  
import base64
import datetime
import io
import plotly.graph_objs as go
from dash_extensions import Lottie

import base64
from io import BytesIO


url_coonections = "https://assets2.lottiefiles.com/packages/lf20_rD5POs.json"
options = dict(loop=True, autoplay=True, rendererSettings=dict(preserveAspectRatio='xMidYMid slice'))



image_component = html.Img(id = 'structure-image')
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


df_with_smiles = pd.read_csv(r'D:\df_with_smiles.csv',encoding='ISO-8859-1', delimiter=';')
df_with_smiles['Factor'] = df_with_smiles.Area * df_with_smiles.Fuerza
tipo_olor_total = df_with_smiles.groupby('Tipo_olor').agg(frecuency = ('Tipo_olor', 'count')).reset_index()
suma_factor = df_with_smiles.groupby('Tipo_olor').agg(suma_factor = ('Factor', 'sum')).reset_index()
df_with_smiles_total_sum = pd.concat([tipo_olor_total.iloc[:,0:2], suma_factor.iloc[:,1:2].round(3)], axis =1)

variables = df_with_smiles.columns
variables3 = df_with_smiles_total_sum.columns

app.layout = html.Div([
    dbc.Container([
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                    html.H1('Análisis de señales cromatográficas en fragancias'),
                    dbc.CardHeader(Lottie(options=options, width="20%", height="20%", url=url_coonections)),

                ],  style={'textAlign':'center'})
                    ],color="info"),

        ],width=8),
    ])
    ],className ='mb-2 mr-1'),

#####################################
    html.Div([
        dbc.CardBody([
            html.H6('Seleccione variables del eje Y')
        ],style={'textAlign':'left'}),
        dcc.Dropdown(
            id ='ejex',
            options = [{'label':i, 'value':i} for i in variables.unique()],
            value = 'Presion_vapor'
        )

                ], style = {'width': '48%', 'display': 'inline-block'}),

    html.Div([
    dbc.CardBody([
        html.H6('Seleccione variables del eje X')
    ],style={'textAlign':'right'}),
        dcc.Dropdown(
            id = 'ejey',
            options = [{'label':i, 'value':i} for i in variables],
            value = 'Area')
                ],style = {'width': '48%', 'float':'right', 'display': 'inline-block'}),

    html.Div([
        dcc.Graph(id = 'Grafico_var')],style={'width': '100%'}), # Para identificar el grafico

    html.Div([image_component]),

###################################################################################################################################################################################

#### HTML DIV PARA LOS GRAFICOS DE TIPO DE DISCRIPTOR POR TIPO DE OLOR
    html.Div([

            dcc.Dropdown(
                id = 'ejeyyy',
                options = [{'label':'Tipo_olor', 'value':'Tipo_olor'}],
                value = 'Tipo_olor')
                ],style = {'width': '10%'}),

#### HTML DIV PARA LOS GRAFICOS DE TIPO DE DISCRIPTOR
    html.Div([
            dcc.Graph(id='smiles_total_sum')
            ],style={'width': '60%', 'display': 'inline-block'}),

################## interactividad grafico de barras

    html.Div([
        dcc.Graph(id='barplot_descriptor')
        ],style={'width': '40%', 'display':'inline-block'}),

################# Mi firma
        dbc.Container([
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardImg(src='/assets/firmaJDM.png')
                        ],style={'textAlign':'right'} ),
            ],width=8),
        ])
        ],className ='mr-1'),


], style = {'padding':10}) 

    ###################################################################################################################################################################################


## Creacion del grafico e interactividad
# El callback para actualizar grafico en funcion de los 2 dropdown

@app.callback(
    Output('Grafico_var', 'figure'),
    [Input('ejex', 'value'),
    Input('ejey', 'value')]
)
def actualizar_graf(nombre_ejex, nombre_ejey):

        fig = px.scatter(df_with_smiles,
                    x = nombre_ejex,
                    y =  nombre_ejey,
                    color = df_with_smiles['Tipo_olor'],
                    size = df_with_smiles['Fuerza'],
                    hover_data=variables[:-1],

                    #log_x=True,
                    size_max=20,
                    opacity =1)
        
        return fig

### Definir las imagenes de las moleculas
@ app.callback(
    Output('structure-image', 'src'),
    [Input('Grafico_var', 'selectedData')])


def display_selected_data(selectedData):

    max_structs = 18
    structs_per_row = 6
    empty_plot = "data:image/gif;base64,R0lGODlhAQABAAAAACwAAAAAAQABAAA="

    if selectedData:
        if len(selectedData['points']) == 0:
            return empty_plot
        match_idx = [x['x'] for x in selectedData['points']]
        smiles_list = [Chem.MolFromSmiles(x) for x in list(df_with_smiles.iloc[match_idx].smiles)]
        name_list = list(df_with_smiles.iloc[match_idx].nomenclature)

        #active_list = list(df_with_smiles.iloc[match_idx].Tipo_olor)
        #name_list = [x + ' ' + str(y) for (x,y) in zip(name_list)]
        img = MolsToGridImage(smiles_list[0:max_structs], molsPerRow = structs_per_row, legends = name_list)
        buffered = BytesIO()
        img.save(buffered, format="JPEG")
        encoded_image = base64.b64encode(buffered.getvalue())
        src_str = 'data:image/png;base64,{}'.format(encoded_image.decode())
    else:
        return empty_plot
    return src_str

###################################################################################################################################################################################



@app.callback(
    Output('smiles_total_sum', 'figure'),
    [Input('ejeyyy', 'value')]
)

def actualizar_graf(nombre_ejeyyy):

    fig = px.scatter(df_with_smiles_total_sum, x = nombre_ejeyyy,
                            y = 'frecuency',
                            color = df_with_smiles_total_sum['Tipo_olor'],
                            size = 'suma_factor',
                            hover_data= variables3 ,
                            size_max=30,
                            opacity =0.9)
    fig.update_layout(hovermode='x')
    fig.update_layout(title ='Total tipos de olores', yaxis_title="Frecuency")

    return fig

@app.callback(
    Output('barplot_descriptor', 'figure'),
    [Input('smiles_total_sum', 'clickData')]) 

def actualizar_graph_cat(clickData):
    v_index = clickData['points'][0]['x']

    descriptor = df_with_smiles.groupby('Tipo_olor').get_group(v_index)

    fig = go.Figure()
    fig.add_trace(go.Histogram(histfunc="count", x=descriptor['Descriptor_1'], y=descriptor['Tipo_olor'], name="Descriptor_1"))
    fig.add_trace(go.Histogram(histfunc="count", x=descriptor['Descriptor_2'], y=descriptor['Tipo_olor'], name='Descriptor_2'))
    fig.update_layout(barmode='stack')
    fig.update_layout(hovermode='x')
    fig.update_layout(title = v_index, yaxis_title="Frecuency")


    return fig



if __name__ == '__main__':
    app.run_server()
