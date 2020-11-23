""" LOAD MODULES AND PROCES DATA """
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

import numpy as np
import pandas as pd
import pathlib
import pickle
import dash_katex
import time
import seaborn as sn
import flask
import os
import dash_auth

import Figures
import LinearRegression

#remove regression warnings
np.seterr(all='ignore')

""" TEMPORARY PASSWORD """
#VALID_USERNAME_PASSWORD_PAIRS = [ ['itz', 'itz'] ]
server = flask.Flask(__name__)
server.secret_key = os.environ.get('secret_key', 'secret')

""" CREATE APP """
app = dash.Dash(
    __name__,
    server = server,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    url_base_pathname='/gunicorn/',
)

app.config.suppress_callback_exceptions = True
#auth = dash_auth.BasicAuth(app, VALID_USERNAME_PASSWORD_PAIRS)

""" LOAD DATA """
# Path
BASE_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()

#load gene data
#gene_data = shelve.open("data/gene_data")
gene_data = pickle.load( open( "data/gene_data_dic.p", "rb" ))
gene_list = sorted(gene_data.keys())

#load atger data (already filtered to contain only common genes)
gene_data_atg = pickle.load( open( "data/gene_data_atg_dic.p", "rb" ))
gene_list_atg = sorted(gene_data_atg.keys())

#load 2D regressions
dic_reg = pickle.load( open( "data/dic_reg.p", "rb" ))

#load colors
colorscale = sn.color_palette("GnBu_d",8).as_hex()

""" CREATE TAB STYLE """
tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    #'padding': '18px',
    'fontWeight': 'normal',
    #'font-size' : 'large'
}

tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': colorscale[3],
    'color': 'white',
    'padding': '18px',
    #'font-size' : 'large'
}

tab_style_reg = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '10px',
    'fontWeight': 'normal',
}

tab_selected_style_reg = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': colorscale[3],
    'padding': '10px',
    'color' : 'white'
}

""" CREATE IMPORTANT APP PARTS """
def description_card():
    """
    :return: A Div containing dashboard title & descriptions for tab 1.
    """
    return html.Div(
        id="description-card",
        children = [
            html.H3("Welcome to CZviz"),
            dcc.Markdown('''This app allows you to explore the full dataset and \
                          analysis from the study by _Droin & al_ published in \
                          _Nature Metabolism_.'''),
        ],
    )

def description_card2():
    """
    :return: A Div containing dashboard title & descriptions for tab 2.
    """
    return html.Div(
            id="intro2",
            children = [
                html.Div(id = 'paragraph-selected-gene-2', style = {'text-align':'left'}),
                #dcc.Markdown("Please select the type of analysis you're interested in.")
                ],
    )


def generate_control_card():
    """
    :return: A Div containing controls for graphs.
    """
    return html.Div(
        id="control-card",
        children=[
            html.P("Please type gene name"),
            dcc.Dropdown(
                id="gene-select",
                options=[{"label": gene_name.capitalize(), "value": gene_name} for gene_name in gene_list],
                value='cry1',
            ),
            html.Br(),
        ],
    )

def generate_control_card2():
    """
    :return: A Div containing controls for analysis.
    """
    return html.Div(
        id="control-card2",
        children=[
            #html.P("Please type gene name"),
            #dcc.Dropdown(
            #    id="gene-select2",
            #    options=[{"label": gene_name.capitalize(), "value": gene_name} for gene_name in gene_list],
            #    value='pecr',
            #),
            #html.Br(),
            html.P("Please select analysis"),
            dcc.Dropdown(
                id="analysis-select",
                options=[{"label": "Spatial analysis", "value": 'spatial'},
                         {"label": "Temporal analysis", "value": 'temporal'},
                         {"label": "Spatiotemporal analysis", "value": 'spatiotemporal'},
                         {"label": "Validation", "value": 'validation'}],
                value='spatial',
            ),
            html.Br(),
        ],
    )



def generate_analysis_card1():
    """
    :return: A Div containing analysis 1.
    """
    return html.Div(
        id="analysis-card1",
        children=[
            html.H4('Zonation: analysis timepoint per timepoint', style = {'textAlign' : 'center'}),
            dcc.Tabs(id = 'tabs-time',
                children = [
                    dcc.Tab(
                        label ='t = ' + str((i-1)*6) + 'h',
                        value = str(i),
                        style = tab_style_reg,
                        selected_style = tab_selected_style_reg
                        )
                        for i in range(1,5)
                    ],
                    value = '1',
                    vertical = False,
                    style = {'background-color': 'Grey '}
                ),
            dcc.Loading(
                id="sub-loading-1",
                type="circle",
                color = colorscale[3],
                children=[
                    html.Div(
                        id = 'summary-space-regression',
                        children = [
                            dcc.Graph(
                                id ='graph-stat-space',
                                config = {
                                    'displayModeBar': False,
                                    'modeBarButtonsToRemove': [],
                                    'displaylogo' : False
                                    }
                                )
                            ],
                        ),
                    ]
                )
            ]
        )

def generate_analysis_card2():
    """
    :return: A Div containing analysis 2.
    """
    return html.Div(
        id="analysis-card2",
        children=[
            html.Div(
                id = 'initial-title-tab3',
                children = [
                    html.H4(
                        children = 'Rhythmicity: analysis layer per layer',
                        style = {'textAlign' : 'center'}
                        ),
                    ]
                ),
            html.Div(
                children = [
                    dcc.Graph(
                        id = 'graph-polar',
                        style = {
                            'width' : '24vw',
                            'zIndex': '4'
                            },
                        config = {
                            'displayModeBar' : False,
                            'modeBarButtonsToRemove': [],
                            'displaylogo' : False
                            }
                        ),
                    dcc.Graph(
                        id = 'graph-mean',
                        style = {
                            'width' : '32vw',
                            'height' : '40vh'
                            },
                        config={
                            'displayModeBar': False,
                            'modeBarButtonsToRemove': [],
                            'displaylogo' : False
                            }
                        ),
                    ],
                style = {
                    'display' : 'flex',
                    'flex-direction': 'rows',
                    'flex-wrap' : 'wrap',
                    'align-items' : 'flex-start',
                    'justify-content' : 'space-between'
                    }
                ),
            html.Div(
                id = 'second-title-tab3',
                children = [
                    dcc.Tabs(
                        id = 'tabs-space',
                        children = [
                            dcc.Tab(
                                label='x = ' + str(i-1),
                                value = str(i),
                                style = tab_style_reg,
                                selected_style=tab_selected_style_reg) for i in range(1,9)
                            ],
                        value = '1',
                        vertical = False,
                        style = {'background-color': 'Grey '}
                    ),                                            ],
                style = {'zIndex': '10'}
                ),
            html.Div(
                style = {'height': '40vw'}, #to avoid skipping when changing tab
                children = [
                    dcc.Loading(
                        id="sub-loading-2",
                        type="circle",
                        color = colorscale[3],
                        children=[
                            html.Div(
                                id = 'summary-time-regression',
                                children =
                                    dcc.Graph(
                                        id='graph-stat-time',
                                        config={
                                            'displayModeBar': False,
                                            'modeBarButtonsToRemove': [],
                                            'displaylogo' : False
                                            }
                                        ),
                                ),
                            ],
                        ),
                    ]
                )
            ]
        )


def generate_analysis_card3():
    """
    :return: A Div containing analysis 3.
    """
    return html.Div(
        id="analysis-card3",
        children = [
            html.Div(
                id = 'fourth-title-tab3',
                children = [
                    html.H4(
                        children = 'Rythmic zonation: analysis of all datapoints together',
                        style = {'textAlign' : 'center'}
                        ),
                    ]
                ),
            html.Div(
                id = 'summary-2D-regression',
                children = [
                    dcc.Graph(id='graph-polar-fit'),
                    dcc.Graph(id='graph-fit-3D')
                    ],
                ),
            html.Div(
                style = {'margin' : '0px 0px 300px 0px'}
                )
            ],
            #style = {'display':'none'},
        )

def generate_analysis_card4():
    """
    :return: A Div containing analysis 4.
    """
    return html.Div(
        id="analysis-card4",
        children = [
            html.Div(
                id = 'title-validation',
                children = [
                    html.H4(
                        children = 'Comparison with the dataset from Atger & al.',
                        style = {'textAlign' : 'center'}
                        ),
                    ]
                ),
            html.Div(children = [
                    dcc.Graph(
                        id='graph-comparison',
                        style = {
                            'width' : '40vw',
                            #'height' : '40vh'
                            },
                        config={'displayModeBar': False,
                                'modeBarButtonsToRemove': [],
                                'displaylogo' : False }
                                ),

                    dcc.Graph(
                        id='graph-polar-comparison',
                        style = {'width' : '25vw'},
                        config={'displayModeBar': False,
                                'modeBarButtonsToRemove': [],
                                'displaylogo' : False }
                        )
                    ],
                style = {
                    'display' : 'flex',
                    'flex-direction': 'rows',
                    'flex-wrap' : 'wrap',
                    'align-items' : 'flex-start',
                    'justify-content' : 'space-between'
                    }
                ),
            html.Div(id = 'no-gene', style = {'text-align': 'center'}),
            ],
        )


""" APP MAIN LAYOUT """
app.layout = html.Div(
    id="app-container",
    children=[
        dcc.Tabs(
            id = "tabs",
            value = 'main-tab-1',
            children=[
                dcc.Tab(
                    label='Gene selection',
                    value = 'main-tab-1',
                    style=tab_style, selected_style=tab_selected_style,
                    children=[
                        # Left column
                        html.Div(
                            id="left-column",
                            className="two columns",
                            children=[description_card(), generate_control_card()]
                        ),
                        #html.H4(id='wait-text', children = 'Please wait while the data is being processed...'),
                        dcc.Loading(
                            id="loading1",
                            type="circle",
                            color = colorscale[3],
                            children=[
                        # Middle column
                        html.Div(
                            id="middle-column",
                            className="five columns",
                            children=[
                                html.Div(
                                    id="graph_tab_1_card",
                                    children = dcc.Graph(
                                        id='graph-space',
                                        config={'displayModeBar': False,
                                                'modeBarButtonsToRemove': [],
                                                'displaylogo' : False
                                                },
                                        #style = {'margin': 'auto'},
                                            ),
                                ),
                                html.Div(
                                    id="graph_tab_3_card",
                                    children = dcc.Graph(
                                        id='graph-3d',
                                        config={
                                            'modeBarButtonsToRemove': [
                                                'sendDataToCloud',
                                                'resetCameraLastSave3d',
                                                'hoverClosest3d',
                                                'zoom3d',
                                                'toImage'
                                                ],
                                            'displaylogo' : False,
                                            'scrollZoom': False
                                            },
                                        style = {
                                            'border-width':'1px',
                                            'border-style':'solid',
                                            'border-color':'#e8e8e8'
                                            }
                                        ),
                                ),
                            ],
                        ),

                        # Right column
                        html.Div(
                            id="right-column",
                            className="five columns",
                            children=[
                                html.Div(
                                    id="graph_tab_2_card",
                                    children = dcc.Graph(
                                        id='graph-time',
                                        config={'displayModeBar': False,
                                                'modeBarButtonsToRemove': [],
                                                'displaylogo' : False
                                                },
                                        #style = {'margin': 'auto'},
                                            ),
                                ),
                                html.Div(
                                    id="data_card",
                                    children=[
                                        html.H6("Raw data"),
                                        html.Br(),
                                        html.Div(id = 'div-selected-genes', children = []),
                                    ],
                                ),
                            ],
                        ), ]),
                    ],
                ),
                dcc.Tab(
                    label='Statistical analysis',
                    value = 'main-tab-2',
                    style=tab_style, selected_style=tab_selected_style,
                    children=[
                        # Left column
                        html.Div(
                            id="left-column2",
                            className="three columns",
                            children=[description_card2(), generate_control_card2()]
                        ),
                        dcc.Loading(
                            id="loading2",
                            color = colorscale[3],
                            type="circle",
                            children=[


                        # Right column
                        html.Div(
                            id="right-column2",
                            className="nine columns",
                            children=[
                                html.Div(
                                    id = 'analysis-id',
                                    children=[
                                    generate_analysis_card1(),
                                    generate_analysis_card2(),
                                    generate_analysis_card3(),
                                    generate_analysis_card4()
                                    ]
                                ),
                            ],
                        ), ]),
                    ],
                ),
            ],
        ),
    html.Div(children = "©2020 Naef lab",
            style = {'position':'fixed',
                     'bottom':'0',
                     'right':'0',
                     'left':'0',
                     'background':colorscale[3],
                     'padding':'10px',
                     'box-sizing':'border-box',
                     'color':'white',
                     }
            )
    ],
)

""" CALLBACK FOR TAB 1 """
@app.callback(
    [Output('graph-space', 'figure'),
     Output('graph-time', 'figure'),
     Output('graph-3d', 'figure'),
     Output('div-selected-genes', 'children'),
     ],
    [Input('gene-select', 'value')])#, Input('yaxis-type', 'value'), Input('yaxis-scale', 'value'), Input('data-type', 'value')])
def update_figure_time(gene_name):#, yaxis_type, yaxis_scale, data_type, value_tab ):
    if gene_name is None:
        raise Exception()
    else:

        fig_space = Figures.compute_figure_space(gene_data[gene_name])
        fig_time = Figures.compute_figure_time(gene_data[gene_name])
        fig_3D = Figures.compute_figure_3D(gene_data[gene_name])

        data = gene_data[gene_name]

        array_gene = data['rep1']
        array_gene_std = data['rep1_std']
        array_gene_2 = data['rep2']
        array_gene_std_2 = data['rep2_std']
        array_gene_3 = data['rep3']
        array_gene_std_3 = data['rep3_std']

        l_tr = [ [ Styled_th('      ', {'background-color': colorscale[5]}) ] + [ Styled_th('x = ' + str(x)) for x in range(8) ] ]
        for idx, i in enumerate(range(0,24,6)):
            l_th  = [ Styled_th('t = ' + str(i) + 'h', {'background-color': colorscale[5]}) ]
            for j in range(0,8,1):
                if i==0 or i==12:
                    l_th.append( Styled_th( format(round(array_gene[j][idx],2)) + ', ' + format(round(array_gene_2[j][idx],2)) + ', ' + format(round(array_gene_3[j][idx],2))  , small = True ) )
                else:
                    l_th.append( Styled_th( format(round(array_gene[j][idx],2)) + ', ' + format(round(array_gene_2[j][idx],2)) , small = True )  )
            l_tr.append(l_th)
        table = html.Table( [html.Tr(l, style = { 'background-color': colorscale[5]}) if m==0 else html.Tr(l) for m,l in enumerate(l_tr)], style = {'border-collapse': 'collapse' , 'width': '100%'})

        return fig_space, fig_time, fig_3D, table



@app.callback(
     Output('paragraph-selected-gene-2', 'children'),
    [Input('gene-select', 'value')])
def update_figure_time(gene_name):
    if gene_name is None:
        raise Exception()
    else:
        return dcc.Markdown('**The gene you selected is: ' + gene_name.capitalize() + '**')

""" CALLBACK FOR TAB 2 """
@app.callback([
     Output('analysis-card1', 'style'),
     Output('analysis-card2', 'style'),
     Output('analysis-card3', 'style'),
     Output('analysis-card4', 'style'),
     ],
     [Input('analysis-select', 'value')])
def update_figure_time(value):
    if value =='spatial':
        return {}, {'display' : 'none'}, {'display' : 'none'}, {'display' : 'none'}
    elif value == 'temporal':
        return {'display' : 'none'}, {}, {'display' : 'none'}, {'display' : 'none'}
    elif value == 'validation':
        return {'display' : 'none'}, {'display' : 'none'}, {'display' : 'none'}, {}
    else:
        return {'display' : 'none'}, {'display' : 'none'}, {}, {'display' : 'none'}


@app.callback(
    Output('summary-space-regression', 'children'),
    [Input('tabs-time', 'value'), Input('tabs', 'value'), Input('analysis-select', 'value')],
    [State('gene-select', 'value')])#, Input('yaxis-type', 'value'), Input('yaxis-scale', 'value'),Input('tabs', 'value')])
def update_figure_fits(value, value_main_tab, value_analysis, gene_name): #yaxis_type, yaxis_scale, value_tab ):
    if gene_name is None:
        raise Exception()
    else:
        if value_main_tab == 'main-tab-2' and value_analysis == 'spatial':
            #correct value
            value = int(value)-1
            t = int(value)


            array_gene_space = np.concatenate( (gene_data[gene_name]['rep1'], gene_data[gene_name]['rep2'], gene_data[gene_name]['rep3']), axis = 0)

            if t==0 or t==2:
                selected, B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred = LinearRegression.make_space_regression(array_gene_space[:,t], predict= True)
            else:
                selected, B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred = LinearRegression.make_space_regression(array_gene_space[:16,t], predict= True)

            if len(selected) == 1:
                str_param = dash_katex.DashKatex(id='katex_a0', expression=' \mu_0')
            else:
                str_param = dash_katex.DashKatex(id='katex_other_parameters', expression = return_str_list_param(selected) )


            space_domain = np.concatenate((np.linspace(0,8,8, endpoint = False),np.linspace(0,8,8, endpoint = False), np.linspace(0,8,8, endpoint = False)))

            if t==0 or t==2:
                figure =  Figures.compute_figure_space_tab_3(space_domain, array_gene_space[:,t], X_pred, Y_pred)#, yaxis_type, yaxis_scale)
            else:
                figure = Figures.compute_figure_space_tab_3(space_domain[:16], array_gene_space[:16,t], X_pred, Y_pred)#, yaxis_type, yaxis_scale)

            return [html.Div(children = [
                                        html.P('Retained parameters: '),
                                        html.Div(style = {'width' : '5px'}),str_param]
                                        ,
                            style = {'display' : 'flex', 'justify-content':'center'}
                            ),
                    dcc.Graph(id='graph-stat-space',
                              figure = figure,
                              config={'displayModeBar': False, 'modeBarButtonsToRemove': [], 'displaylogo' : False },
                              style = {'width' : '60vw'}
                            )
                    ]
        else:
            raise PreventUpdate


@app.callback(
    [Output('graph-polar', 'figure'),Output('graph-mean', 'figure')],
    [Input('tabs', 'value'), Input('analysis-select', 'value')],
    [State('gene-select', 'value')])
def update_figure_polar(value_main_tab, value_analysis, gene_name):
    if gene_name is None:
        raise Exception()
    else:
        if value_main_tab == 'main-tab-2' and value_analysis == 'temporal':
            array_gene_time =np.concatenate( (gene_data[gene_name]['rep1'], gene_data[gene_name]['rep2'], gene_data[gene_name]['rep3'][:,[0,2]]), axis = 1)

            l_time_reg = []
            for x in range(8):
                l_time_reg.append(LinearRegression.make_time_regression(array_gene_time[x,:], simple = False, predict= True))
            l_time_reg_simple = []
            for x in range(8):
                l_time_reg_simple.append(LinearRegression.make_time_regression(array_gene_time[x,:], simple = True, predict= False))

            figure_polar = Figures.compute_figure_polar_tab_3(l_time_reg)
            figure_mean = Figures.compute_figure_mean_tab_3(l_time_reg)#, yaxis_type, yaxis_scale)

            return figure_polar, figure_mean
        else:
            raise PreventUpdate

@app.callback(
    Output('summary-time-regression', 'children'),
    [Input('tabs-space', 'value'), Input('tabs', 'value'), Input('analysis-select', 'value')],
    [State('gene-select', 'value')])
def update_figure_polar(value, value_main_tab, value_analysis, gene_name):
    if gene_name is None:
        raise Exception()
    else:
        if value_main_tab == 'main-tab-2' and value_analysis == 'temporal':
            array_gene_time =np.concatenate( (gene_data[gene_name]['rep1'], gene_data[gene_name]['rep2'], gene_data[gene_name]['rep3'][:,[0,2]]), axis = 1)

            l_time_reg = []
            for x in range(8):
                l_time_reg.append(LinearRegression.make_time_regression(array_gene_time[x,:], simple = False, predict= True))
            l_time_reg_simple = []
            for x in range(8):
                l_time_reg_simple.append(LinearRegression.make_time_regression(array_gene_time[x,:], simple = True, predict= False))

            #correct value
            value = int(value)-1

            B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred = l_time_reg[value]
            [mu_1, a_1, b_1] = B.flatten()
            [std_mu_1, std_a_1, std_b_1] = np.diagonal(SE)
            bic_1 = bic
            aic_1 = aic
            r2_1 = adj_r2

            B_simple, SE_simple, adj_r2_simple, aic_simple, bic_simple, pv_simple = l_time_reg_simple[value]
            [mu_2] =  B_simple.flatten()
            [std_mu_2] = np.diagonal(SE_simple)
            bic_2 = bic_simple
            aic_2 = aic_simple
            r2_2 = adj_r2_simple

            table_model_1 = html.Table([html.Tr([Styled_th('Parameter'), Styled_th('Mean'), Styled_th('SE')], style = { 'background-color': colorscale[5]}),
                                        html.Tr([Styled_th(dash_katex.DashKatex(expression='\mu'), { 'background-color': colorscale[5]}), Styled_th('{:.2e}'.format(mu_1)), Styled_th('{:.2e}'.format(std_mu_1))]),
                                        html.Tr([Styled_th(dash_katex.DashKatex(expression='a'), { 'background-color': colorscale[5]}), Styled_th('{:.2e}'.format(a_1)), Styled_th('{:.2e}'.format(std_a_1))]),
                                        html.Tr([Styled_th(dash_katex.DashKatex(expression='b'), { 'background-color': colorscale[5]}), Styled_th('{:.2e}'.format(b_1)), Styled_th('{:.2e}'.format(std_b_1))])
                                        ])

            table_model_2 = html.Table([html.Tr([Styled_th('Parameter'), Styled_th('Mean'), Styled_th('SE')], style = { 'background-color': colorscale[5]}),
                                        html.Tr([Styled_th(dash_katex.DashKatex(expression='\mu'), { 'background-color': colorscale[5]}), Styled_th('{:.2e}'.format(mu_2)), Styled_th('{:.2e}'.format(std_mu_2))])
                                        ])
            table_comparison = html.Table([html.Tr([Styled_th('Model'), Styled_th('BIC'), Styled_th('AIC'), Styled_th(dash_katex.DashKatex(expression='\\text{Adj. } R^2') )], style = { 'background-color': colorscale[5]}),
                                        html.Tr([Styled_th('Intercept-only', { 'background-color': colorscale[5]}), Styled_th('{:.2e}'.format(bic_2)), Styled_th('{:.2e}'.format(aic_2)), Styled_th('{:.2e}'.format(r2_2))]),
                                        html.Tr([Styled_th('Oscillatory', { 'background-color': colorscale[5]}), Styled_th('{:.2e}'.format(bic_1)), Styled_th('{:.2e}'.format(aic_1)), Styled_th('{:.2e}'.format(r2_1))])
                                        ])


            time_domain = np.concatenate((np.linspace(0,24,4, endpoint = False),np.linspace(0,24,4, endpoint = False), np.linspace(0,24,2, endpoint = False)))
            x = int(value)
            B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred = l_time_reg[x]
            figure = Figures.compute_figure_time_tab_3(time_domain, array_gene_time[x,:], X_pred, Y_pred)#, yaxis_type, yaxis_scale)

            return [html.Div(children = [
                        html.Div(children = [html.H6('Intercept-only model', style = {'textAlign' : 'center'}) , table_model_2]),
                        html.Div(children = [html.H6('Oscillatory model', style = {'textAlign' : 'center'}), table_model_1]),
                        html.Div(children = [html.H6('Models comparison', style = {'textAlign' : 'center'}), table_comparison, html.P('P-value associated with the oscillatory model (ANOVA): ' + str(pv))], style = {'display' : 'flex', 'flex-direction': 'column'}),
                    ],
                    style = {'display' : 'flex', 'flex-direction': 'row', 'justify-content' : 'space-around', 'flex-wrap' : 'wrap', 'flex-align':'baseline'}),
                    dcc.Graph(id='graph-stat-time', figure = figure, config={'displayModeBar': False, 'modeBarButtonsToRemove': [], 'displaylogo' : False }, style = {'width' : '60vw'})
                    ]
        else:
            raise PreventUpdate


@app.callback(
    Output('summary-2D-regression', 'children'),
    [Input('tabs', 'value'), Input('analysis-select', 'value')],
    [State('gene-select', 'value')])#, Input('yaxis-type', 'value'), Input('yaxis-scale', 'value'),Input('tabs', 'value')])
def update_figure_fits(value_main_tab, value_analysis, gene_name):#, yaxis_type, yaxis_scale, value_tab ):
    if gene_name is None:
        raise Exception()
    else:
        if value_main_tab == 'main-tab-2' and value_analysis == 'spatiotemporal':
            array_gene_time = np.concatenate( (gene_data[gene_name]['rep1'], gene_data[gene_name]['rep2'], gene_data[gene_name]['rep3'][:,[0,2]]), axis = 1)
            fig_3D = Figures.compute_figure_3D_tab_3(dic_reg[gene_name], array_gene_time)#, yaxis_type, yaxis_scale)

            selected = dic_reg[gene_name][0]
            pv = dic_reg[gene_name][6]
            set_selected = set(selected)
            if len(selected) == 1:
                str_param = dash_katex.DashKatex(expression='\\text{ } \mu_0\\text{. }')
                str_param2 = 'This corresponds to the flat model.'
            else:
                str_param = dash_katex.DashKatex(expression='\\text{ }'+return_str_list_param(selected)+'\\text{. }')

                if set_selected == set(['mu0', 'a0','b0']) or set_selected ==set(['mu0', 'a0']) or set_selected == set(['mu0', 'b0']):
                    str_param2= "This corresponds to the rhythmic model."

                elif 'a0' not in selected and 'a1' not in selected and 'a2' not in selected and 'b0' not in selected and 'b1' not in selected and 'b2' not in selected:
                    if 'mu1' in selected or 'mu2' in selected:
                        str_param2 = "This corresponds to the zonated model."
                else:
                    str_param2 = "This corresponds to the rhythmic-zonated model."



            l_time_reg = []
            for x in range(8):
                l_time_reg.append(LinearRegression.make_time_regression(array_gene_time[x,:], simple = False, predict= True))

            fig_polar = Figures.compute_figure_polar_fit_tab_3(dic_reg[gene_name],l_time_reg)
            if fig_polar!=None:
                return [html.Div(children=[html.P("Retained parameters: "), str_param, html.P(str_param2)], style = {'display' : 'flex', 'justify-content':'center'}),
                        html.Div(children = [
                                            dcc.Graph(id='graph-polar-fit',
                                                      figure =  fig_polar,
                                                      config={'displayModeBar': False, 'modeBarButtonsToRemove': [], 'displaylogo' : False },
                                                      style = {'width' : '25vw'}
                                                      ),
                                            dcc.Graph(id='graph-fit-3D',
                                                      figure = fig_3D,
                                                      config={'modeBarButtonsToRemove': ['sendDataToCloud', 'resetCameraLastSave3d', 'hoverClosest3d', 'zoom3d', 'toImage'], 'displaylogo' : False },
                                                      style = {'border-width':'1px', 'border-style':'solid', 'border-color':'#e8e8e8', 'width' : '45vw'} )
                                            ],
                                style = {'display' : 'flex', 'flex-direction': 'row', 'justify-content' : 'space-between', 'height': '500px'},
                                ),
                        ]
            else:
                return [html.Div(children=[html.P("Retained parameters: "), str_param, html.P(str_param2)], style = {'display' : 'flex', 'justify-content':'center'}),
                        html.Div(children = [
                                            dcc.Graph(id='graph-fit-3D', figure = fig_3D, config={'modeBarButtonsToRemove': ['sendDataToCloud', 'resetCameraLastSave3d', 'hoverClosest3d', 'zoom3d', 'toImage'], 'displaylogo' : False },
                                                    style = {'border-width':'1px', 'border-style':'solid', 'border-color':'#e8e8e8', 'width' : '60vw'} )
                                            ],
                                style = {'display' : 'flex', 'flex-direction': 'row', 'justify-content' : 'space-around'}
                                ),
                        ]
        else:
            raise PreventUpdate

@app.callback(
    Output('graph-comparison', 'figure'),
    [Input('tabs', 'value'), Input('analysis-select', 'value')],
    [State('gene-select', 'value')])
def make_graph(value_main_tab, value_analysis, gene_name):#, yaxis_type, yaxis_scale, data_type, value_tab):
    if gene_name is None or gene_name not in gene_data_atg:
        raise Exception()
    else:
        if value_main_tab == 'main-tab-2' and value_analysis == 'validation':

            data_atg = gene_data_atg[gene_name]
            data_itz =  gene_data[gene_name]

            array_atg = np.nanmean( data_atg, axis = 1)
            array_itz =  np.nanmean(np.nanmean( [data_itz['rep1'], data_itz['rep2'], data_itz['rep3']], axis = 0), axis = 0)

            return Figures.compute_figure_comparison(array_atg, array_itz)
        else:
            raise PreventUpdate

@app.callback(
    Output('graph-polar-comparison', 'figure'),
    [Input('tabs', 'value'), Input('analysis-select', 'value')],
    [State('gene-select', 'value')])
def make_graph(value_main_tab, value_analysis, gene_name):#, yaxis_type, yaxis_scale, data_type, value_tab):
    if gene_name is None or gene_name not in gene_data_atg:
        raise Exception()
    else:
        if value_main_tab == 'main-tab-2' and value_analysis == 'validation':
            data_atg = gene_data_atg[gene_name]
            data_itz =  gene_data[gene_name]

            array_atg = np.nanmean( data_atg, axis = 1)
            array_itz =  np.nanmean(np.nanmean( [data_itz['rep1'], data_itz['rep2'], data_itz['rep3']], axis = 0), axis = 0)

            return Figures.compute_figure_polar_comparison(array_atg, array_itz)
        else:
            raise PreventUpdate

@app.callback(
    [Output('graph-comparison', 'style'),
     Output('graph-polar-comparison', 'style')
     ],
    [Input('tabs', 'value'), Input('analysis-select', 'value')],
    [State('gene-select', 'value')])
def make_graph(value_main_tab, value_analysis, gene_name):#, yaxis_type, yaxis_scale, data_type, value_tab):
    if value_main_tab == 'main-tab-2' and value_analysis == 'validation':
        if gene_name is None or gene_name not in gene_data_atg:
            return {'display' : 'none'}, {'display' : 'none'}
        else:
            return {'width' : '40vw'}, {'width' : '25vw'}
    else:
        raise PreventUpdate

@app.callback(
    Output('no-gene', 'children'),
    [Input('tabs', 'value'), Input('analysis-select', 'value')],
    [State('gene-select', 'value')])
def make_graph(value_main_tab, value_analysis, gene_name):#, yaxis_type, yaxis_scale, data_type, value_tab):
    if value_main_tab == 'main-tab-2' and value_analysis == 'validation':
        if gene_name is None or gene_name not in gene_data_atg:
            return dcc.Markdown('This gene is not available in *Atger & al dataset*')
        else:
            return ''
    else:
        raise PreventUpdate


""" IMPORTANT FUNCTIONS """
def Styled_th(x, supp = {}, small = False):
    if small:
        style = { 'border': '1px solid #dddddd', 'text-align': 'center', 'padding': '4px', 'font-size' : 'small', 'font-weight' : 'lighter'}
    else:
        style = { 'border': '1px solid #dddddd', 'text-align': 'center', 'padding': '4px', 'font-weight' : 'lighter'}
    for key, val in supp.items():
        style[key] = val
    return  html.Th( x , style = style)

def return_str_list_param(l_param, b_sorted = False):
    str_p = ''
    if b_sorted:
        l_param = sorted(l_param)
    for param in l_param:
        if len(param)>3:
            p1, p2 = param.split('+')
            if len(p1)>=3:
                p1 = '\\' + p1
            str_p += p1[:-1]+'_'+p1[-1] + ', ' + p2[:-1]+ '_'+p2[1]+ ', '
        else:
            if len(param)>=3:
                param = '\\' + param
            str_p += param[:-1] + '_' +  param[-1]  + ', '
    return str_p[:-2]


# Run the server
if __name__ == "__main__":
    app.run_server(debug=False)
