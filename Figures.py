import numpy as np
import plotly.graph_objs as go
import LinearRegression
import operator
import seaborn as sn

colorscale = [ 'rgb(165,0,38)', 'rgb(215,48,39)', 'rgb(244,109,67)',
             'rgb(253,174,97)', 'rgb(254,224,144)', 'rgb(224,243,248)',
             'rgb(171,217,233)', 'rgb(116,173,209)']

colorscale_1 = sn.color_palette("husl", 4).as_hex()
colorscale_2 = sn.color_palette("GnBu_d",8).as_hex()

def compute_figure_time(data_gene):#, yaxis_type, yaxis_scale, data_type):
    data_type = "Log"
    yaxis_scale = "rel"
    yaxis_type = "Linear"

    array_gene_1 = data_gene['rep1']
    array_gene_2 = data_gene['rep2']
    array_gene_3 = data_gene['rep3']

    if data_type == 'Linear':
        array_gene_1 = np.exp(array_gene_1)
        array_gene_2 = np.exp(array_gene_2)
        array_gene_3 = np.exp(array_gene_3)

    traces = []
    for i in range(8):
            traces.append(go.Scatter(x = list(range(0,24,6)),
                                     y = array_gene_1[i,:],
                                     text = "x = " +str(i),
                                     mode = 'markers',
                                     marker = dict( color = colorscale_2[i]),
                                     showlegend = False
                                     )
                         )

            traces.append(go.Scatter(x = list(range(0,24,6)),
                                     y = array_gene_2[i,:],
                                     text = "x = " +str(i),
                                     mode = 'markers',
                                     marker = dict( color = colorscale_2[i]),
                                     showlegend = False
                                     )
                        )
            avg = [0]*4
            for j in range(4):
                if j==0 or j==2:
                    avg[j]=1/3 * array_gene_1[i,j]+1/3 * array_gene_2[i,j]+1/3 * array_gene_3[i,j]
                else:
                    avg[j] = 0.5 * array_gene_1[i,j]+0.5 * array_gene_2[i,j]

            traces.append(go.Scatter(x = list(range(0,24,6)),
                                     y = avg ,
                                     text = "x = " +str(i),
                                     name = "x = " +str(i),
                                     mode = 'lines',
                                     line = dict( color = colorscale_2[i])
                                     )
                        )

    if yaxis_scale=='abs':
        if yaxis_type == 'Log':
            yrange = [min(0.,np.log10(np.amin(array_gene_1)*1.1), np.log10(np.amin(array_gene_2)*1.1)),
                      max(np.log10(np.amax(array_gene_1)*1.1),np.log10(np.amax(array_gene_2)*1.1))
                      ]
        else:
            yrange = [0, max(np.amax(array_gene_1),np.amax(array_gene_2)) *1.1]
    else:
        yrange = None

    figure = go.Figure(data=traces,
                      layout = go.Layout(xaxis={'type': 'linear', 'title': 'Time (h)'},
                                         yaxis={'title': 'Expression [Log2]',
                                                'type': 'linear' if yaxis_type == 'Linear' else 'log',
                                                'range' : yrange },
                                         hovermode='closest',
                                         title=dict(x=0.5, text = 'Temporal profile'),
                                         template='plotly_white'
                                        )
                        )
    return figure


def compute_figure_space(data_gene):#, yaxis_type, yaxis_scale, data_type):

    data_type = "Log"
    yaxis_scale = "rel"
    yaxis_type = "Linear"

    array_gene_1 = data_gene['rep1']
    array_gene_2 = data_gene['rep2']
    array_gene_3 = data_gene['rep3']

    if data_type == 'Linear':
        array_gene_1 = np.exp(array_gene_1)
        array_gene_2 = np.exp(array_gene_2)
        array_gene_3 = np.exp(array_gene_3)

    traces = []
    for i in range(4):
        traces.append(go.Scatter(x = list(range(8)),
                                 y = array_gene_1[:,i],
                                 text = "t = " +str(i*6) + 'h', mode = 'markers',
                                 marker = dict( color = colorscale_1[i]),  showlegend = False
                                 )
                    )
        traces.append(go.Scatter(x = list(range(8)),
                                 y = array_gene_2[:,i],
                                 text = "t = " +str(i*6) + 'h',
                                 mode = 'markers',
                                 marker = dict( color = colorscale_1[i]),
                                 showlegend = False
                                 )
                    )
        traces.append(go.Scatter(x = list(range(8)),
                                 y = array_gene_3[:,i],
                                 text = "t = " +str(i*6) + 'h',
                                 mode = 'markers',
                                 marker = dict( color = colorscale_1[i]),
                                 showlegend = False
                                 )
                     )

        if i==0 or i==2:
            avg = 1/3*array_gene_1[:,i]+1/3*array_gene_2[:,i]+1/3*array_gene_3[:,i]
        else:
            avg = 0.5*array_gene_1[:,i]+0.5*array_gene_2[:,i]
        traces.append(go.Scatter(x = list(range(8)),
                                 y = avg,
                                 text = "t = " +str(i*6) + 'h',
                                 name = "t = " +str(i*6) + 'h',
                                 mode = 'lines',
                                 marker = dict( color = colorscale_1[i])
                                 )
                     )

    if yaxis_scale=='abs':
        if yaxis_type == 'Log':
            yrange = [min(0.,np.log10(np.amin(array_gene_1)*1.1), np.log10(np.amin(array_gene_2)*1.1) ),
                      max(np.log10(np.amax(array_gene_1)*1.1),np.log10(np.amax(array_gene_2)*1.1))
                     ]
        else:
            yrange = [0, max(np.amax(array_gene_1),np.amax(array_gene_2)) *1.1]
    else:
        yrange = None

    figure = go.Figure(data = traces,
                       layout = go.Layout(
                                        xaxis={'type': 'linear', 'title': 'Spatial position (layer)'},
                                        yaxis={'title': 'Expression [Log2]',
                                               'type': 'linear' if yaxis_type == 'Linear' else 'log',
                                               'range' : yrange },
                                        #margin={'l': 100, 'b': 40, 't': 0, 'r': 40},
                                        hovermode='closest',
                                        #autosize = True,
                                        title=dict(x=0.5, text = 'Spatial profile'),
                                        template='plotly_white'
                                        )
                        )
    return figure

def compute_figure_3D(data_gene):#, yaxis_type, yaxis_scale, data_type):

    data_type = "Log"
    yaxis_scale = "rel"
    yaxis_type = "Linear"

    array_gene_1 = data_gene['rep1']
    array_gene_2 = data_gene['rep2']
    array_gene_3 = data_gene['rep3']

    if data_type == 'Linear':
        array_gene_1 = np.exp(array_gene_1)
        array_gene_2 = np.exp(array_gene_2)
        array_gene_3 = np.exp(array_gene_2)

    #surface parametrization
    x_surf = np.array(list(range(8)))
    y_surf = np.array(list(range(0,24,6)))
    z_surf = np.zeros((8,4))
    for i in range(8):
        for j in range(4):
            if j==0 or j==2:
                z_surf[i,j] = 1/3*array_gene_1[i,j]+1/3*array_gene_2[i,j]+1/3*array_gene_3[i,j]
            else:
                z_surf[i,j] = 0.5*array_gene_1[i,j]+0.5*array_gene_2[i,j]

    #points parametrization
    x = [i for i in range(8) for j in range(4)]
    y = [j for i in range(8) for j in range(0,24,6)]
    z_1 = [array_gene_1[i][j] for i in range(8) for j in range(4)]
    z_2 = [array_gene_2[i][j] for i in range(8) for j in range(4)]
    z_3 = [array_gene_3[i][j] for i in range(8) for j in range(4)]

    if yaxis_scale=='abs':
        if yaxis_type == 'Log':
            zrange = [min(0.,np.log10(np.amin(array_gene_1)*1.1), np.log10(np.amin(array_gene_2)*1.1) ),
                      max(np.log10(np.amax(array_gene_1)*1.1),np.log10(np.amax(array_gene_2)*1.1))
                     ]
        else:
            zrange = [0, max(np.amax(array_gene_1),np.amax(array_gene_2)) *1.1]
    else:
        zrange = None

    points_1 = go.Scatter3d( x = y+y+y,
                             y = x+x+x,
                             z = z_1+z_2+z_3,
                             mode = 'markers',
                             marker=dict(color = colorscale_2[1],
                                         size=4,
                                         symbol='circle',
                                         line=dict( color = colorscale_2[3], width=1 ),
                                         opacity=0.8
                                         )
                            )
    fit = go.Surface( x = y_surf,
                      y = x_surf,
                      z = z_surf ,
                      name = "Linear fit",
                      opacity = 0.7,
                      showscale = False,
                      hoverinfo = 'none',
                      colorscale="YlGnBu"
                      )



    zoom_factor = 6.5
    camera = dict( up=dict(x=0, y=0, z=1),
                   center=dict(x=0, y=0, z=0),
                   eye=dict(x=0.25*zoom_factor,
                            y=0.01*zoom_factor,
                            z=0.02*zoom_factor
                            )
                 )
    layout = go.Layout( scene = dict(
                                    xaxis = dict(title='Time (h)'),
                                    yaxis = dict(title='Spatial position (layer)'),
                                    zaxis = dict(title='Expression [Log2]',
                                                 type = 'linear' if yaxis_type == 'Linear' else 'log',
                                                 range = zrange ),
                                    camera = camera,
                                    aspectmode = 'manual',
                                    aspectratio =dict(x = 1, y = 2, z = 0.5)
                                    ),
                        title=dict(x=0.5, text = '3D profile'),
                        margin=dict(r=0, l=0, b=0, t=70),
                        showlegend=False,
                        )

    figure = go.Figure(data=[ fit, points_1], layout=layout)

    return figure

def compute_figure_space_tab_3(X, Y, X_pred, Y_pred):#, yaxis_type, yaxis_scale):

    #if yaxis_scale=='abs':
    #    if yaxis_type == 'Log':
    #        yrange = [min(0.,np.log10(np.amin(Y)*1.1) ), np.log10(np.amax(Y)*1.1)]
    #    else:
    #        yrange = [0, np.amax(Y)*1.1]
    #else:
    #    yrange = None

    data_type = "Log"
    yaxis_scale = "rel"
    yaxis_type = "Linear"
    yrange = None

    #plot gene vs prediction
    traces = [ go.Scatter(x = X, y = Y,  name = "Experimental signal", mode = 'markers', marker = dict( color = colorscale_2[3])),
               go.Scatter(x = X_pred, y = Y_pred ,  name = "Polynomial regression", mode = 'lines', marker = dict( color = colorscale_2[3]))]


    layout = go.Layout( xaxis={'type': 'linear', 'title': 'Layer'},
                        yaxis={'title': 'Expression [Log2]',
                               'type': 'linear' if yaxis_type == 'Linear' else 'log',
                               'range' : yrange
                               },
                        margin=dict(r=0, l=50, b=100, t=50),
                        #title = 'Spatial regression',
                        template='plotly_white',
                        )

    figure = go.Figure(data=traces, layout=layout)

    return figure

def compute_figure_time_tab_3(X, Y, X_pred, Y_pred):#, yaxis_type, yaxis_scale):
    #if yaxis_scale=='abs':
    #    if yaxis_type == 'Log':
    #        yrange = [min(0.,np.log10(np.amin(Y)*1.1) ), np.log10(np.amax(Y)*1.1)]
    #    else:
    #        yrange = [0, np.amax(Y)*1.1]
    #else:
    #    yrange = None

    data_type = "Log"
    yaxis_scale = "rel"
    yaxis_type = "Linear"
    yrange = None


    #plot gene vs prediction
    traces = [ go.Scatter(x = X,
                          y = Y,
                          name = "Experimental signal",
                          mode = 'markers',
                          marker = dict( color = colorscale_2[3]),
                          ),
               go.Scatter(x = X_pred,
                          y = Y_pred ,
                          name = "Prediction with 3 parameters",
                          mode = 'lines',
                          marker = dict( color = colorscale_2[1])
                          ) ,
               go.Scatter(x = X_pred,
                          y = [np.mean(Y)]*len(X_pred) ,
                          name = "Prediction with 1 parameter",
                          mode = 'lines',
                          marker = dict( color = colorscale_2[5])
                          )
                ]

    layout = go.Layout( xaxis={'type': 'linear', 'title': 'Time (h)'},
                        yaxis={'title': 'Expression [Log2]',
                               'type': 'linear' if yaxis_type == 'Linear' else 'log',
                               'range' : yrange},
                        margin=dict(r=0, l=50, b=100, t=50),
                        #title='Regression: mean',
                        template='plotly_white'
                        )
    figure = go.Figure(data=traces, layout=layout)

    return figure


def compute_figure_polar_tab_3(l_full_reg):
    max_r = 0
    min_r = 10000

    traces  = []
    for x in range(8):
        B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred = l_full_reg[x]
        [mu, a, b] = B.flatten()
        [std_mu, std_a, std_b] = np.diagonal(SE)
        l_lambda, P = np.linalg.eig(np.nan_to_num(SE[1:,1:]))
        angle = np.arctan2(b,a)/(2*np.pi)*360
        r = np.sqrt(a**2+b**2)
        if r>max_r:
            max_r = r
        if r<min_r:
            min_r = r
        ellipse_domain = np.linspace(0,2*np.pi,80, endpoint = True)
        ellipse_x = [a + l_lambda[0] * np.cos(d) for d in ellipse_domain]
        ellipse_y = [b + l_lambda[1] * np.sin(d) for d in ellipse_domain]
        ellipse_project = np.matmul(np.vstack((ellipse_x, ellipse_y)).T,P)
        ellipse_r = [np.sqrt(x**2+y**2) for x,y in zip(ellipse_project[:,0], ellipse_project[:,1])]
        ellipse_t = [np.arctan2(y,x)/(2*np.pi)*360 for x,y in zip(ellipse_project[:,0], ellipse_project[:,1])]
        traces.append( [ r, angle, a, b, std_a, std_b, ellipse_r, ellipse_t, x   ]      )
        #print(r, angle, a, b, std_a, std_b, ellipse_r, ellipse_t, x   ] )
    traces = sorted(traces, key=operator.itemgetter(4), reverse = True)


    data_traces = []
    for trace in traces:
        data_traces += [ go.Scatterpolar( r = [trace[0]],
                                          theta = [trace[1]],
                                          thetaunit = "degrees",
                                          mode='markers',
                                          name='x = ' + str(trace[8]),
                                          marker=dict( size=20, color=colorscale_2[trace[8]], opacity=0.9),
                                          showlegend=False
                                         )
                        ]
        data_traces += [ go.Scatterpolar( r = trace[6],
                                          theta = trace[7],
                                          thetaunit = "degrees",
                                          mode='lines',
                                          name = 'x_std = ' + str(trace[8]),
                                          line=dict( width=10, color=colorscale_2[trace[8]]),
                                          opacity = 0.7,
                                          showlegend=False
                                          )
                        ]

    #create legend separately to have the good order
    for x in range(8):
        data_traces += [  go.Scatterpolar( r = [1000],
                                           theta = [0],
                                           mode='markers',
                                           name='x = ' + str(x),
                                           marker=dict( size=20, color=colorscale_2[x], opacity=0.9),
                                           showlegend=True,
                                           visible = True
                                           )
                        ]

    layout = go.Layout( title=dict(x=0.5, text = 'Regression: amplitude and phase'),
                        plot_bgcolor='rgb(223, 223, 223)',
                        polar = dict( angularaxis=dict( tickcolor='rgb(253,253,253)',
                                                        direction  = "clockwise",
                                                        thetaunit = "degrees",
                                                        #categoryarray = [str(x) for x in np.linspace(0,24,4, endpoint = False)]
                                                        ),
                                      radialaxis=dict( range =[ 0,  max_r*1.5 ]),
                                     ),
                        showlegend = True,
                        autosize = True,
                        margin=dict(r=0, l=100, b=50, t=100),
                        template='plotly_white'
                        )

    layout.polar.angularaxis.tickvals = [hr2angle(hr) for hr in range(0,24,3)]
    layout.polar.angularaxis.ticktext = [hr_str(hr) for hr in range(0,24,3)]

    figure = go.Figure(data=data_traces, layout=layout)

    return figure



def compute_figure_mean_tab_3(l_full_reg):#, yaxis_type, yaxis_scale):
    l_mu = []
    l_std_mu = []
    for x in range(8):
        B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred = l_full_reg[x]
        [mu, a, b] = B.flatten()
        [std_mu, std_a, std_b] = np.diagonal(SE)
        l_mu.append(mu)
        l_std_mu.append(std_mu)

    #if yaxis_scale=='abs':
    #    if yaxis_type == 'Log':
    #        yrange = [min(0.,np.log10(np.amin(l_mu)*1.1) ), np.log10(np.amax(l_mu)*1.1)]
    #    else:
    #        yrange = [0, np.amax(l_mu)*1.1]
    #else:
    #    yrange = None
    data_type = "Log"
    yaxis_scale = "rel"
    yaxis_type = "Linear"
    yrange = None

    data = [ go.Scatter(x = list(range(8)),
                        y = l_mu,
                        mode='lines+markers',
                        error_y=dict(type='data',
                        array=l_std_mu,  visible=True  ),
                        marker = dict( color = colorscale_2[3])
                        )
            ]
    layout = go.Layout( title=dict(x=0.5, text = 'Regression: mean  [Log2]'),
                        autosize = True,
                        #margin=dict(r=100, l=200, b=100, t=100),
                        xaxis={'type': 'linear', 'title': 'Spatial position (layer)'},
                        yaxis={'title': 'Expression',
                               'type': 'linear' if yaxis_type == 'Linear' else 'log',
                               'range' : yrange},
                        template='plotly_white'
                        )
    figure = go.Figure(data=data, layout=layout)
    return figure

def compute_figure_3D_tab_3(reg_2D, array_gene_time):#, yaxis_type, yaxis_scale):
    [selected, B, SE, bic, l_schwartz, Xx_pred, Xt_pred, Y_pred, var_exp, var_exp_re] = reg_2D

    #print('other B', B)
    y_formatted = np.zeros((80,40))
    for idx, val in enumerate(Y_pred):
        y_formatted[int(idx/40), int(idx%40)] = val

    #plot gene vs prediction
    x = [i for i in range(8) for j in range(4)]
    y = [j for i in range(8) for j in range(0,24,6)]
    z_1 = [array_gene_time[i][j] for i in range(8) for j in range(4)]
    z_2 = [array_gene_time[i][j] for i in range(8) for j in range(4,8)]
    array_gene_time = np.insert(array_gene_time, 9, np.nan, axis = 1)
    array_gene_time = np.insert(array_gene_time,11, np.nan, axis = 1)
    z_3 = [array_gene_time[i][j] for i in range(8) for j in range(8,12)]


    points_1 = go.Scatter3d( x = y,
                             y = x, z = z_1 ,
                             name = "Experimental points",
                             mode = 'markers',
                             marker = dict(color=colorscale_2[1],
                                         size=4,
                                         symbol='circle',
                                         line=dict( color=colorscale_2[3],  width=1 ),
                                         opacity=0.8
                                         )
                            )
    points_2 = go.Scatter3d( x = y,
                             y = x,
                             z = z_2 ,
                             name = "Experimental points",
                             mode = 'markers',
                             marker = dict(color=colorscale_2[1],
                                           size=4,
                                           symbol='circle',
                                           line=dict( color=colorscale_2[3],  width=1 ),
                                           opacity=0.8
                                           )
                            )
    points_3 = go.Scatter3d( x = y,
                             y = x,
                             z = z_3 ,
                             name = "Experimental points",
                             mode = 'markers',
                             marker = dict(color=colorscale_2[1],
                                         size=4,
                                         symbol='circle',
                                         line=dict( color=colorscale_2[3],  width=1 ),
                                         opacity=0.8
                                         )
                            )
    fit = go.Surface( x = Xt_pred,
                      y = Xx_pred,
                      z = y_formatted ,
                      name = "Linear fit",
                      opacity = 0.7,
                      showscale = False,
                      colorscale = 'YlGnBu',
                      hoverinfo = 'none',)

    #if yaxis_scale=='abs':
    #    if yaxis_type == 'Log':
    #        zrange = [min(0.,np.log10(np.amin(array_gene_time)*1.1) ), np.log10(np.amax(array_gene_time)*1.1)]
    #    else:
    #        zrange = [0, np.amax(array_gene_time)*1.1]
    #else:
    #    zrange = None

    data_type = "Log"
    yaxis_scale = "rel"
    yaxis_type = "Linear"
    yrange = None
    zrange = None

    zoom_factor = 6.5
    camera = dict( up=dict(x=0, y=0, z=1),
                   center=dict(x=0, y=0, z=0),
                   eye=dict(x=0.25*zoom_factor,
                            y=0.01*zoom_factor,
                            z=0.02*zoom_factor
                            )
                   )
    layout = go.Layout( #hovermode='closest' ,
                        #autosize = True,
                        #width=1000,
                        #height = 100,
                        scene = dict(
                        xaxis = dict(title='Time (h)'),
                        yaxis = dict(title='Spatial position (layer)'),
                        zaxis = dict(title='Expression [Log2]',
                                     type = 'linear' if yaxis_type == 'Linear' else 'log',
                                     range = zrange),
                        camera = camera,
                        aspectmode = 'manual',
                        aspectratio =dict(x = 1, y = 2, z = 0.5)
                        ),
                        #title = '<b>Spatiotemporal fit</b>',
                        margin=dict(r=0, l=100, b=0, t=0),
                        #width=80%
                        showlegend=False
                        )

    figure = go.Figure(data=[ fit, points_1, points_2, points_3], layout=layout)

    return figure

def compute_figure_polar_fit_tab_3(reg_2D, l_full_reg):


    #experimental points
    traces  = []
    for x in range(8):
        B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred = l_full_reg[x]
        [mu, a, b] = B.flatten()
        [std_mu, std_a, std_b] = np.diagonal(SE)
        angle = np.arctan2(b,a)/(2*np.pi)*360
        r = np.sqrt(a**2+b**2)
        traces.append( [ r, angle, std_a, x   ]      )
    traces = sorted(traces, key=operator.itemgetter(3), reverse = False)

    data_traces = []
    for trace in traces:
        data_traces += [ go.Scatterpolargl( r = [trace[0]],
                                            theta = [trace[1]],
                                            thetaunit = "degrees",
                                            mode='markers',
                                            name='x = ' + str(trace[3]),
                                            marker=dict( size=20, color=colorscale_2[trace[3]], opacity=0.9),
                                            showlegend=True
                                            )
                        ]

    #fit
    [selected, B, SE, bic, l_schwartz, Xx_pred, Xt_pred, Y_pred, var_exp, var_exp_re] = reg_2D
    if ('a1' in selected and 'b1' in selected) or ('a2' in selected and 'b2' in selected) or ('a3' in selected and 'b3' in selected):
        dic_param = {'a0' : 0, 'a1' : 0, 'a2' : 0, 'b0' : 0, 'b1' : 0, 'b2' : 0, 'mu0' : 0, 'mu1' : 0, 'mu2' : 0}
        param_name = return_str_list_param(selected, False).replace(" ", "").split(',')
        param_value = B
        #print("B", B)
        for par, name_par in zip(param_value, param_name):
            dic_param[name_par] = par

        #print('dic param', dic_param)
        X = np.linspace(0,7,100)
        l_a = dic_param['a0'] + dic_param['a1']*X + dic_param['a2']*0.5*(3*X**2-1)
        l_b = dic_param['b0'] + dic_param['b1']*X + dic_param['b2']*0.5*(3*X**2-1)
        l_r = [(a**2 + b**2)**0.5 for a,b in zip(list(l_a), list(l_b))]

        max_r = np.max(l_r)

        trace = go.Scatterpolargl(r = l_r,
                                  theta = np.arctan2(l_b, l_a)/(2*np.pi)*360,
                                  thetaunit = "degrees",
                                  mode='lines',
                                  name='fit',
                                  line=dict( color='peru', width=5 )
                                  )
        data = [trace]

        layout = go.Layout( #title='Polar fit',
                            font=dict( size=15 ),  plot_bgcolor='rgb(223, 223, 223)',
                            polar = dict( angularaxis = dict( tickcolor='rgb(253,253,253)',
                                                              direction  = "clockwise",
                                                              thetaunit = "degrees",
                                                             ),
                                           #radialaxis=dict( range =[ 0,  max_r*1.5  ] ),
                                       ),
                            showlegend = True,
                            autosize = True,
                            margin=dict(r=40, l=60, b=20, t=50)
                            )

        layout.polar.angularaxis.tickvals = [hr2angle(hr) for hr in range(0,24,3)]
        layout.polar.angularaxis.ticktext = [hr_str(hr) for hr in range(0,24,3)]

        fig_polar = go.Figure(data= data+data_traces, layout=layout)
    else:
        fig_polar = None

    return fig_polar


def compute_figure_comparison(array_atg, array_itz):

    #standardize
    array_itz = (array_itz-np.nanmean(array_itz))/np.nanstd(array_itz)
    array_atg = (array_atg-np.nanmean(array_atg))/np.nanstd(array_atg)

    traces = [go.Scatter(x = list(range(0,24,6)),
                        y = array_itz,
                        text = "Loaded dataset" ,
                        name = "Loaded dataset",
                        mode = 'lines+markers',
                        marker = dict( color = colorscale_2[0]),
                        line = dict( color = colorscale_2[0])
                        ),
                go.Scatter(x = list(range(0,24,2)),
                           y = array_atg,
                           text = "Atger & al. dataset" ,
                           name = "Atger & al. dataset",
                           mode = 'lines+markers',
                           marker = dict( color = colorscale_2[4]),
                           line = dict( color = colorscale_2[4])
                           )
            ]
    data = traces
    layout = go.Layout( xaxis={'type': 'linear', 'title': 'Time (h)'},
                        yaxis={'title': 'Normalized expression',
                               'type': 'linear' ,
                               'zeroline': False},
                        #margin=dict(r=20, l=50, b=20, t=20),
                        title=dict(x=0.5, text = 'Profile comparison'),
                        template='plotly_white'
                       )

    fig = go.Figure(data=data,layout=layout)
    return fig



def compute_figure_polar_comparison(array_atg, array_itz):
    try:
        #standardize
        array_itz = (array_itz-np.nanmean(array_itz))/np.nanstd(array_itz)
        array_atg = (array_atg-np.nanmean(array_atg))/np.nanstd(array_atg)

        max_r = 0
        min_r = 10000
        w = 2*np.pi/24

        ### Regression for itzkovtiz
        B, SE, adj_r2, aic, bic = LinearRegression.make_time_regression_no_replicates(Y = array_itz)
        B = list(B)
        SE = list(np.diagonal(SE))
        mu_itz, a_itz, b_itz = B[0], B[1], B[2]
        std_mu_itz, std_a_itz, std_b_itz = SE[0], SE[1], SE[2]
        angle_itz = np.arctan2(b_itz,a_itz)/(2*np.pi)*360
        r_itz = np.sqrt(a_itz**2+b_itz**2)
        if r_itz>max_r:
            max_r = r_itz
        if r_itz<min_r:
            min_r = r_itz

        ### Regression for atger
        Xt_atg = np.linspace(0,24,12,endpoint = False)
        B, SE, adj_r2, aic, bic = LinearRegression.make_time_regression_no_replicates(Y = array_atg, domain = Xt_atg)
        B = list(B)
        SE = list(np.diagonal(SE))
        #model with cos and sin
        mu_atg, a_atg, b_atg = B[0], B[1], B[2]
        std_mu_atg, std_a_atg, std_b_atg = SE[0], SE[1], SE[2]
        angle_atg = np.arctan2(b_atg,a_atg)/(2*np.pi)*360
        r_atg = np.sqrt(a_atg**2+b_atg**2)
        if r_atg>max_r:
            max_r = r_atg
        if r_atg<min_r:
            min_r = r_atg


        data_traces = [
                        go.Scatterpolargl( r = [1000],
                                           theta = [0],
                                           mode='markers',
                                           name='Loaded dataset',
                                           marker=dict(size=15, opacity=0.7, color = colorscale_2[0]),
                                           showlegend=True,
                                          ),
                        go.Scatterpolargl( r = [1000],
                                           theta = [0],
                                           mode='markers',
                                           name='Atger & al. dataset',
                                           marker=dict( size=15, opacity=0.7, color = colorscale_2[4]),
                                           showlegend=True,
                                          ),

                        go.Scatterpolargl( r = [float(r_itz)],
                                           theta = [float(angle_itz)],
                                           thetaunit = "degrees",
                                           mode ='markers',
                                           name ='Loaded dataset',
                                           marker = dict(size=float(std_a_itz*200/max_r),
                                           opacity=0.7, color = colorscale_2[0]),
                                           showlegend=False
                                           ),
                        go.Scatterpolargl( r = [float(r_atg)],
                                           theta = [float(angle_atg)],
                                           thetaunit = "degrees",
                                           mode='markers',
                                           name='Atger & al. dataset',
                                           marker=dict( size=float(std_a_atg*200/max_r),
                                                        opacity=0.7,
                                                        color = colorscale_2[4]
                                                       ),
                                           showlegend=False
                                        ) ,

                         ]

        layout = go.Layout( polar = dict(radialaxis = dict( visible = True,
                                                            range =[ 0,  float(max_r*2)]
                                                           ),
                                         angularaxis=dict( tickcolor='rgb(253,253,253)',
                                                           direction  = "clockwise",
                                                           thetaunit = "degrees"
                                                           )
                                        ),
                            showlegend = True,
                            autosize = True,
                            #margin=dict(r=0, l=50, b=20, t=50),
                            #margin=dict(r=20, l=20, b=20, t=20),
                            title=dict(x=0.5, text = 'Regressed phase'),
                            template='plotly_white'
                            )

        layout.polar.angularaxis.tickvals = [hr2angle(hr) for hr in range(0,24,3)]
        layout.polar.angularaxis.ticktext = [hr_str(hr) for hr in range(0,24,3)]


        fig = go.Figure(data=data_traces,layout=layout)
        return fig
    except:
        return 'bug'


def return_str_list_param(l_param, b_sorted = False):
    str_p = ''
    if b_sorted:
        l_param = sorted(l_param)
    for param in l_param:
        if len(param)>3:
            p1, p2 = param.split('+')
            str_p += p1 + ', ' + p2+ ', '
        else:
            str_p += param + ', '
    return str_p[:-2]



def angle2hr(angle):
    return (((angle - 15) % 360) + 15) / 15

def hr2angle(hr):
    return (hr * 15) % 360

def hr_str(hr):
    # Normalize hr to be between 1 and 12
    hr_str = str(((hr) % 24) )
    #suffix = ' AM' if (hr % 24) < 12 else ' PM'
    prefix = 'ZT'
    return prefix + hr_str

if __name__ == "__main__":
    pass
