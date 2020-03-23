import numpy as np
from statsmodels.formula.api import ols
import pandas as pd
from statsmodels.stats.anova import anova_lm
import copy
import statsmodels.api as sm
import warnings

np.seterr(divide='raise')
w = 2*np.pi/24

class LinearRegression():
    def __init__(self, Y, dic_X):
        """
        Y: observations vector
        dictionnary of predictors: design matrix
        """
        self.Y = Y[:,np.newaxis]
        self.dic_X = dic_X
        self.n = len(Y)

    def make_complete_regression(self, chosen_predictors, no_intercept = False):
        p = len(chosen_predictors)
        #create design matrix
        X = np.ones(self.n)[:,None].T

        for pred in chosen_predictors:
            X = np.vstack( (X, self.dic_X[pred]) )
        if no_intercept:
            X = np.delete(X,0,0)
        #make X a matrix with predictors as horizontal vectors
        X = X.T
        X2_inv = np.linalg.inv(X.T @ X)
        #get coef matrix
        B =  X2_inv @ X.T @ self.Y
        #get predictions
        Y_pred = X@B
        #compute s
        s2 =  (1/(self.n-p)) * ((self.Y - Y_pred).T @ (self.Y - Y_pred))
        #compute SE Parameters
        SE = (s2 * X2_inv)**0.5
        #compute Y_mean
        Y_mean = np.mean(self.Y)
        #compute r2
        r2 = np.squeeze(((Y_pred -Y_mean).T @ (Y_pred -Y_mean) ) / (   (self.Y-Y_mean).T @ (self.Y -Y_mean) ))
        ##compute adjusted r2
        adj_r2 = 1-(1-r2) * (self.n -1)/(self.n-p-1)
        #compute log-likelihood (cf. https://stats.stackexchange.com/questions/87345/calculating-aic-by-hand-in-r)
        try:
            minus_two_ll = self.n*(np.log(2*np.pi)+1+np.log((np.sum((self.Y - Y_pred).T @ (self.Y - Y_pred))/self.n)))
        except:
            minus_two_ll = np.nan
        #compute AIC (cf https://stackoverflow.com/questions/35131450/calculating-bic-manually-for-lm-object)
        aic =  minus_two_ll+(p+1)*2
        #compute BIC (cf https://stackoverflow.com/questions/35131450/calculating-bic-manually-for-lm-object)
        bic = minus_two_ll+np.log(self.n)*(p+1)

        return X, B, SE, adj_r2, aic, bic

def make_time_regression_no_replicates(Y, domain = None):
    if domain is None:
        Xt = np.linspace(0,24,4,endpoint = False)
    else:
        Xt = domain
    dic_X = {'x': np.cos(w*Xt) , 'y': np.sin(w*Xt)}
    pred = ['x','y']
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            Xmat, B, SE, adj_r2, aic, bic = LinearRegression(Y = Y, dic_X = dic_X).make_complete_regression(pred)
    except:
        Xmat, B, SE, adj_r2, aic, bic = None, [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], 0,0,0

    return B, SE, adj_r2, aic, bic

def make_time_regression(Y, simple = False, predict=False):
    Xt = np.concatenate( (np.linspace(0,24,4,endpoint = False), np.linspace(0,24,4,endpoint = False), np.linspace(0,24,2,endpoint = False)))
    if simple :
        dic_X = {}
        pred = []
    else:
        dic_X = {'x': np.cos(w*Xt) , 'y': np.sin(w*Xt)}
        pred = ['x','y']

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            Xmat, B, SE, adj_r2,  aic, bic = LinearRegression(Y = Y, dic_X = dic_X).make_complete_regression(pred)
    except:
        B = np.empty((3,1,))
        B[:]=np.nan
        SE = np.empty((3,3,))
        SE[:]=np.nan
        Xmat, B, SE, adj_r2, aic, bic = None, B , SE , 0,0,0

    try:
        data = pd.DataFrame({'x': np.cos(w*Xt) , 'y': np.sin(w*Xt), 'z': Y})
        model_1 = ols("z ~ x + y", data).fit()
        model_2 = ols("z ~ 1", data).fit()
        anovaResults = anova_lm(model_2, model_1)
        pv = anovaResults.loc[1,'Pr(>F)']
    except:
        pv = np.nan

    if predict:
        X_pred = np.linspace(0,24,100)
        try:
            Y_pred = B[0]*X_pred**0 + B[1]*np.cos(w*X_pred) + B[2]*np.sin(w*X_pred)
        except:
            X_pred = [np.nan]
            Y_pred = [np.nan]
        return B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred
    else:
        return B, SE, adj_r2, aic, bic, pv

def make_space_regression(Y, predict):
    if len(Y)==16:
        X =np.concatenate( (np.linspace(0,8,8,endpoint = False), np.linspace(0,8,8,endpoint = False) ))
        vect_structure =np.array( [0]*8+[1]*8)
    else:
        X =np.concatenate( (np.linspace(0,8,8,endpoint = False), np.linspace(0,8,8,endpoint = False),np.linspace(0,8,8,endpoint = False) ))
        vect_structure =np.array( [0]*8+[1]*8+[2]*8)
    data = pd.DataFrame({'y': Y, 'mu1': X, 'mu2': 0.5*(3*X**2-1)})

    #model selection
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model, selected = forward_selected(data, 'y', vect_structure = vect_structure)
    selected = [c.split('+') for c in selected if c is not '']
    selected = [item for sublist in selected for item in sublist]
    #space_selected = [mu for mu in selected if mu=='mu1' or mu=='mu2']
    space_selected = []

    #get results
    B = model._results.params
    SE = model._results.bse
    adj_r2 = np.nan
    aic = -2*model.llf +2*( len(selected)+1 + (len(space_selected)+1)**2  )
    bic = -2*model.llf +np.log(len(vect_structure))*( len(selected)+1 + (len(space_selected)+1)**2  )



    #comparison
    #anovaResults = anova_lm(model_0, model)
    pv = np.nan#anovaResults.loc[1,'Pr(>F)']

    selected = ['mu0'] + selected
    if predict:
        dic_par = {'mu0' : 0, 'mu1' : 0, 'mu2' : 0}
        for par, val in zip(selected,B):
            dic_par[par] = val
        X_pred = np.linspace(0,8,100,endpoint = False)
        Y_pred = dic_par['mu0']*X_pred**0 + dic_par['mu1']*X_pred + dic_par['mu2']*0.5*(3*X_pred**2-1)
        return selected, B, SE, adj_r2, aic, bic, pv, X_pred, Y_pred
    else:
        return selected, B, SE, adj_r2, aic, bic, pv


def make_atger_regression(dic_atger_full):
    Xt = np.array(list(range(0,24,2))*4)
    dic_X = {'x': np.cos(w*Xt) , 'y': np.sin(w*Xt)}

    dic_reg_atg = {}
    for gene, val in dic_atger_full.items():
        Y = val.flatten('F')
        try:
            #compare oscillating model with flat model
            chosen_predictors =['x','y']
            X_osc, B_osc, SE_osc, adj_r2_osc, aic_osc, bic_osc = LinearRegression(Y, dic_X).make_complete_regression(chosen_predictors)
            chosen_predictors =[]
            X_flat, B_flat, SE_flat, adj_r2_flat, aic_flat, bic_flat = LinearRegression(Y, dic_X).make_complete_regression(chosen_predictors)
            if aic_flat<aic_osc:
                dic_reg_atg[gene] = ['flat', X_flat, B_flat, SE_flat, adj_r2_flat, aic_flat, bic_flat]
            else:
                dic_reg_atg[gene] = ['osc', X_osc, B_osc, SE_osc, adj_r2_osc, aic_osc, bic_osc]
        except:
            print('Gene ', gene, ' was discarded of the analysis in the Atger dataset')
    return dic_reg_atg

def forward_selected(data, response, time_together = False, vect_structure = None):
    remaining = set(data.columns)
    remaining.remove(response)
    selected = []

    #group candidates such that their last index is the same
    if time_together:
        for i in range(3):
            remaining.remove('a' + str(i))
            remaining.remove('b' + str(i))
            remaining.add(('a' + str(i) + '+b' + str(i)))

    try:
        flat_model = sm.MixedLM.from_formula('y~1', data, groups=vect_structure).fit(method='nm')
        #score = -2*flat_model.llf +np.log(len(vect_structure))*( 1  +1**2  )
        score = -2*flat_model.llf  +2*( 1  +1**2  )
        score = -score
    except:
        score = -10000

    current_score, best_new_score = score, score
    l_failed = []
    while remaining and current_score == best_new_score:
        scores_with_candidates = []
        for candidate in remaining:
            if (candidate=='a1+b1' or candidate=='a2+b2') and not 'a0+b0' in selected:
                candidate='a0+b0' + candidate
            formula = "{} ~ {} + 1".format(response, ' + '.join(selected + [candidate]))
            if vect_structure is None:
                score = ols(formula, data).fit().bic
            else:

                #space_selected = [mu for mu in selected + [candidate] if mu=='mu1' or mu=='mu2']
                #re_formula="{} + 1".format( ' + '.join(space_selected))
                space_selected = []
                re_formula='1'

                try:
                    model = sm.MixedLM.from_formula(formula, data, re_formula=re_formula, groups=vect_structure).fit(method='nm')
                    #score = -2*model.llf +np.log(len(vect_structure))*( len(selected + [candidate])+1 + (len(space_selected)+1)**2  )
                    score = -2*model.llf +2*( len(selected + [candidate])+1 + (len(space_selected)+1)**2  )
                    score = -score
                    if model.converged is False:
                        l_failed.append((formula, score))
                    #print("SCORE", score)
                except:
                    score = -1000000


            scores_with_candidates.append((score, candidate))
        scores_with_candidates.sort()
        best_new_score, best_candidate = scores_with_candidates.pop()
        if current_score < best_new_score:
            remaining.remove(best_candidate)
            selected.append(best_candidate)
            current_score = best_new_score

    formula = "{} ~ {} + 1".format(response, ' + '.join(selected))

    if vect_structure is None:
        model = ols(formula, data).fit()
    else:
        #space_selected = [mu for mu in selected + [candidate] if mu=='mu1' or mu=='mu2']
        space_selected = []
        re_formula='1'
        #re_formula="{} + 1".format( ' + '.join(space_selected))
        model = sm.MixedLM.from_formula(formula, data, re_formula=re_formula, groups=vect_structure).fit()

    return model, selected


if __name__ == '__main__':
    X = np.linspace(0,47,100)
    Y = X*2 + X**2 + np.random.random((100))*1000
    dic_X = {'X': X, 'X2': X**2}

    Reg = LinearRegression(Y, dic_X)

    chosen_predictors = {'X', 'X2'}
    X_d, B, SE, adj_r2, aic, bic = Reg.make_complete_regression(chosen_predictors)
    print('B',B)
    print('SE', SE)
    print('r2' ,adj_r2)
    print('aic', aic)
    print('bic', bic)
    #plt.plot(X,Y, '.')
    #plt.plot(X,X_d@B)
    #plt.show()

    dic_X['y'] = Y
    data = pd.DataFrame(dic_X)
    model = ols("y ~ X + X2", data).fit()
    mu, a, b = model._results.params
    std_mu, std_a, std_b = model._results.bse
    r2 = model.rsquared_adj
    aic = model.aic
    bic = model.bic
    print(mu, a, b)
    print(std_mu, std_a, std_b)
    print(r2)
    print('aic', aic)
    print('bic', bic)
