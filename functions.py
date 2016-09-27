import numpy as np
import pandas as pd

def autocorr_test(_xdata, _ydata):
    from statsmodels.stats.diagnostic import acorr_ljungbox
    from statsmodels.tsa.stattools import acf
    
#all statst need regularly spaced, continuous time series - just y variable
#Durbin-Watson statistics:
# calculated correctly with missing data
# but no significance level. Apparently critical values for DW are not implemented in any python library
#ACF:
# crashes on missing data
# Ljung-Box:
# crashes on missing data too
    _ydata=np.ma.masked_invalid(_ydata)
    #autocorrelation in residuals
    #this is acf function that does not allow nans
#    print "\nautocorrelation for first three lags:", acf(_ydata)[1:4]
    #this is from pandas, is nan agnostic
    pdf=pd.Series(_ydata, index=_xdata, copy=True)
    print "autocorrelation for first three lags:", [pdf.autocorr(i) for i in range(1,4)]
    #durbin-watson
    a=_ydata[:-1].astype('float')
    b=_ydata[1:].astype('float')
    _stat=np.nansum((b-a)**2)/np.nansum(_ydata**2)
    print "Durbin-Watson statistic (close to 2 if no autocorrelation):", _stat
    _stat, _pvalue=acorr_ljungbox(_ydata, lags=1, boxpierce=False)    
    print "Ljung-Box p-value on lag 1 autocorrelation:", _pvalue
    print ""

def trend_CI(x_var, y_var, n_boot=1000, ci=95, trendtype="linreg", q=0.5, frac=0.6, it=3, autocorr=None, CItype="bootstrap"):
    """calculates bootstrap confidence interval and significance level for trend, ignoring autocorrelation or accounting for it
    Parameters
    ----------
    x_var : list
      independent variable
    y_var : list
      dependent variable, same length as x_var
    q : int, optional, only if trendtype==quantreg
      quantile for which regression is to be calculated
    n : int, optional
      number of bootstrap samples
    ci : int, optional
      confidence level. Default is for 95% confidence interval
    frac : int, optional, only if trendtype==lowess
      lowess parameter (fraction of time period length used in local regression)
    it : int, optional, only if trendtype==lowess
      lowess parameter (numbre of iterations)
    autocorr : str, optional
      way of accounting for autocorrelation, possible values: None, "bootstrap"
    trendtype : str, optional
      method of trend derivation, possible values: lowess, linreg, quantreg, TheilSen
    CItype : str, optional
      method of CI derivation, possible values: "analytical" and "bootstrap". 
      if trendtype is "lowess", CItype will be set to None
      if CItype is "analytical": autocorrelation will be set to None
      

    Results
    -------
    returns library with following elements:
    slope - slope of the trend
    CI_high - CI on the slope value
    CI_low - as above
    pvalue - trend's significance level
    trend - trend line, or rather its y values for all x_var
    trendCI_high - confidence interval for each value of y
    trendCI_low - as above

    Remarks
    -------
    the fit function ocassionally crashes on resampled data. The workaround is to use try statement
    """
    #for linreg
    import statsmodels.api as sm
    from statsmodels.regression.linear_model import OLS
    #for arima
    import statsmodels.tsa as tsa
    #for quantreg
    import statsmodels.formula.api as smf
    from statsmodels.regression.quantile_regression import QuantReg
    #for lowess
    import statsmodels.nonparametric.api as npsm
    #other
    from statsmodels.distributions.empirical_distribution import ECDF
    from scipy.stats import mstats, mannwhitneyu, t, kendalltau
    from arch.bootstrap import StationaryBootstrap, IIDBootstrap

    #preparing data
    if CItype=="analytical" and trendtype=="TheilSen":
        CItype="bootstrap"
    x_var=np.array(x_var)
    y_var=np.ma.masked_invalid(y_var)
    n_data=len(y_var)
    ci_low=(100-ci)/2
    ci_high=100-ci_low
    
    #setting bootstrapping function
    if autocorr=="bootstrap":
        bs=StationaryBootstrap(3, np.array(range(len(y_var))))
    else:
        bs=IIDBootstrap(np.array(range(len(y_var))))
    
    if trendtype=="quantreg":
        print "Quantile regression, CI type: "+CItype+", autocorrelation adjustment: "+str(autocorr)+"\n"
        xydata=pd.DataFrame(np.column_stack([x_var, y_var]), columns=['X', 'Y'])
        model=smf.quantreg('Y ~ X', xydata)
        res=model.fit(q=q)
        intcpt=res.params.Intercept
        slope=res.params.X
        pvalue=res.pvalues[1]
        CI_low=res.conf_int()[0]['X']
        CI_high=res.conf_int()[1]['X']
        y_pred=res.predict(xydata)
        #calculating residuals
        resids=y_var-y_pred
        #calculate autocorrelation indices
        autocorr_test(x_var, resids)
            
        if CItype=="bootstrap":
            #bootstrapping
            bs_trends=np.copy(y_pred).reshape(-1,1)
            bs_slopes=[]
            bs_intcpts=[]
            for data in bs.bootstrap(n_boot):
                ind=data[0][0]
                model = smf.quantreg('Y ~ X', xydata.ix[ind,:])
                try:
                    res = model.fit(q=q)
                    bs_slopes=bs_slopes+[res.params.X]
                    bs_intcpts=bs_intcpts+[res.params.Intercept]
                    bs_trends=np.append(bs_trends,res.predict(xydata).reshape(-1,1), 1)
                except:
                    goingdownquietly=1
    if trendtype=="linreg":
        print "Linear regression, CI type: "+CItype+", autocorrelation adjustment: "+str(autocorr)+"\n"
        x_varOLS = sm.add_constant(x_var)
        model = sm.OLS(y_var, x_varOLS, hasconst=True, missing='drop')
        res = model.fit()
        intcpt,slope=res.params
        pvalue=res.pvalues[1]
        CI_low,CI_high=res.conf_int()[1]
        y_pred=res.predict(x_varOLS)
        #calculating residuals
        resids=y_var-y_pred
        #calculate autocorrelation indices
        autocorr_test(x_var, resids)
        
        if CItype=="bootstrap":        
            #bootstrapping for confidence intervals
            bs_slopes=[]
            bs_intcpts=[]
            bs_trends=np.copy(y_pred).reshape(-1,1)
            for data in bs.bootstrap(n_boot):
                ind=data[0][0]
                model = sm.OLS(y_var[ind], x_varOLS[ind,:], hasconst=True, missing='drop')
                try:
                    res = model.fit()
                    bs_slopes=bs_slopes+[res.params[1]]
                    bs_intcpts=bs_intcpts+[res.params[0]]
                    bs_trends=np.append(bs_trends,res.predict(x_varOLS).reshape(-1,1), 1)
                except:
                    goingdownquietly=1
                    
    if trendtype=="TheilSen":
#        print "Theil-Sen slope, CI type: "+CItype+", autocorrelation adjustment: "+str(autocorr)+"\n"
        #significance of MK tau
        tau,pvalue=kendalltau(x_var, y_var)
#        print "raw MK tau:", tau, "raw MK pvalue:", pvalue
        #TS slope and confidence intervals
        slope,intercept,CI_low,CI_high=mstats.theilslopes(y_var, x_var, alpha=0.95)        
        #getting slope line's y values
        y_pred=intercept+slope*x_var
        #calculating residuals
        resids=y_var-y_pred
        #calculate autocorrelation indices
        autocorr_test(x_var, resids)
                    
        if CItype=="bootstrap":
            #bootstrapping for confidence intervals
            bs_slopes=[]
            bs_intcpts=[]
            bs_trends=np.copy(y_pred).reshape(-1,1)
            for data in bs.bootstrap(n_boot):
                ind=data[0][0]
                res=mstats.theilslopes(y_var[ind], x_var[ind], alpha=0.95)
                bs_slopes=bs_slopes+[res[0]]
                bs_intcpts=bs_intcpts+[res[1]]
                bs_trends=np.append(bs_trends, (res[1]+res[0]*x_var).reshape(-1,1), 1)

    if trendtype=="lowess":
        print "Lowess\n"
        temp=dict(npsm.lowess(y_var, x_var, frac=frac, it=it, missing="drop"))
        y_pred=np.array(map(temp.get, x_var)).astype("float").reshape(-1,1)
        bs_trends=np.copy(y_pred)
        
        for data in bs.bootstrap(n_boot):
            ind=data[0][0]
            try:
                temp = dict(npsm.lowess(y_var[ind], x_var[ind], frac=frac, it=it, missing="drop"))
                temp=np.array(map(temp.get, x_var)).astype("float").reshape(-1,1)
                pred=pd.DataFrame(temp, index=x_var)
                temp_interp=pred.interpolate().values
                bs_trends=np.append(bs_trends, temp_interp, 1)
            except:
                goingdownquietly=1


    #calculating final values of CI and p-value

    #skipping when lowess
    if trendtype=="lowess":
        CI_low=np.nan
        CI_high=np.nan
        slope=np.nan
        intcpt=np.nan
        pvalue=np.nan
        confint=np.nanpercentile(bs_trends, [ci_low,ci_high], 1)
        trendCI_low=confint[:,0]
        trendCI_high=confint[:,1]
    else:
        if CItype=="bootstrap":
            #values for slope, intercept and trend can be obtained as medians of bootstrap distributions, but normally analytical parameters are used instead
            # it the bootstrap bias (difference between analytical values and bootstap median) is strong, it might be better to use bootstrap values. 
            # These three lines would need to be uncommented then
#            slope=np.median(bs_slopes)
#            intcpt=np.median(bs_intcpts)
#            trend=intcpt+slope*x_var
            #these are from bootstrap too, but needs to be used for this accounts for autocorrelation, which is the point of this script
            CI_low,CI_high=np.percentile(bs_slopes, [5, 95])                
            ecdf=ECDF(bs_slopes)
            pvalue=ecdf(0)
            #this makes sure we are calculating p-value on the correct side of the distribution. That will be one-sided pvalue
            if pvalue>0.5:
                pvalue=1-pvalue
            confint=np.nanpercentile(bs_trends, [ci_low,ci_high], 1)
            print "bs_trends:", bs_trends.shape, confint.shape
            trendCI_low=confint[:,0]
            trendCI_high=confint[:,1]
        else:
            #this is for analytical calculation of trend confidence interval
            #it happens in the same way for each of the trend types, thus it is done here, not under the trendtype subroutines
            #making sure x are floats
            xtemp=np.array(x_var)*1.0
            #squared anomaly
            squanom=(xtemp-np.mean(xtemp))**2
            temp=((1./len(x_var))+(squanom/sum(squanom)))**0.5
            #standard error of estmation
            see=(np.nansum((np.array(y_var)-np.nanmean(y_pred))**2)/len(x_var))**0.5
            #adjusting ci
            ci_adj=1-((1-ci/100.)/2)
            #accounting for uncertainty in mean through student's t
            tcomp=t.ppf(ci_adj, len(x_var)-2)
            #confidence interval
            cint=tcomp*see*temp
            #for trend only
            trendCI_high=y_pred+cint
            trendCI_low=y_pred-cint

        print trendtype, "slope:",slope, "pvalue (one sided):", pvalue, "conf interval:", CI_low, CI_high, "autocorrelation adjustment:", autocorr, "\n"
    output={"slope":slope, "CI_high":CI_high, "CI_low":CI_high, "pvalue":pvalue, "trend": y_pred, "trendCI_low":trendCI_low, "trendCI_high":trendCI_high}
    return output

def get_TheilSen(_x, what, _nboot, _y):
    #the x y are weird, it appears that apply passes the dataframe column as last element
    from arch.bootstrap import StationaryBootstrap, IIDBootstrap
    from scipy.stats import mstats, mannwhitneyu, t, kendalltau
    from statsmodels.distributions.empirical_distribution import ECDF
    try:
        if what=="slope":
            return mstats.theilslopes(np.ma.masked_invalid(_y.values), _x)[0]*86400*365*1000000000
        elif what=="pval_tau":
            return kendalltau(_x, _y)[1]/2
        elif what=="pval_autocorr":            
            res0=mstats.theilslopes(_y, _x, alpha=0.95)[0]
            bs=StationaryBootstrap(3, np.array(range(len(_y))))
            bs_slopes=[]
            for data in bs.bootstrap(_nboot):
                ind=data[0][0]
                res=mstats.theilslopes(_y[ind], _x, alpha=0.95)
                bs_slopes=bs_slopes+[res[0]]
            ecdf=ECDF(bs_slopes)
            pvalue=ecdf(res0)
            if pvalue>0.5:
                pvalue=1-pvalue
#            print pvalue
            return pvalue
        elif what=="pval":
            bs=IIDBootstrap(np.array(range(len(_y))))
            bs_slopes=[]
            for data in bs.bootstrap(_nboot):
                ind=data[0][0]
                res=mstats.theilslopes(_y[ind], _x, alpha=0.95)
                bs_slopes=bs_slopes+[res[0]]
            ecdf=ECDF(bs_slopes)
            pvalue=ecdf(0)
            if pvalue>0.5:
                pvalue=1-pvalue
#            print pvalue
            return pvalue
    except:
        return np.nan


