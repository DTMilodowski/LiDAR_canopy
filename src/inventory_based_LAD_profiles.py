# This contains code to estimate LAD profiles based on field measurements of
# tree height and crown dimensions, in addition to crown depth based on a
# regional allometric equation.
import numpy as np
from scipy import stats

pi=np.pi

# linear regression with confidence intervals and prediction intervals
def linear_regression(x_,y_,conf=0.95):
    mask = np.all((np.isfinite(x_),np.isfinite(y_)),axis=0)
    x = x_[mask]
    y = y_[mask]

    # regression to find power law exponents D = a.H^b
    m, c, r, p, serr = stats.linregress(x,y)
    x_i = np.arange(x.min(),x.max(),(x.max()-x.min())/1000.)
    y_i = m*x_i + c
    PI,PI_u,PI_l = calculate_prediction_interval(x_i,x,y,m,c,conf)
    CI,CI_u,CI_l = calculate_confidence_interval(x_i,x,y,m,c,conf)

    model_y = m*x + c
    error = y-model_y
    MSE = np.mean(error**2)

    return m, c, r**2, p, x_i, y_i, CI_u, CI_l, PI_u, PI_l

# log-log regression with confidence intervals
def log_log_linear_regression(x,y,conf=0.95):
    mask = np.all((np.isfinite(x),np.isfinite(y)),axis=0)
    logx = np.log(x[mask])
    logy = np.log(y[mask])

    # regression to find power law exponents D = a.H^b
    b, loga, r, p, serr = stats.linregress(logx,logy)
    logx_i = np.arange(logx.min(),logx.max(),(logx.max()-logx.min())/1000.)
    PI,PI_upper,PI_lower = calculate_prediction_interval(logx_i,logx,logy,b,loga,conf)

    x_i = np.exp(logx_i)
    PI_u = np.exp(PI_upper)
    PI_l = np.exp(PI_lower)

    model_logy = b*logx + loga
    error = logy-model_logy
    MSE = np.mean(error**2)
    CF = np.exp(MSE/2) # Correction factor due to fitting regression in log-space (Baskerville, 1972)
    a = np.exp(loga)
    return a, b, CF, r**2, p, x_i, PI_u, PI_l

#=================================
# ANALYTICAL CONFIDENCE AND PREDICTION INTERVALS

# Calculate confidence intervals analytically (assumes normal distribution)
# x_i = x location at which to calculate the confidence interval
# x_obs = observed x values used to fit model
# y_obs = corresponding y values
# m = gradient
# c = constant
# conf = confidence interval
# returns dy - the confidence interval
def calculate_confidence_interval(x_i,x_obs,y_obs,m,c,conf):

    alpha = 1.-conf

    n = x_obs.size
    y_mod = m*x_obs+c
    se =  np.sqrt(np.sum((y_mod-y_obs)**2/(n-2)))
    x_mean = x_obs.mean()

    # Quantile of Student's t distribution for p=1-alpha/2
    q=stats.t.ppf(1.-alpha/2.,n-2)
    dy = q*se*np.sqrt(1/float(n)+((x_i-x_mean)**2)/np.sum((x_obs-x_mean)**2))
    y_exp=m*x_i+c
    upper = y_exp+abs(dy)
    lower = y_exp-abs(dy)
    return dy,upper,lower

# Calculate prediction intervals analytically (assumes normal distribution)
# x_i = x location at which to calculate the prediction interval
# x_obs = observed x values used to fit model
# y_obs = corresponding y values
# m = gradient
# c = constant
# conf = confidence interval
# returns dy - the prediction interval
def calculate_prediction_interval(x_i,x_obs,y_obs,m,c,conf):

    alpha = 1.-conf

    n = x_obs.size
    y_mod = m*x_obs+c
    se =  np.sqrt(np.sum((y_mod-y_obs)**2/(n-2)))
    x_mean = x_obs.mean()

    # Quantile of Student's t distribution for p=1-alpha/2
    q=stats.t.ppf(1.-alpha/2.,n-2)
    dy = q*se*np.sqrt(1+1/float(n)+((x_i-x_mean)**2)/np.sum((x_obs-x_mean)**2))
    y_exp=m*x_i+c
    upper = y_exp+abs(dy)
    lower = y_exp-abs(dy)
    return dy,upper,lower

# Calculate a prediction based on a linear regression model
# As above, but this time randomly sampling from prediction interval
# m = regression slope
# c = regression interval
def random_sample_from_regression_model_prediction_interval(x_i,x_obs,y_obs,m,c,array=False):
    mask = np.all((np.isfinite(x_obs),np.isfinite(y_obs)),axis=0)
    x_obs=x_obs[mask]
    y_obs=y_obs[mask]
    n = x_obs.size
    y_mod = m*x_obs+c
    se =  np.sqrt(np.sum((y_mod-y_obs)**2/(n-2)))
    y_exp = x_i*m+c # expected value of y from model
    x_mean = x_obs.mean()
    # randomly draw quantile from t distribution (n-2 degrees of freedom for linear regression)
    if array:
        q = np.random.standard_t(n-2,size=x_i.size)
    else:
        q = np.random.standard_t(n-2)
    dy = q*se*np.sqrt(1+1/float(n)+((x_i-x_mean)**2)/np.sum((x_obs-x_mean)**2))
    y_i = y_exp+dy

    return y_i

# as above, but using log-log space (i.e. power law functions)
# a = scalar
# b = exponent
def random_sample_from_powerlaw_prediction_interval(x_i,x_obs,y_obs,a,b,array=False):
    if array:
        logy_i = random_sample_from_regression_model_prediction_interval(np.log(x_i),np.log(x_obs),np.log(y_obs),b,np.log(a),array=True)
    else:
        logy_i = random_sample_from_regression_model_prediction_interval(np.log(x_i),np.log(x_obs),np.log(y_obs),b,np.log(a))
    y_i = np.exp(logy_i)
    return y_i


#=================================
# BOOTSTRAP TOOLS

# Calculate prediction intervals through bootstrapping and resampling from residuals.
# The bootstrap model accounts for parameter uncertainty
# The residual resampling accounts for uncertainty in the residual - i.e. effects not
# accounted for by the regression model
# Inputs:
# - x_i = x location(s) at which to calculate the prediction interval
#   This should be either a numpy array or single value. nD arrays will be
#   converted to 1D arrays
# - x_obs = the observed x values used to fit model
# - y_obs = corresponding y values
# - conf = confidence interval, as fraction
# - niter = number of iterations over which to bootstrap
# - n_i = number of locations x_i (default is
# Returns:
# - ll and ul = the upper and lower bounds of the confidence interval
def calculate_prediction_interval_bootstrap_resampling_residuals(x_i,x_obs,y_obs,conf,niter):

    from matplotlib import pyplot as plt
    # some fiddles to account for likely possible data types for x_i
    n=0
    if np.isscalar(x_i):
        n=1
    else:
        try:
            n=x_i.size # deal with numpy arrays
            if x_i.ndim > 1: # linearize multidimensional arrays
                x_i=x_i.reshape(n)
        except TypeError:
            print("Sorry, not a valid type for this function")

    y_i = np.zeros((n,niter))*np.nan
    # Bootstrapping
    for ii in range(0,niter):
        # resample observations (with replacement)
        ix = np.random.choice(x_obs.size, size=n,replace=True)
        x_boot = np.take(x_obs,ix)
        y_boot = np.take(y_obs,ix)

        # regression model
        m, c, r, p, serr = stats.linregress(x_boot,y_boot)

        # randomly sample from residuals with replacement
        res = np.random.choice((y_boot-(m*x_boot + c)),size = n,replace=True)

        # estimate y based on model and randomly sampled residuals
        y_i[:,ii] = m*x_i + c + res

    # confidence intervals simply derived from the distribution of y
    ll=np.percentile(y_i,100*(1-conf)/2.,axis=1)
    ul=np.percentile(y_i,100*(conf+(1-conf)/2.),axis=1)

    return ll,ul

# equivalent to above but for log-log space prediction
def calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(x_i,x_obs,y_obs,
                                                                    conf=.9,niter=1000):
    log_ll,log_ul = calculate_prediction_interval_bootstrap_resampling_residuals(np.log(x_i),np.log(x_obs),np.log(y_obs),conf,niter)
    return np.exp(log_ll),np.exp(log_ul)

# Calculate a prediction based on a linear regression model
# As above, but this time randomly sampling from prediction interval
# calculated using random sampling from residuals.
# Note that this is intended to be used within a montecarlo framework
# m = regression slope
# c = regression interval
def random_sample_from_bootstrap_linear_regression_prediction_interval(x_i,x_obs,y_obs):

    # some fiddles to account for likely possible data types for x_i
    n=0
    if np.isscalar(x_i):
        n=1
    else:
        try:
            n=x_i.size # deal with numpy arrays
            if x_i.ndim > 1: # linearize multidimensional arrays
                x_i=x_i.reshape(n)
        except TypeError:
            print("Sorry, not a valid type for this function")

    # resample observations (with replacement) i.e. one iteration of bootstrap procedure
    ix = np.random.choice(x_obs.size, size=n,replace=True)
    x_boot = np.take(x_obs,ix)
    y_boot = np.take(y_obs,ix)

    # regression model
    m, c, r, p, serr = stats.linregress(x_boot,y_boot)

    # randomly sample from residuals with replacement
    res = np.random.choice((y_boot-(m*x_boot + c)),size = n,replace=True)

    # estimate y based on model and randomly sampled residuals
    y_i = m*x_i + c + res

    return y_i

# as above but fitting relationship in log space
def random_sample_from_bootstrap_powerlaw_prediction_interval(x_i,x_obs,y_obs):
    logy_i = random_sample_from_bootstrap_linear_regression_prediction_interval(np.log(x_i),np.log(x_obs),np.log(y_obs))
    y_i = np.exp(logy_i)
    return y_i


#================================
# INVENTORY BASED PROFILES

# This function reads in the crown allometry data from the database: Falster et al,. 2015; BAAD: a Biomass And Allometry Database for woody plants. Ecology, 96: 1445. doi: 10.1890/14-1889.1
def retrieve_crown_allometry(filename,conf=0.9):
    datatype = {'names': ('ID', 'Ref', 'Location', 'Lat', 'Long', 'Species', 'Family','Diameter','Height','CrownArea','CrownDepth'), 'formats': ('int_','S32','S256','f','f','S32','S32','f','f','f','f')}
    data = np.genfromtxt(filename, skip_header = 1, delimiter = ',',dtype=datatype)

    mask = np.all((~np.isnan(data['Height']),~np.isnan(data['CrownDepth'])),axis=0)
    H = data['Height'][mask]
    D = data['CrownDepth'][mask]
    logH = np.log(H)
    logD = np.log(D)

    # regression to find power law exponents D = a.H^b
    b, loga, r, p, serr = stats.linregress(logH,logD)
    logH_i = np.arange(logH.min(),logH.max(),(logH.max()-logH.min())/1000.)
    PI,PI_upper,PI_lower = calculate_prediction_interval(logH_i,logH,logD,b,loga,conf=conf)

    H_i = np.exp(logH_i)
    PI_u = np.exp(PI_upper)
    PI_l = np.exp(PI_lower)

    model_logD = b*logH + loga
    error = logD-model_logD
    MSE = np.mean(error**2)
    CF = np.exp(MSE/2.) # Correction factor due to fitting regression in log-space (Baskerville, 1972)
    a = np.exp(loga)

    return a, b, CF, r**2, p, H, D, H_i, PI_u, PI_l

def load_BAAD_crown_allometry_data(filename):
    datatype = {'names': ('ID', 'Ref', 'Location', 'Lat', 'Long', 'Species', 'Family',
                'Diameter','Height','CrownArea','CrownDepth'),
                'formats': ('int_','S32','S256','f','f','S32','S32','f','f','f','f')}
    data = np.genfromtxt(filename, skip_header = 1, delimiter = ',',dtype=datatype)
    mask = np.all((~np.isnan(data['Diameter']),~np.isnan(data['CrownDepth'])),axis=0)
    H = data['Height'][mask]
    D = data['CrownDepth'][mask]
    DBH = data['Diameter'][mask]
    return DBH, H, D

def load_BAAD_allometry_data(filename, filter_TropRF = True, filter_nodata = True,
                            filter_dbh = True,filter_status='None'):
    datatype = ['S10','S100','f','f','S6','f','f', 'S200','f','S100','S100','S100',
                'S4','S4','S4','S100','f','f','f','f','f','f','f','f','f','f','f',
                'f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f',
                'f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']
    data = np.genfromtxt(filename, names=True, delimiter = ',',dtype=datatype)
    if filter_TropRF:
        mask = data['vegetation']=='TropRF'
        data = data[mask]
    if filter_nodata:
        mask = np.all((np.isfinite(data['cd']),np.isfinite(data['ht'])),axis=0)
        data = data[mask]
    if filter_dbh:
        mask = np.any((data['dbh']>=0.1,data['dba']>=0.1,data['ht']>=2,np.all((data['ht']>=5,np.isnan(data['dbh']),np.isnan(data['dba'])),axis=0)),axis=0)
        data = data[mask]
    if filter_status != 'None':
        mask = data['status']==filter_status
        data=data[mask]
    return data


# Derive local allometric relationship between DBH and height -> fill gaps in census data
def load_crown_survey_data(census_file):
    datatype = {'names': ('plot','subplot','date','observers','tag','DBH','H_DBH','Height','flag','alive','C1','C2','subplotX','subplotY','density','spp','cmap_date','Xfield','Yfield','Zfield','DBH_field','Height_field','CrownArea','C3','dead_flag1','dead_flag2','dead_flag3','brokenlive_flag'), 'formats': ('S16','i8','S10','S32','i8','f','f','f','S8','i8','S132','S132','f','f','f','S64','S10','f','f','f','f','f','f','S132','i8','i8','i8','i8')}
    data = np.genfromtxt(census_file, skip_header = 1, delimiter = ',',dtype=datatype)
    return data

def calculate_allometric_equations_from_survey(data):
    # first do tree heights
    mask = np.all((~np.isnan(data['DBH_field']),~np.isnan(data['Height_field'])),axis=0)
    H = data['Height_field'][mask]
    DBH = data['DBH_field'][mask]
    logH = np.log(H)
    logDBH = np.log(DBH)

    # regression to find power law exponents H = a.DBH^b
    b_ht, loga, r, p, serr = stats.linregress(logDBH,logH)

    model_logH = b_ht*logDBH + loga
    error = logH-model_logH
    MSE = np.mean(error**2)
    CF_ht = np.exp(MSE/2) # Correction factor due to fitting regression in log-space (Baskerville, 1972)
    a_ht = np.exp(loga)

    # now do crown areas
    mask = np.all((~np.isnan(data['CrownArea']),~np.isnan(data['DBH_field']),~np.isnan(data['Height'])),axis=0)
    DBH = data['DBH_field'][mask]
    A = data['CrownArea'][mask]
    logDBH = np.log(DBH)
    logA = np.log(A)

    # regression to find power law exponents A = a.DBH^b
    b_A, loga, r, p, serr = stats.linregress(logDBH,logA)

    model_logA = b_A*logDBH + loga
    error = logA-model_logA
    MSE = np.mean(error**2)
    CF_A = np.exp(MSE/2) # Correction factor due to fitting regression in log-space (Baskerville, 1972)
    a_A = np.exp(loga)

    return a_ht, b_ht, CF_ht, a_A, b_A, CF_A


# Apply power law allometric models to estimate crown depths from heights
def calculate_crown_dimensions(DBH,Ht,Area, a_ht, b_ht, CF_ht, a_area, b_area, CF_area, a_depth, b_depth, CF_depth):
    # Gapfill record with local allometry
    # Heights
    mask = np.isnan(Ht)
    Ht[mask] = CF_ht*a_ht*DBH[mask]**b_ht
    #Crown areas
    mask = np.isnan(Area)
    Area[mask] = CF_area*a_area*DBH[mask]**b_area

    # Apply canopy depth model
    Depth = CF_depth*a_depth*Ht**b_depth

    # Remove any existing nodata values (brought forwards from input data
    mask = np.all((~np.isnan(Depth),~np.isnan(Ht),~np.isnan(Area)),axis=0)
    Depth = Depth[mask]
    Ht = Ht[mask]
    Area = Area[mask]
    return Ht, Area, Depth

# As above, but randomly sampling from prediction intervals for gapfilling
# Note that allometries are taken from local data (ref1_...) except for crown depth
# which comes from an allometic database (ref2_...).
def calculate_crown_dimensions_mc(DBH,Ht,Area,ref1_DBH,ref1_Ht,ref1_Area,ref2_DBH,ref2_Ht,ref2_D, a_ht, b_ht, a_area, b_area, a_depth, b_depth):
    # Gapfill record with local allometry
    # Heights
    mask = np.isnan(Ht)
    Ht[mask] = random_sample_from_powerlaw_prediction_interval(DBH[mask],ref1_DBH,ref1_Ht,b_ht,a_ht,array=True)
    #Ht[mask] = CF_ht*a_ht*DBH[mask]**b_ht

    #Crown areas
    mask = np.isnan(Area)
    Area[mask] = random_sample_from_powerlaw_prediction_interval(DBH[mask],ref1_DBH,ref1_Area,b_area,a_area,array=True)

    # Apply canopy depth model (from regional database)
    #Depth = random_sample_from_powerlaw_prediction_interval(DBH,ref_DBH,ref_D,b_depth,a_depth,array=True)
    Depth = random_sample_from_powerlaw_prediction_interval(Ht,ref2_Ht,ref2_D,b_depth,a_depth,array=True)

    # Remove any existing nodata values (brought forwards from input data
    mask = np.all((np.isfinite(Depth),np.isfinite(Ht),np.isfinite(Area)),axis=0)
    Depth = Depth[mask]
    Ht = Ht[mask]
    Area = Area[mask]

    return Ht, Area, Depth

# a, b, c = principal axes of the ellipse.  Initially assume circular horizontal plane i.e. a = b
# z0 = vertical position of origin of ellipse
def construct_ellipses_for_subplot(Ht, Area, Depth):
    z0 = Ht-Depth/2.
    a = np.sqrt(Area/np.pi)
    b = a.copy()
    c = Depth/2.
    #c= a.copy()
    #z0 = Ht-c

    return a, b, c, z0

# Retrieve canopy profiles based on an ellipsoidal lollipop model
# The provided ht_u vector contains the upper boundary of the canopy layers
def calculate_LAD_profiles_ellipsoid(canopy_layers, a, b, c, z0, plot_area, leafA_per_unitV=1.):
    layer_thickness = np.abs(canopy_layers[1]-canopy_layers[0])
    N_layers = canopy_layers.size
    CanopyV = np.zeros(N_layers)
    pi=np.pi
    zeros = np.zeros(a.size)
    # Formula for volume of ellipsoidal cap: V = pi*a*b*x**2*(3c-x)/c**2 where x is the vertical distance from the top of the sphere along axis c.
    # Formula for volume of ellipsoid: V = 4/3*pi*a*b*c
    for i in range(0,N_layers):
        ht_u = canopy_layers[i]
        ht_l = ht_u-layer_thickness
        mask = np.all((z0+c>=ht_l,z0-c<=ht_u),axis=0)
        x1 = np.max((z0+c-ht_u,zeros),axis=0)
        x2 = np.min((z0+c-ht_l,2*c),axis=0)
        CanopyV[i]+= np.sum(pi/3.*a[mask]*b[mask]/c[mask]**2 *(x2[mask]**2.*(3.*c[mask]-x2[mask]) - x1[mask]**2.*(3.*c[mask]-x1[mask])))

    # sanity check
    TestV = np.nansum(4*pi*a*b*c/3)
    print(CanopyV.sum(),TestV)
    LAD = CanopyV*leafA_per_unitV/plot_area
    return LAD, CanopyV

# An alternative model providing more generic canopy shapes - currently assume radial symmetry around trunk.  The crown
# volume in a given layer is determined by the volume of revolution of the function r = a*D^b
def calculate_LAD_profiles_generic(canopy_layers, Area, D, Ht, beta, plot_area, leafA_per_unitV=1.):
    r_max = np.sqrt(Area/np.pi)
    layer_thickness = np.abs(canopy_layers[1]-canopy_layers[0])
    N_layers = canopy_layers.size
    CanopyV = np.zeros(N_layers)
    pi=np.pi
    zeros = np.zeros(Ht.size)
    # Formula for volume of revolution of power law function r = alpha*D^beta:
    #                 V = pi*(r_max/D_max^beta)^2/(2*beta+1) * (D2^(2beta+1) - D1^(2beta+1))
    #                 where alpha = (r_max/D_max^beta)^2
    for i in range(0,N_layers):
        ht_u = canopy_layers[i]
        ht_l = ht_u-layer_thickness
        mask = np.all((Ht>=ht_l,Ht-D<=ht_u),axis=0)
        d1 = np.max((Ht-ht_u,zeros),axis=0)
        d2 = np.min((Ht-ht_l,D),axis=0)
        CanopyV[i]+= np.sum( pi*(r_max[mask]/D[mask]**beta)**2/(2*beta+1) * (d2[mask]**(2*beta+1) - d1[mask]**(2*beta+1.)) )

    # sanity check
    TestV = np.nansum(pi*D*r_max**2/(2*beta+1.))
    precision_requirement = 10**-8
    if CanopyV.sum() <= TestV - precision_requirement:
        print("Issue - sanity check fail: ", CanopyV.sum(),TestV)
    LAD = CanopyV*leafA_per_unitV/plot_area
    return LAD, CanopyV

def calculate_LAD_profiles_generic_mc(canopy_layers, Area, D, Ht, beta_min, beta_max,
                                        plot_area, leafA_per_unitV=1.):
    r_max = np.sqrt(Area/np.pi)
    layer_thickness = np.abs(canopy_layers[1]-canopy_layers[0])
    N_layers = canopy_layers.size
    CanopyV = np.zeros(N_layers)
    pi=np.pi
    zeros = np.zeros(Ht.size)
    D[D>Ht]=Ht[D>Ht]
    # Formula for volume of revolution of power law function r = alpha*D^beta:
    #                 V = pi*(r_max/D_max^beta)^2/(2*beta+1) * (D2^(2beta+1) - D1^(2beta+1))
    #                 where alpha = (r_max/D_max^beta)^2
    n_trees = Ht.size
    beta=np.random.rand(n_trees)*(beta_max-beta_min)+beta_min
    beta[:]=0.6
    for i in range(0,N_layers):
        ht_u = canopy_layers[i]
        ht_l = ht_u-layer_thickness
        mask = np.all((Ht>=ht_l,Ht-D<=ht_u),axis=0)
        d1 = np.max((Ht-ht_u,zeros),axis=0)
        d2 = np.min((Ht-ht_l,D),axis=0)
        CanopyV[i]+= np.sum( pi*(r_max[mask]/D[mask]**beta[mask])**2/(2*beta[mask]+1) * (d2[mask]**(2*beta[mask]+1) - d1[mask]**(2*beta[mask]+1.)) )

    # sanity check
    TestV = np.nansum(pi*D*r_max**2/(2*beta+1.))
    precision_requirement = 10**-8
    if CanopyV.sum() <= TestV - precision_requirement:
        print("Issue - sanity check fail: ", CanopyV.sum(),TestV)
    LAD = CanopyV*leafA_per_unitV/plot_area
    return LAD, CanopyV

#---------------------------------------
# SMALL SAFE PLOTS DATA
# some code to load in the survey data on small stems from the small SAFE plots.
# Numbers of stems are given as a  per-metre density
# SAFE small plots are 25 m x 25 m
def load_SAFE_small_plot_data(filename, plot_area=625.):
    datatype = {'names': ('Block', 'Plot', 'TreeID', 'DBH', 'Height', 'CrownRadius'), 'formats': ('S3','i8','S8','f16','f16','f16')}
    data = np.genfromtxt(filename, skip_header = 1, delimiter = ',',dtype=datatype)

    block_dict = {}
    blocks = np.unique(data['Block'])

    # some basic params
    bin_width = 0.5

    for bb in range(0,blocks.size):
        mask = data['Block']==blocks[bb]
        plots = np.unique(data['Plot'][mask])

        DBH = np.arange(0.,10.,bin_width)+bin_width/2.

        n_stems = np.zeros((DBH.size,plots.size))
        sum_height = np.zeros((DBH.size,plots.size))
        n_height = np.zeros((DBH.size,plots.size))
        sum_area = np.zeros((DBH.size,plots.size))
        n_area = np.zeros((DBH.size,plots.size))

        for pp in range(0,plots.size):
            mask2 = np.all((data['Block']==blocks[bb],data['Plot']==plots[pp]),axis=0)
            n_trees = mask2.sum()

            # loop through trees and only look at trees < 10 cm DBH
            for tt in range(n_trees):
                dbh = data['DBH'][mask2][tt]
                if dbh<10.:
                    ii = np.floor(dbh/bin_width)
                    n_stems[ii,pp]+=1.
                    ht = data['Height'][mask2][tt]
                    if np.isfinite(ht):
                        sum_height[ii,pp]+=ht
                        n_height[ii,pp]+=1

        # get plot means
        mean_height = sum_height/n_height

        # plot
        crown_geometry = {}
        crown_geometry['dbh'] = DBH[1:]
        crown_geometry['stem_density'] = np.mean(n_stems,axis=1)[1:]/plot_area
        crown_geometry['height'] = np.mean(mean_height,axis=1)[1:]
        block_dict[blocks[bb]]= crown_geometry

    return block_dict

# some code to get crown dimensions for stem size distributions
def calculate_crown_dimensions_small_plots(DBH,Ht,stem_density,a_ht, b_ht, CF_ht,
                                            a_area, b_area, CF_area,
                                            a_depth, b_depth, CF_depth):

    # Get rid of bins with no trees
    Ht = Ht[stem_density>0]
    DBH = DBH[stem_density>0]
    stem_density = stem_density[stem_density>0]

    # Gapfill record with local allometry
    # Heights
    mask = np.isnan(Ht)
    Ht[mask] = CF_ht*a_ht*DBH[mask]**b_ht
    #Crown areas
    Area = CF_area*a_area*DBH**b_area
    # Apply canopy depth model
    Depth = CF_depth*a_depth*Ht**b_depth

    # Remove any existing nodata values (brought forwards from input data
    mask = np.all((~np.isnan(Depth),~np.isnan(Ht),~np.isnan(Area)),axis=0)
    Depth = Depth[mask]
    Ht = Ht[mask]
    Area = Area[mask]
    return Ht, Area, Depth, stem_density

# calculate the crown profiles from the stem size distributions
def calculate_LAD_profiles_from_stem_size_distributions(canopy_layers, Area, D, Ht,
                                        stem_density, beta, leafA_per_unitV=1.):
    r_max = np.sqrt(Area/pi)
    layer_thickness = np.abs(canopy_layers[1]-canopy_layers[0])
    N_layers = canopy_layers.size
    CanopyV = np.zeros(N_layers)
    zeros = np.zeros(Ht.size)
    # Formula for volume of revolution of power law function r = alpha*D^beta:
    #                 V = pi*(r_max/D_max^beta)^2/(2*beta+1) * (D2^(2beta+1) - D1^(2beta+1))
    #                 where alpha = (r_max/D_max^beta)^2
    for i in range(0,N_layers):
        ht_u = canopy_layers[i]
        ht_l = ht_u-layer_thickness
        mask = np.all((Ht>=ht_l,Ht-D<=ht_u),axis=0)
        d1 = np.max((Ht-ht_u,zeros),axis=0)
        d2 = np.min((Ht-ht_l,D),axis=0)
        CanopyV[i]+= np.sum( pi*(r_max[mask]/D[mask]**beta)**2/(2*beta+1) *
                            (d2[mask]**(2*beta+1) - d1[mask]**(2*beta+1))*stem_density[mask] )

    LAD = CanopyV*leafA_per_unitV
    return LAD

# as above for ellipsoidal geometry
def calculate_LAD_profiles_ellipsoid_from_stem_size_distributions(canopy_layers,
                                    Area, D, Ht, stem_density, leafA_per_unitV=1.):
    # ellipsoid dims
    a = np.sqrt(Area/pi)
    b=a.copy()
    c = D/2.
    z0 = Ht-c

    layer_thickness = np.abs(canopy_layers[1]-canopy_layers[0])
    N_layers = canopy_layers.size
    CanopyV = np.zeros(N_layers)
    zeros = np.zeros(a.size)
    # Formula for volume of ellipsoidal cap:
    #       V = 1/3*pi*a*b*x**2*(3c-x)/c**2
    # where x is the vertical distance from the top of the ellipsoid along axis c.
    # Therefore ellipsoidal slice between x1 and x2:
    #       V = (1/3*pi*a*b*/c**2)*(x2**2*(3c-x2)-x2**2*(3c-x2))
    for i in range(0,N_layers):
        ht_u = canopy_layers[i]
        ht_l = ht_u-layer_thickness
        mask = np.all((z0+c>=ht_l,z0-c<=ht_u),axis=0)
        x1 = np.max((z0+c-ht_u,zeros),axis=0)
        x2 = np.min((z0+c-ht_l,2*c),axis=0)
        CanopyV[i]+= np.sum( pi/3.*a[mask]*b[mask]/c[mask]**2 *
                            (x2[mask]**2.*(3.*c[mask]-x2[mask])-
                            x1[mask]**2.*(3.*c[mask]-x1[mask]))*stem_density[mask] )

    LAD = CanopyV*leafA_per_unitV
    return LAD

#=====================================================================================================================
# Load in data for SAFE detailed subplot census - all trees >2 cm  -
# only interested in a subset of the fields
def load_SAFE_small_stem_census(filename, sp_area=20.**2, N_subplots = 25):

    datatype = {'names': ('Plot', 'Subplot', 'Date', 'Obs','tag', 'DBH_ns',
                         'DBH_ew', 'H_POM', 'Height', 'Flag', 'Notes'),
                'formats': ('S8','i8','S10','S64','int_','f',
                        'f','f','f','S4','S32')}
    data = np.genfromtxt(filename, skip_header = 1, usecols=np.arange(0,11),
                        delimiter = ',',dtype=datatype)

    data['DBH_ns']/=10. # convert from mm to cm
    data['DBH_ew']/=10. # convert from mm to cm

    # loop through data & remove lianas
    N = data['Plot'].size

    mask = np.ones(N,dtype='bool')
    for i in range(0,N):
        if(b'liana' in data['Notes'][i]):
            mask[i] = False
        elif(b'Liana' in data['Notes'][i]):
            mask[i] = False

    # Remove lianas and dead trees
    data = data[mask]
    data = data[data['Flag']!='dead']
    data = data[data['Flag']!='dead, broken']
    data = data[data['Flag']!='liana']

    plot_dict = {}
    plots = np.unique(data['Plot'])

    for pp in range(plots.size):
        plot_data = data[data['Plot']==plots[pp]]
        subplots = np.unique(plot_data['Subplot'])

        # some basic params
        bin_width = 0.5
        DBH = np.arange(0.,10.,bin_width)+bin_width/2.
        n_stems_i = np.zeros((DBH.size,subplots.size))
        n_stems = np.zeros((DBH.size,N_subplots))


        stem_dict = {}
        sp_present = np.zeros(N_subplots,dtype='bool')

        for ss in range(0,subplots.size):
            sp_present[subplots[ss]-1] = True

            mask = plot_data['Subplot']==subplots[ss]
            n_trees = mask.sum()

            sp_data = plot_data[mask]

            # loop through trees and only look at trees < 10 cm DBH
            for tt in range(n_trees):
                dbh = (sp_data['DBH_ns'][tt]+sp_data['DBH_ew'][tt])/2.
                if dbh<10.:
                    ii = np.floor(dbh/bin_width).astype('int')
                    n_stems_i[ii,ss]+=1.
                    n_stems[ii,subplots[ss]-1]+=1.

        average = np.mean(n_stems_i,axis=1)
        for ss in range(0,N_subplots):
            if ~sp_present[ss]:
                n_stems[:,ss] = average

        stem_dict['dbh'] = DBH[1:]
        stem_dict['stem_density'] = n_stems[1:,:]/sp_area
        plot_dict[plots[pp]]=stem_dict

    return plot_dict

#=====================================================================================================================
# Load in data for Danum detailed census - every subplot censused
def load_Danum_stem_census(filename, sp_area=20.**2):

    datatype = {'names': ('Plot', 'Subplot', 'Date', 'x', 'y','tag','spp', 'Nstem', 'DBH1', 'Codes', 'Notes', 'DBH2', 'DBH3', 'DBH4'), 'formats': ('S6','int_','S10','int_','int_','S6','S6','int_','f','S8','S32','f','f','f')}
    data = np.genfromtxt(filename, skip_header = 1, delimiter = ',',dtype=datatype)
    data['DBH1']/=10. # convert from mm to cm
    data['DBH2']/=10. # convert from mm to cm
    data['DBH3']/=10. # convert from mm to cm
    data['DBH4']/=10. # convert from mm to cm
    subplots = np.unique(data['Subplot'])
    N_subplots = subplots.size

    stem_dict = {}

    # some basic params
    bin_width = 0.5
    DBH = np.arange(0.,10.,bin_width)+bin_width/2.
    n_stems = np.zeros((DBH.size,N_subplots))

    for ss in range(0,N_subplots):
        mask = data['Subplot']==subplots[ss]
        n_trees = mask.sum()

        sp_data = data[mask]

        # loop through trees and only look at trees < 10 cm DBH
        for tt in range(n_trees):
            dbh = sp_data['DBH1'][tt]
            if dbh<10.:
                ii = np.floor(dbh/bin_width).astype('int')
                n_stems[ii,subplots[ss]-1]+=1.

        # loop through 2nd stems
        m2 = sp_data['DBH2']>0
        n_trees2 = np.sum(m2)
        for tt in range(n_trees2):
            dbh = sp_data['DBH2'][m2][tt]
            if dbh<10.:
                ii = np.floor(dbh/bin_width).astype('int')
                n_stems[ii,subplots[ss]-1]+=1.

        # loop through 3rd stems
        m3 = sp_data['DBH3']>0
        n_trees3 = np.sum(m3)
        for tt in range(n_trees3):
            dbh = sp_data['DBH3'][m3][tt]
            if dbh<10.:
                ii = np.floor(dbh/bin_width).astype('int')
                n_stems[ii,subplots[ss]-1]+=1.

        # loop through 4th stems
        m4 = sp_data['DBH4']>0
        n_trees4 = np.sum(m4)
        for tt in range(n_trees4):
            dbh = sp_data['DBH4'][m4][tt]
            if dbh<10.:
                ii = np.floor(dbh/bin_width).astype('int')
                n_stems[ii,subplots[ss]-1]+=1.


    stem_dict['dbh'] = DBH[1:]
    stem_dict['stem_density'] = n_stems[1:,:]/sp_area
    return stem_dict


# calculate crown geometries based on stem distributions only
def calculate_crown_dimensions_for_stem_distributions(DBH,stem_density,a_ht, b_ht, CF_ht,
                                    a_area, b_area, CF_area, a_depth, b_depth, CF_depth):

    # Get rid of bins with no trees
    DBH = DBH[stem_density>0]
    stem_density = stem_density[stem_density>0]

    # Gapfill record with local allometry
    # Heights
    Ht = CF_ht*a_ht*DBH**b_ht
    #Crown areas
    Area = CF_area*a_area*DBH**b_area
    # Apply canopy depth model
    Depth = CF_depth*a_depth*Ht**b_depth

    # Remove any existing nodata values (brought forwards from input data)
    mask = np.all((~np.isnan(Depth),~np.isnan(Ht),~np.isnan(Area)),axis=0)
    Depth = Depth[mask]
    Ht = Ht[mask]
    Area = Area[mask]
    return Ht, Area, Depth, stem_density


#=============================================================================================================================
# 3-dimensional crown models
# create 3D model of individual crown within a specified environment.
# Voxels are included if the centre is contained within the crown.
# This imposes the following constraints:
# - The voxel is in crown if:
#              1) 0 <= z' <= Zmax
#              2) r <= Rmax/Zmax^beta * z'^beta;
#   where r = sqrt((x-x0)**2+(y-y0)**2)
#
# Note that we do not have data on crowns for trees outwith each plot
# that overlap into the plot area. To compensate, we include in our
# plot-level estimates the full crown volume within each tree inside
# the plot, even if these overlap outside the plot bounds. This assumes
# that overlap from trees outwith the plot more or less equals overlap
# from trees inside the plot beyond the plot footprint.
#
# Input variables:
# - canopy_matrix (the three dimensional matrix representing the canopy space - dimensions x,y,z)
# - x,y,z (vectors containing the centre coordinates of each voxel
# - x0,y0 (horizontal coordinates of stem)
# - Z0 (the depth from the top of the domain to the tree top
# - Zmax (the maximum crown depth)
# - Rmax (the maximum radius of the tree)
# - beta (the exponent controlling the canopy morphology)
# Returns:
# - 0 (the canopy_matrix array is updated with the new crown)
def generate_3D_crown(canopy_matrix,x,y,z,x0,y0,Z0,Zmax,Rmax,beta):
    #from matplotlib import pyplot as plt

    # generate masks for tree crown
    xm,ym,zm = np.meshgrid(x,y,z)
    r = np.sqrt((xm-x0)**2+(ym-y0)**2)
    con = np.all((zm>=Z0,zm<=Zmax+Z0,r <= ((zm-Z0)**beta) * Rmax/(Zmax**beta)),axis=0)
    #crown = np.zeros(canopy_matrix.shape)
    canopy_matrix[con] = 1
    return 0
    #return crown

# alternative is to use ellipsoid crowns, a similar approach has been used by other
# researchers. Voxels can readily be assigned to ellpsoidal crowns based on the
# inequality:
#      1 >= ((x-x0)/a)^2 + ((y-y0)/b)^2 + ((z-z0)/c)^2
# where:
# - x,y,z = coordinates of voxel
# - x0,y0 = trunk location in x and y
# - z0 = Ht-CrownDepth/2
# - a=b = radius of crown, R
# - c = CrownDepth/2
#
# Input variables:
# - canopy_matrix (the three dimensional matrix representing the canopy space - dimensions x,y,z)
# - x,y,z (vectors containing the centre coordinates of each voxel
# - x0,y0 (horizontal coordinates of stem)
# - H (the height of the top of the tree)
# - D (the maximum crown depth)
# - R (the maximum radius of the tree crown)
# Returns:
# - 0 (the canopy_matrix array is updated with the new crown)

def generate_3D_ellipsoid_crown(canopy_matrix,xm,ym,zm,x0,y0,H,D,R):
    z0 = H-D/2.
    # generate masks for tree crown
    #xm,ym,zm = np.meshgrid(x,y,z)
    #con = (((xm-x0)/R)**2+((ym-y0)/R)**2+((zm-z0)/(D/2.))**2) <= 1
    #canopy_matrix[con]=1
    canopy_matrix += ((((xm-x0)/R)**2+((ym-y0)/R)**2+((zm-z0)/(D/2.))**2) <= 1)
#

# 3-D canopy
# This function creates 3D crown model by aggregating individual crowns
# constructed using the above function. It takes in the following input
# variables:
# - x,y,z (vectors containing the centre coordinates of each voxel)
# - x0,y0 (vectors indicating the relative horizontal coordinates of the
#   surveyed stems)
# - Z0 (vector containing the depth from the top of the domain to the each
#   tree top
# - Zmax (vector containing the maximum crown depth for each tree surveyed)
# - Rmax (vector containing the maximum radius of the tree)
# - beta (vector containing the exponents controlling the canopy morphology)
# - plot_mask (an optional mask that accounts for non-square plot geometries)
# - buffer (an optional argument that by default is zero, but should be
#   increased so that it is sufficient to account for crown overlap
def generate_3D_canopy(x,y,z,x0,y0,Z0,Zmax,Rmax,beta):
    #from matplotlib import pyplot as plt
    # first create buffer
    n_trees = x0.size
    dx = x[1]-x[0]
    dy = y[1]-y[0]

    # now create buffered canopy matrix
    #crowns = np.zeros((n_trees,y.size,x.size,z.size),dtype='float')
    canopy = np.zeros((y.size,x.size,z.size),dtype='float')
    # Now loop through the trees. For each tree, calculate the crown volume,
    # then add to the crown map.
    for tt in range(0,n_trees):
        generate_3D_crown(canopy,xm,ym,zm,x0[tt],y0[tt],Z0[tt],Zmax[tt],Rmax[tt],beta[tt])

    #plt.imshow(np.transpose(np.sum(canopy,axis=1)),origin='lower');plt.colorbar();plt.show()

    return canopy

# alternative function that this time uses ellipsoid crowns
# It takes in the following input variables:
# - x,y,z (vectors containing the centre coordinates of each voxel)
# - x0,y0 (vectors indicating the relative horizontal coordinates of the
#   surveyed stems)
# - H (vector containing the Heights of each tree)
# - D (vector containing the crown depth for each tree surveyed)
# - R (vector containing the maximum radius of the tree)
# - plot_mask (an optional mask that accounts for non-square plot geometries)
# - buffer (an optional argument that by default is zero, but should be
#   increased so that it is sufficient to account for crown overlap
def generate_3D_ellipsoid_canopy(x,y,z,x0,y0,H,D,R):
    # first create buffer
    n_trees = x0.size
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    canopy = np.zeros((y.size,x.size,z.size),dtype='float')
    xm,ym,zm = np.meshgrid(x,y,z)

    # Now loop through the trees. For each tree, calculate the crown volume,
    # then add to the crown map.
    for tt in range(0,n_trees):
        if(tt%50==0):
            print('\tprocessing tree %i from %i' % (tt,n_trees))
        generate_3D_ellipsoid_crown(canopy,xm,ym,zm,x0[tt],y0[tt],H[tt],D[tt],R[tt])
    canopy[canopy>1]=1
    return canopy


#===============================================================================
# MONTE CARLO ROUTINES FOR CROWN CONTSTRUCTION
# First version samples from the error distribution simulated from the
# allometric rlationships underpinning the modelled crown geometry (ellipsoid)
def calculate_crown_volume_profiles_mc(x,y,z,x0,y0,Ht,DBH,Area,
                                        a_ht,b_ht,a_A,b_A,a_D,b_D,
                                        field_data,BAAD_data,n_iter=10):
    profiles = np.zeros((n_iter,z.size))
    for ii in range(0,n_iter):
        print('iteration %i out of %i' % (ii,n_iter))
        # now get field inventory estimate
        # Note that we only deal with the 1ha plot level estimates as errors relating stem based
        # vs. area based are problematic at subplot level
        Ht[np.isnan(Ht)] = random_sample_from_powerlaw_prediction_interval(DBH[np.isnan(Ht)],
                                                field_data['DBH_field'],field_data['Height_field'],
                                                a_ht,b_ht,array=True)
        Area[np.isnan(Area)] = random_sample_from_powerlaw_prediction_interval(DBH[np.isnan(Area)],
                                                field_data['DBH_field'],field_data['CrownArea'],
                                                a_A,b_A,array=True)
        Depth = random_sample_from_powerlaw_prediction_interval(Ht,BAAD_data['Ht'],BAAD_data['D'],
                                                a_D,b_D,array=True)
        Rmax = np.sqrt(Area/np.pi)

        crown_model = generate_3D_ellipsoid_canopy(x,y,z,x0,y0,Ht,Depth,Rmax)

        profiles[ii,:] = np.sum(np.sum(crown_model,axis=1),axis=0)/10.**4

    return profiles

# second version adds measurement error. Errors are two part lists or arrays
# indicating bias and random error (expressed as an estimated fraction)
def calculate_crown_volume_profiles_mc_with_measurement_error(x,y,z,x0,y0,Ht_,DBH_,Area_,
                                        a_ht,b_ht,a_A,b_A,a_D,b_D,
                                        error,field_data,BAAD_data,n_iter=10):
    profiles = np.zeros((n_iter,z.size))
    for ii in range(0,n_iter):
        # combine random and systematic errors to the observations as fractions
        err_Ht = (np.random.normal(scale = error['Ht'][1],size = Ht.size)+error['Ht'][0])
        err_DBH = (np.random.normal(scale = error['DBH'][1],size = DBH.size)+error['DBH'][0])
        err_Area = (np.random.normal(scale = error['Area'][1],size = Area.size)+error['Area'][0])

        # apply errors
        DBH = DBH_*(1+err_DBH)
        Ht = Ht_*(1+err_Ht_)
        DBH = Area_*(1+err_Area)

        # Randomly sample allometrics from simulated error distribution
        Ht[np.isnan(Ht)] = random_sample_from_powerlaw_prediction_interval(DBH[np.isnan(Ht)],
                                                field_data['DBH_field'],field_data['Height_field'],
                                                a_ht,b_ht,array=True)
        Area[np.isnan(Area)] = random_sample_from_powerlaw_prediction_interval(DBH[np.isnan(Area)],
                                                field_data['DBH_field'],field_data['CrownArea'],
                                                a_A,b_A,array=True)
        Depth = random_sample_from_powerlaw_prediction_interval(Ht,BAAD_data['Ht'],BAAD_data['D'],
                                                a,b,array=True)
        Rmax = np.sqrt(Area/np.pi)
        crown_model = field.generate_3D_ellipsoid_canopy(x,y,z,x0,y0,Ht,Depth,Rmax)

        profiles[ii,:] = np.sum(np.sum(crown_model,axis=1),axis=0)/10.**4

    return profiles
