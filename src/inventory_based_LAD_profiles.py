# This contains code to estimate LAD profiles based on field measurements of tree height and crown dimensions, in addition to crown depth based on a regional allometric equation.

import numpy as np
from scipy import stats

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
# ANALYTICAL PREDICTION INTERVALS

# Calculate prediction intervals analytically (assumes normal distribution)
# x_i = x location at which to calculate the prediction interval
# x_obs = observed x values used to fit model
# y_obs = corresponding y values
# m = gradient
# c = constant
# conf = confidence interval
# returns dy - the confidence interval
def calculate_prediction_interval(x_i,x_obs,y_obs,m,c,conf):
    
    alpha = 1.-conf
    
    n = x_obs.size
    y_mod = m*x_obs+c
    se =  np.sqrt(np.sum((y_mod-y_obs)**2/(n-2)))
    x_mean = x_obs.mean()
    
    # Quantile of Student's t distribution for p=1-alpha/2
    q=stats.t.ppf(1.-alpha/2.,n-2)
    dy = q*se*np.sqrt(1+1/n+((x_i-x_mean)**2)/np.sum((x_obs-x_mean)**2))
    y_exp=m*x_i+c
    upper = y_exp+abs(dy)
    lower = y_exp-abs(dy)
    return dy,upper,lower 

# Calculate a prediction based on a linear regression model
# As above, but this time randomly sampling from prediction interval
# m = regression slope
# c = regression interval
def random_sample_from_regression_model_prediction_interval(x_i,x_obs,y_obs,m,c,array=False):
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
    dy = q*se*np.sqrt(1+1/n+((x_i-x_mean)**2)/np.sum((x_obs-x_mean)**2))
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
            print "Sorry, not a valid type for this function"

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

        # estimate y based on model and randomly sampled residuals
        y_i[:,ii] = m*x_i + c + res

    # confidence intervals simply derived from the distribution of y
    ll=np.percentile(y_i,100*(1-conf)/2.,axis=1)
    ul=np.percentile(y_i,100*(conf+(1-conf)/2.),axis=1)

    return ll,ul

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
            print "Sorry, not a valid type for this function"

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
    datatype = {'names': ('ID', 'Ref', 'Location', 'Lat', 'Long', 'Species', 'Family','Diameter','Height','CrownArea','CrownDepth'), 'formats': ('i8','S32','S256','f16','f16','S32','S32','f16','f16','f16','f16')}
    data = np.genfromtxt(filename, skiprows = 1, delimiter = ',',dtype=datatype)
    
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
    CF = np.exp(MSE/2) # Correction factor due to fitting regression in log-space (Baskerville, 1972)
    a = np.exp(loga)

    return a, b, CF, r**2, p, H, D, H_i, PI_u, PI_l

def load_BAAD_crown_allometry_data(filename):
    datatype = {'names': ('ID', 'Ref', 'Location', 'Lat', 'Long', 'Species', 'Family','Diameter','Height','CrownArea','CrownDepth'), 'formats': ('i8','S32','S256','f16','f16','S32','S32','f16','f16','f16','f16')}
    data = np.genfromtxt(filename, skiprows = 1, delimiter = ',',dtype=datatype)
    mask = np.all((~np.isnan(data['Diameter']),~np.isnan(data['CrownDepth'])),axis=0)
    H = data['Height'][mask]
    D = data['CrownDepth'][mask]
    DBH = data['Diameter'][mask]
    return DBH, H, D

# Derive local allometric relationship between DBH and height -> fill gaps in census data
def load_crown_survey_data(census_file):
    datatype = {'names': ('plot','subplot','date','observers','tag','DBH','H_DBH','Height','flag','alive','C1','C2','subplotX','subplotY','density','spp','cmap_date','Xfield','Yfield','Zfield','DBH_field','Height_field','CrownArea','C3','dead_flag1','dead_flag2','dead_flag3','brokenlive_flag'), 'formats': ('S16','i8','S10','S32','i8','f8','f8','f8','S8','i8','S132','S132','f8','f8','f8','S64','S10','f8','f8','f8','f8','f8','f8','S132','i8','i8','i8','i8')}
    data = np.genfromtxt(census_file, skiprows = 1, delimiter = ',',dtype=datatype)
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
def calculate_crown_dimensions_mc(DBH,Ht,Area,ref_DBH,ref_Ht,ref_D, a_ht, b_ht, a_area, b_area, a_depth, b_depth):
    # Gapfill record with local allometry
    # Heights
    mask = np.isnan(Ht)
    Ht[mask] = random_sample_from_powerlaw_prediction_interval(DBH[mask],DBH,Ht,b_ht,a_ht,array=True)
    #Ht[mask] = CF_ht*a_ht*DBH[mask]**b_ht
    #Crown areas
    mask = np.isnan(Area)
    Area[mask] = random_sample_from_powerlaw_prediction_interval(DBH[mask],DBH,Area,b_area,a_area,array=True)

    # Apply canopy depth model (from regional database)
    #Depth = random_sample_from_powerlaw_prediction_interval(DBH,ref_DBH,ref_D,b_depth,a_depth,array=True)
    Depth = random_sample_from_powerlaw_prediction_interval(Ht,ref_Ht,ref_D,b_depth,a_depth,array=True)
    
    # Remove any existing nodata values (brought forwards from input data
    mask = np.all((~np.isnan(Depth),~np.isnan(Ht),~np.isnan(Area)),axis=0)
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

    #CanopyV = CanopyV1+CanopyV2+CanopyV3+CanopyV4

    # sanity check
    TestV = np.nansum(4*pi*a*b*c/3)
    print CanopyV.sum(),TestV
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
        print "Issue - sanity check fail: ", CanopyV.sum(),TestV
    LAD = CanopyV*leafA_per_unitV/plot_area
    return LAD, CanopyV

def calculate_LAD_profiles_generic_mc(canopy_layers, Area, D, Ht, beta_min, beta_max, plot_area, leafA_per_unitV=1.):
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
        print "Issue - sanity check fail: ", CanopyV.sum(),TestV
    LAD = CanopyV*leafA_per_unitV/plot_area
    return LAD, CanopyV



#---------------------------------------
# SMALL SAFE PLOTS DATA
# some code to load in the survey data on small stems from the small SAFE plots.
# Numbers of stems are given as a  per-metre density
def load_SAFE_small_plot_data(filename, plot_area=25.*25.):
    datatype = {'names': ('Block', 'Plot', 'TreeID', 'DBH', 'Height', 'CrownRadius'), 'formats': ('S3','i8','S8','f16','f16','f16')}
    data = np.genfromtxt(filename, skiprows = 1, delimiter = ',',dtype=datatype)
    
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
                    #rad = data['CrownRadius'][mask2][tt]
                    #if np.isfinite(rad):
                    #    sum_area[ii,pp]+=np.pi*rad**2
                    #    n_area[ii,pp]+=1.

        # get plot means
        mean_height = sum_height/n_height
        #mean_area = sum_area/n_area

        # plot
        crown_geometry = {}
        crown_geometry['dbh'] = DBH[1:]
        crown_geometry['stem_density'] = np.mean(n_stems,axis=1)[1:]/plot_area
        crown_geometry['height'] = np.mean(mean_height,axis=1)[1:]
        #crown_geometry['area'] = np.mean(mean_area,axis=1)[1:]
        block_dict[blocks[bb]]= crown_geometry
        
    return block_dict

# some code to get crown dimensions for stem size distributions 
def calculate_crown_dimensions_small_plots(DBH,Ht,stem_density,a_ht, b_ht, CF_ht, a_area, b_area, CF_area, a_depth, b_depth, CF_depth):

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
def calculate_LAD_profiles_from_stem_size_distributions(canopy_layers, Area, D, Ht, stem_density, beta, leafA_per_unitV=1.):

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
        CanopyV[i]+= np.sum( pi*(r_max[mask]/D[mask]**beta)**2/(2*beta+1) * (d2[mask]**(2*beta+1) - d1[mask]**(2*beta+1))*stem_density[mask] )

    LAD = CanopyV*leafA_per_unitV
    return LAD


#=====================================================================================================================
# Load in data for SAFE detailed subplot census - all trees >2 cm  - only interested in a subset of the fields
def load_SAFE_small_stem_census(filename, sp_area=20.**2, N_subplots = 25):

    datatype = {'names': ('Plot', 'Subplot', 'Date', 'Obs','tag', 'DBH_ns', 'DBH_ew', 'H_POM','Height', 'Flag', 'Notes'), 'formats': ('S8','i8','S10','S64','i16','f32','f32','f32','f32','S4','S32')}
    data = np.genfromtxt(filename, skiprows = 1, usecols=np.arange(0,11), delimiter = ',',dtype=datatype)
              
    data['DBH_ns']/=10. # convert from mm to cm
    data['DBH_ew']/=10. # convert from mm to cm
    
    # loop through data & remove lianas
    N = data['Plot'].size

    mask = np.ones(N,dtype='bool')
    for i in range(0,N):
        if 'liana' in data['Notes'][i]:
            mask[i] = False
        elif 'Liana' in data['Notes'][i]:
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
                    ii = np.floor(dbh/bin_width)
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

    datatype = {'names': ('Plot', 'Subplot', 'Date', 'x', 'y','tag','spp', 'Nstem', 'DBH1', 'Codes', 'Notes', 'DBH2', 'DBH3', 'DBH4'), 'formats': ('S6','i8','S10','i8','i8','S6','S6','i8','f32','S8','S32','f32','f32','f32')}
    data = np.genfromtxt(filename, skiprows = 1, delimiter = ',',dtype=datatype)
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
                ii = np.floor(dbh/bin_width)
                n_stems[ii,subplots[ss]-1]+=1.

        # loop through 2nd stems
        m2 = sp_data['DBH2']>0
        n_trees2 = np.sum(m2)
        for tt in range(n_trees2):
            dbh = sp_data['DBH2'][m2][tt]
            if dbh<10.:
                ii = np.floor(dbh/bin_width)
                n_stems[ii,subplots[ss]-1]+=1.
        
        # loop through 3rd stems
        m3 = sp_data['DBH3']>0
        n_trees3 = np.sum(m3)
        for tt in range(n_trees3):
            dbh = sp_data['DBH3'][m3][tt]
            if dbh<10.:
                ii = np.floor(dbh/bin_width)
                n_stems[ii,subplots[ss]-1]+=1.

        # loop through 4th stems
        m4 = sp_data['DBH4']>0
        n_trees4 = np.sum(m4)
        for tt in range(n_trees4):
            dbh = sp_data['DBH4'][m4][tt]
            if dbh<10.:
                ii = np.floor(dbh/bin_width)
                n_stems[ii,subplots[ss]-1]+=1.
        
        
    stem_dict['dbh'] = DBH[1:]
    stem_dict['stem_density'] = n_stems[1:,:]/sp_area
    return stem_dict


# calculate crown geometries based on stem distributions only
def calculate_crown_dimensions_for_stem_distributions(DBH,stem_density,a_ht, b_ht, CF_ht, a_area, b_area, CF_area, a_depth, b_depth, CF_depth):

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
#   where r = sqrt((x-x0)
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
# - a 3D matrix containing the calculated crown for this tree, located within the plot domain
def generate_3D_crown(canopy_matrix,x,y,z,x0,y0,Z0,Zmax,Rmax,beta):
    # generate masks for tree crown
    z_ = z-Z0
    xm,ym,z_ = np.meshgrid(x,y,z_)
    r = np.sqrt((xm-x0)**2+(ym-y0)**2)
    con1 = np.all((z<0,z>Zmax),axis=0)
    con2 = r > (z_**beta) * Rmax/(Zmax**beta)

    crown = np.ones(canopy_matrix.shape)
    crown[con1] = 0
    crown[con2] = 0
    
    return crown

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
def generate_3D_canopy(x,y,z,x0,y0,Z0,Zmax,Rmax,beta,plot_mask=np.empty([]), buffer=0):

    # first create buffer
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    x_buff = np.arange(x[0]-buffer,x[-1]+buffer,dx)
    y_buff = np.arange(y[0]-buffer,y[-1]+buffer,dy)

    # now create buffered canopy matrix
    canopy_matrix = np.zeros((x_buff.size,y_buff.size,z.size),dtype='float')

    # impose plot mask
    plot_mask_buff = np.zeros(x_buff.size,y_buff.size)
    plot_mask_buff[int(buffer/dx):x.size,int(buffer/dy):y.size]=[plot_mask.copy()]
    if plot_mask.size > 0:
        canopy_matrix[plot_mask_buff!=1]=np.nan

    # Now loop through the trees. For each tree, calculate the crown volume,
    # then add to the crown map.
    n_trees = x0.size
    for tt in range(0,n_trees):
        canopy_matrix += generate_3D_crown(canopy_matrix,x_buff,y_buff,z,x0[tt],y0[tt],Z0[tt],Zmax[tt],Rmax[tt],beta[tt])
    
    # reset PAD in overlapping areas so that PAD = 1. m2 m-3 throughout crowns
    canopy_matrix[canopy_matrix>1]=1

    return canopy_matrix
    
