# This contains code to estimate LAD profiles based on field measurements of tree height and crown dimensions, in addition to crown depth based on a regional allometric equation.

import numpy as np
from scipy import stats


# This function reads in the crown allometry data from the database: Falster et al,. 2015; BAAD: a Biomass And Allometry Database for woody plants. Ecology, 96: 1445. doi: 10.1890/14-1889.1
def retrieve_crown_allometry_file(filename):
    datatype = {'names': ('ID', 'Ref', 'Location', 'Lat', 'Long', 'Species', 'Family','Diameter','Height','CrownArea','CrownDepth'), 'formats': ('i8','S32','S256','f16','f16','S32','S32','f16','f16','f16','f16')}
    data = np.genfromtxt(filename, skiprows = 1, delimiter = ',',dtype=datatype)
    
    mask = np.all((~np.isnan(data['Height']),~np.isnan(data['CrownDepth'])),axis=0)
    H = data['Height'][mask]
    D = data['CrownDepth'][mask]
    logH = np.log(H)
    logD = np.log(D)
    
    # regression to find power law exponents D = a.H^b
    b, loga, r, p, serr = stats.linregress(logH,logD)
    
    model_logD = b*logH + loga
    error = logD-model_logD
    MSE = np.mean(error**2)
    CF = np.exp(MSE/2) # Correction factor due to fitting regression in log-space (Baskerville, 1972)
    a = np.exp(loga)

    return a, b, CF, r**2, p, H, D

# Derive local allometric relationship between DBH and height -> fill gaps in census data
def load_crown_survey_data(census_file):
    datatype = {'names': ('plot','subplot','date','observers','tag','DBH','H_DBH','Height','flag','alive','C1','C2','subplotX','subplotY','density','spp','cmap_date','Xfield','Yfield','Zfield','DBH_field','Height_field','CrownArea','C3','dead_flag1','dead_flag2','dead_flag3'), 'formats': ('S16','i8','S10','S32','i8','f8','f8','f8','S8','i8','S132','S132','f8','f8','f8','S64','S10','f8','f8','f8','f8','f8','f8','S132','i8','i8','i8')}
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
    mask = np.all((~np.isnan(data['CrownArea']),~np.isnan(data['DBH_field'])),axis=0)
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
    CanopyV1 = np.zeros(N_layers)
    CanopyV2 = np.zeros(N_layers)
    CanopyV3 = np.zeros(N_layers)
    CanopyV4 = np.zeros(N_layers)
    pi=np.pi
    # Formula for volume of ellipsoidal cap: V = pi*a*b*x**2*(3c-x)/c**2 where x is the vertical distance from the top of the sphere along axis c.
    # Formula for volume of ellipsoid: V = 4/3*pi*a*b*c
    for i in range(0,N_layers):
        ht_u = canopy_layers[i]
        ht_l = ht_u-layer_thickness
        # Case 1: z0+c >= ht_upper && z0-c < ht_upper && z0-c >= ht_lower
        mask1 = np.all((z0+c>=ht_u,z0-c<ht_u,z0-c>=ht_l),axis=0)
        x = z0+c-ht_u
        CanopyV1[i] += np.sum(pi/3.*a[mask1]*b[mask1]*((4*c[mask1]) - (x[mask1]**2.*(3.*c[mask1]-x[mask1])/c[mask1]**2.)))

        # Case 2: z0+c >= ht_upper && z0 - c < ht_lower
        mask2 = np.all((z0+c>=ht_u,z0-c<ht_l),axis=0)
        x1 = z0+c-ht_u
        x2 = z0+c-ht_l
        CanopyV2[i] += np.sum(pi/3.*a[mask2]*b[mask2]/c[mask2]**2*(x2[mask2]**2.*(3.*c[mask2]-x2[mask2]) - x1[mask2]**2.*(3.*c[mask2]-x1[mask2])))

        # Case 3: z0+c < ht_upper, z0+c >= ht_lower, z0-c < ht_lower
        mask3 = np.all((z0+c<ht_u,z0+c>=ht_l,z0-c<ht_l),axis=0)
        x = z0+c-ht_l
        CanopyV3[i] += np.sum(pi/3.*a[mask3]*b[mask3]*(x[mask3]**2.*(3.*c[mask3]-x[mask3])/c[mask3]**2.)) 

        # Case 4 z0+c < ht_upper, z0-c >= ht_lower
        mask4 = np.all((z0+c<ht_u,z0-c>=ht_l),axis=0)
        CanopyV4[i] += np.sum(4*pi*a[mask4]*b[mask4]*c[mask4]/3)

    CanopyV = CanopyV1+CanopyV2+CanopyV3+CanopyV4

    # sanity check
    TestV = np.nansum(4*pi*a*b*c/3)
    print CanopyV.sum(),TestV
    LAD = CanopyV*leafA_per_unitV/plot_area
    return LAD, CanopyV

