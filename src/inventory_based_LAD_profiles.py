# This contains code to estimate LAD profiles based on field measurements of tree height and crown dimensions, in addition to crown depth based on a regional allometric equation.

import numpy as np
from scipy import stats


# This function reads in the crown allometry data from the database: Falster et al,. 2015; BAAD: a Biomass And Allometry Database for woody plants. Ecology, 96: 1445. doi: 10.1890/14-1889.1
def retrieve_crown_allometry(filename):
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
        CanopyV[i]+= np.sum( pi*(r_max[mask]/D[mask]**beta)**2/(2*beta+1) * (d2[mask]**(2*beta+1) - d1[mask]**(2*beta+1)) )
    
    # sanity check
    TestV = np.nansum(pi*D*r_max**2/(2*beta+1))
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
