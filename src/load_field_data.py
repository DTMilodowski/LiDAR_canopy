import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

def load_leaf_chemistry(chem_file):
    chem_dtype={'names':('code','forest_type','batch','labcode','N_perc','15N','C_perc','13C','CN','quality_check'),'formats':('S32','S8','i8','S8','f16','f16','f16','f16','f16','S8')}
    chem_data = np.genfromtxt(chem_file,skip_header=1,delimiter=',',dtype=chem_dtype)

    keep_indices1 = chem_data['quality_check']!='outlier'
    keep_indices2 = np.isfinite(chem_data['N_perc'])
    keep_indices = keep_indices1*keep_indices2

    branch = chem_data['code'][keep_indices]
    N_perc = chem_data['N_perc'][keep_indices]
    C_perc = chem_data['C_perc'][keep_indices]
    CN = chem_data['CN'][keep_indices]

    forest_type = chem_data['forest_type'][keep_indices]
    return forest_type, branch, N_perc, C_perc, CN

def load_species_list(spp_file):
    spp_dtype = {'names':('forest_type','code','plot','type','ID','branch','sampling_date','spp_code','remarks','spp','data_quality'),'formats': ('S8','S32','S8','S2','i8','S8','S10','S8','S128','S128','f8')}
    spp_data = np.genfromtxt(spp_file,skip_header=1,delimiter=',',dtype=spp_dtype)

    keep = spp_data['data_quality']==1
    branch = spp_data['code'][keep]
    spp = spp_data['spp'][keep]

    # get genus from spp
    N = spp.size
    G=[]
    for i in range(N):
        G.append(spp[i].split(' ')[0])
    genus = np.asarray(G)
    return branch,spp, genus
# clean LiCOR photosynthesis data
# Version 0 (default) is my own cleaning
# Version 1 is a very similar cleaning scheme from Sabine
# Version 2 is a conservative cleaning scheme based on liCOR instructions
def clean_photosynthesis_data(Asat,Amax,Rd,version = 0):
    # 1) Rd if photosynthesis is positive when insolation(PARi) is zero, delete
    Rd=Rd[Rd['Photo']<=0]
    # 2) Filter out bad photosynthetic rates for Amax and Asat
    Amax=Amax[Amax['Photo']>=0]
    Asat=Asat[Asat['Photo']>=0]
    # 3) Stable parameters threshold (0.7)
    Rd=Rd[Rd['StableF']>0.7]
    Amax=Amax[Amax['StableF']>0.7]
    Asat=Asat[Asat['StableF']>0.7]

    # 4) Filter bad Ci
    if version == 0:
        Asat = Asat[Asat['Ci']<=Asat['CO2R']]
        Amax = Amax[Amax['Ci']<=Amax['CO2R']]

    elif np.any((version == 1,version ==2)):
        Asat=Asat[Asat['Ci']<=300]
        Asat=Asat[Asat['Ci']>=150]

    Asat = Asat[Asat['Ci']>0]
    Amax = Amax[Amax['Ci']>0]
    Rd = Rd[Rd['Ci']>0]

    if version == 2:
        # 5) Conductance threshold for Amax and Asat
        Amax=Amax[Amax['Cond']>=0.04]
        Asat=Asat[Asat['Cond']>=0.04]

    if version == 0:
        # 7) Check H2O - should be between 15 and 35
        Rd=Rd[Rd['H2OR']>=15]
        Amax=Amax[Amax['H2OR']>=15]
        Asat=Asat[Asat['H2OR']>=15]

        Rd=Rd[Rd['H2OR']<=35]
        Amax=Amax[Amax['H2OR']<=35]
        Asat=Asat[Asat['H2OR']<=35]
    return Asat, Amax, Rd

def load_photosynthesis_data(photo_file,version=0):
    datatype = {'names': ('Licor', 'code', 'Obs', 'Time', 'FTime', 'EBal', 'Photo', 'Cond', 'Ci', 'Trmmol', 'VpdL', 'CTleaf','Area', 'BLC_1', 'StmRat','BLCond','Tair','Tleaf','Tblk','CO2R','CO2S','H2OR','H2OS','RH_R','RH_S','Flow','PARi','PARo','Press','CsMch','HsMch','CsMchSD','HsMchSD','CrMchSD','HrMchSD','StableF','BLCslope','BLCoffst','f_parin','f_parout','alphaK'), 'formats': ('i8','S32','i8','S32','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16')}
    data = np.genfromtxt(photo_file, skip_header = 1, delimiter = ',',dtype=datatype)

    # Do some calculations in preparation for calculating Vcmax and Jmax
    R = 8.314 # Universal gas constant
    Oi = 210 # mmol/mol

    # The data table is laid out as continuous data with no leaf ID other than the observation number ('Obs') which iterates through from 1 to ~30 and represents the number of repetitions.  I'm going to assume that these come in sets of three: Asat, Amax, Rd.  If the end of the dataset does not yield three groups, then kick out an error.
    count = 1
    leaf_count = 0
    leaf_measurements=[1,]
    while count<len(data):
        if data['code'][count]==data['code'][count-1]:
            leaf_measurements[leaf_count]+=1
        else:
            # new leaf
            leaf_count+=1
            leaf_measurements.append(1)

        count+=1
    leaf_count+=1

    # Next up, get Asat and Rd, the Michaelis-Menten constant and CO2 compensation point for the leaf
    Asat=np.zeros(leaf_count)
    Amax=np.zeros(leaf_count)
    Rd=np.zeros(leaf_count)
    Km=np.zeros(leaf_count)
    c2co=np.zeros(leaf_count)
    Jmax=np.zeros(leaf_count)
    Vcmax=np.zeros(leaf_count)
    Ci=np.zeros(leaf_count)
    VPD=np.zeros(leaf_count)
    LeafT=np.zeros(leaf_count)
    gs=np.zeros(leaf_count)
    leaf_ID = []
    BALIplot = []
    Tree_tag = np.zeros(leaf_count)
    leaf_height = np.zeros(leaf_count)
    tree_height = np.zeros(leaf_count)
    leaf_thickness = np.zeros(leaf_count)
    leaf_area = np.zeros(leaf_count)
    #SLA = np.zeros(leaf_count)
    shade_tag = np.zeros(leaf_count)
    branch_tag = []
    forest_type = []

    start_point = 0
    err_count=0
    for i in range(0,leaf_count):

        i1 = start_point+leaf_measurements[i]
        Leaf_data = data[start_point:i1]
        leaf_ID.append(Leaf_data['code'][0])

        # 1) First, split up Asat, Amax and Rd measurements
        Rd_data = Leaf_data[Leaf_data['PARi']==0]
        Asat_data = Leaf_data[np.all((Leaf_data['CO2R']<500,Leaf_data['PARi']>=1900),axis=0)]
        Amax_data = Leaf_data[np.all((Leaf_data['CO2R']>1900,Leaf_data['PARi']>=1900,Leaf_data['CO2R']<2100),axis=0)]

        # post processing
        Asat_data, Amax_data, Rd_data = clean_photosynthesis_data(Asat_data,Amax_data,Rd_data,version)
        """
        # 2) Rd if photosynthesis is positive when insolation(PARi) is zero, delete
        Rd_data=Rd_data[Rd_data['Photo']<=0]
        # 3) Filter out bad photosynthetic rates for Amax and Asat
        Amax_data=Amax_data[Amax_data['Photo']>=0]
        Asat_data=Asat_data[Asat_data['Photo']>=0]
        # 4) Stable parameters threshold (0.7)
        Rd_data=Rd_data[Rd_data['StableF']>0.7]
        Amax_data=Amax_data[Amax_data['StableF']>0.7]
        Asat_data=Asat_data[Asat_data['StableF']>0.7]
        """

        """
        # 5) Conductance threshold for Amax and Asat - v. conservative
        Amax_data=Amax_data[Amax_data['Cond']<=0.04]
        Asat_data=Asat_data[Asat_data['Cond']<=0.04]

        # 6) Filter bad Asat Ci
        Asat_data=Asat_data[Asat_data['Ci']<=300]
        Asat_data=Asat_data[Asat_data['Ci']>=150]
        """

        """
        Asat_data = Asat_data[Asat_data['Ci']<=Asat_data['CO2R']]
        Amax_data = Amax_data[Amax_data['Ci']<=Amax_data['CO2R']]

        Asat_data = Asat_data[Asat_data['Ci']>0]
        Amax_data = Amax_data[Amax_data['Ci']>0]
        Rd_data = Rd_data[Rd_data['Ci']>0]

        # 7) Check H2O - should be between 15 and 35
        Rd_data=Rd_data[Rd_data['H2OR']>=15]
        Amax_data=Amax_data[Amax_data['H2OR']>=15]
        Asat_data=Asat_data[Asat_data['H2OR']>=15]

        Rd_data=Rd_data[Rd_data['H2OR']<=35]
        Amax_data=Amax_data[Amax_data['H2OR']<=35]
        Asat_data=Asat_data[Asat_data['H2OR']<=35]
        """
        #---------------------------------------------------------------------
        # Now post-processing is complete calculate c2co and Km
        # following Crous et al., 2013
        # 1) carbon compensation point
        c2co_Asat=38.892*np.exp(20437*(Asat_data['Tleaf']+273.15-298.15)/(298.15*R*((Asat_data['Tleaf'])+273.15)))
        # 2) Michaelis constant for CO2
        Kc_Asat=404.9*np.exp(79403*(Asat_data['Tleaf']+273.15-298.15)/(298.15*R*((Asat_data['Tleaf'])+273.15)))
        # 3) Michaelis constant for O2
        Ko_Asat=278.4*np.exp(36380*(Asat_data['Tleaf']+273.15-298.15)/(298.15*R*((Asat_data['Tleaf'])+273.15)))
        # 4) Michaelis-Menten constant
        Km_Asat=Kc_Asat*(1+Oi/Ko_Asat)

        # Now get the average values for the leaf
        if np.sum(np.isfinite(Asat_data['Photo']))>0 and np.sum(np.isfinite(Rd_data['Photo']))>0 and np.sum(np.isfinite(Amax_data['Photo']))>0:
            Asat[i]=np.mean(Asat_data['Photo'][np.isfinite(Asat_data['Photo'])])
            Amax[i]=np.mean(Amax_data['Photo'][np.isfinite(Amax_data['Photo'])])
            Rd[i]=-np.mean(Rd_data['Photo'][np.isfinite(Rd_data['Photo'])])
            Km[i]=np.mean(Km_Asat[np.isfinite(Km_Asat)])
            c2co[i]=np.mean(c2co_Asat[np.isfinite(c2co_Asat)])
            Ci[i]=np.mean(Asat_data['Ci'][np.isfinite(Asat_data['Ci'])])
            Vcmax[i] = (Asat[i]+Rd[i])*((Ci[i]+Km[i])/(Ci[i]-c2co[i])) # De Kauwe et al. 2016
            Jmax[i] = 1.197*0.847*Vcmax[i] # Walker et al. 2014
            LeafT[i] = np.mean(Asat_data['Tleaf'][np.isfinite(Asat_data['Tleaf'])])
            gs[i] = np.mean(Asat_data['Cond'][np.isfinite(Asat_data['Cond'])])
            VPD[i] = np.mean(Asat_data['VpdL'][np.isfinite(Asat_data['VpdL'])])
            if np.isnan(Asat[i]):
                print("error")
                err_count+=1
                print(start_point, i1, len(Asat_data), err_count)
        else:
            Asat[i]=np.nan
            Amax[i]=np.nan
            Rd[i]=np.nan
            Km[i]=np.nan
            c2co[i]=np.nan
            Ci[i]=np.nan
            Vcmax[i] = np.nan
            Jmax[i] = np.nan
            LeafT[i] = np.nan
            gs[i] = np.nan
            VPD[i] = np.nan

        # now get plot, tree, branch information
        plot = Leaf_data['code'][0].split('-')[0]
        tree = Leaf_data['code'][0].split('-')[1]
        branch= Leaf_data['code'][0].split('-')[2]
        if branch=='B1S':
            branch='BS'
        elif branch=='B1SH':
            branch='BSH'
        branch_tag.append(plot+'-'+tree+'-'+branch)
        BALIplot.append(plot)
        if branch[-1] == 'H':
            shade_tag[i] = 1
        elif branch[-1] == 'S':
            shade_tag[i] = 0
        else:
            shade_tag[i] = np.nan

        # update index
        start_point=i1

    plot = np.asarray(BALIplot,dtype='S8')
    leafID = np.asarray(leaf_ID,dtype='S16')
    branchID = np.asarray(branch_tag)
    Ftype = np.asarray(forest_type)
    return plot, leafID, branchID, shade_tag, Ci, VPD, LeafT, gs, Km, c2co, Asat, Amax, Rd, Vcmax, Jmax

def load_branch_info(branch_file):
    branch_dtype = {'names': ('date', 'forest_type', 'code', 'plot', 'tree', 'branch', 'tree_height', 'branch_height', 'time_cut'), 'formats': ('S10','S8','S32','S8','S8','S8','f16','f16','f16')}
    branch_data = np.genfromtxt(branch_file, skip_header = 1, delimiter = ',',dtype=branch_dtype)

    forest_type = branch_data['forest_type']
    branch_ID = branch_data['code']
    tree_height = branch_data['tree_height']
    branch_height = branch_data['branch_height']
    N = forest_type.size
    branch_tag = []
    for i in range(0,N):
        # now get plot, tree, branch information
        plot = branch_data['code'][i].split('-')[0]
        tree = branch_data['code'][i].split('-')[1]
        branch = branch_data['code'][i].split('-')[2]
        if branch=='B1S':
            branch='BS'
        elif branch=='B1SH':
            branch='BSH'
        branch_tag.append(plot+'-'+tree+'-'+branch)

    branch_ID=np.asarray(branch_tag)
    return forest_type, branch_ID, tree_height, branch_height

def load_leaf_traits(leaf_file):
    leaf_dtype = {'names': ('date', 'forest_type', 'code', 'plot', 'tree', 'branch', 'leaf' , 'leaf_thickness', 'fresh_wt', 'dry_wt_g', 'dry_wt_mg', 'leafA_cm2', 'leafA_mm2', 'SLA', 'LDMC','q_check'), 'formats': ('S10','S8','S32','S8','S8','S8','S8','f16','f16','f16','f16','f16','f16','f16','f16','S8')}
    leaf_data = np.genfromtxt(leaf_file, skip_header = 1, delimiter = ',',dtype=leaf_dtype)
    leaf_data = leaf_data[leaf_data['q_check']=='ok']

    forest_type = leaf_data['forest_type']
    leaf_ID = leaf_data['code']

    N = leaf_ID.size
    shade_tag = np.zeros(N)
    branch_tag = []
    for i in range(0,N):
        # now get plot, tree, branch information
        plot = leaf_data['code'][i].split('-')[0]
        tree = leaf_data['code'][i].split('-')[1]
        branch = leaf_data['code'][i].split('-')[2]
        if branch=='B1S':
            branch='BS'
        elif branch=='B1SH':
            branch='BSH'
        branch_tag.append(plot+'-'+tree+'-'+branch)
        if branch[-1] == 'H':
            shade_tag[i] = 1
        elif branch[-1] == 'S':
            shade_tag[i] = 0
        else:
            shade_tag[i] = np.nan

    leaf_thickness = leaf_data['leaf_thickness']
    leaf_area = leaf_data['leafA_mm2']
    SLA = leaf_data['SLA'] # SLA in mm^2/mg
    # Convert SLA to LMA
#    LMA = 1/(SLA*10**3*10**-2) #LMA in g/cm^2
    ###Check!!!

    LMA = (1./SLA)*10.**3 #LMA in g/m^2

    branch_ID = np.asarray(branch_tag)
    return forest_type, leaf_ID, branch_ID, shade_tag, leaf_thickness, leaf_area, SLA, LMA



def load_leaf_traits_old(traits_file, branch_file, leaf_file):
    datatype = {'names': ('Licor', 'code', 'Obs', 'Time', 'FTime', 'EBal', 'Photo', 'Cond', 'Ci', 'Trmmol', 'VpdL', 'CTleaf','Area', 'BLC_1', 'StmRat','BLCond','Tair','Tleaf','Tblk','CO2R','CO2S','H2OR','H2OS','RH_R','RH_S','Flow','PARi','PARo','Press','CsMch','HsMch','CsMchSD','HsMchSD','CrMchSD','HrMchSD','StableF','BLCslope','BLCoffst','f_parin','f_parout','alphaK'), 'formats': ('i8','S32','i8','S32','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16','f16')}
    data = np.genfromtxt(traits_file, skip_header = 1, delimiter = ',',dtype=datatype)

    branch_dtype = {'names': ('date', 'forest_type', 'code', 'plot', 'tree', 'branch', 'tree_height', 'branch_height', 'time_cut'), 'formats': ('S10','S8','S32','S8','S8','S8','f16','f16','f16')}
    branch_data = np.genfromtxt(branch_file, skip_header = 1, delimiter = ',',dtype=branch_dtype)

    leaf_dtype = {'names': ('date', 'forest_type', 'code', 'plot', 'tree', 'branch', 'leaf' , 'leaf_thickness', 'fresh_wt', 'dry_wt_g', 'dry_wt_mg', 'leafA_cm2', 'leafA_mm2', 'SLA', 'LDMC','q_check'), 'formats': ('S10','S8','S32','S8','S8','S8','S8','f16','f16','f16','f16','f16','f16','f16','f16','S8')}
    leaf_data = np.genfromtxt(leaf_file, skip_header = 1, delimiter = ',',dtype=leaf_dtype)
    leaf_data = leaf_data[leaf_data['q_check']=='ok']

    # Do some calculations in preparation for calculating Vcmax and Jmax
    R = 8.314 # Universal gas constant
    Oi = 210 # mmol/mol


    # The data table is laid out as continuous data with no leaf ID other than the observation number ('Obs') which iterates through from 1 to ~30 and represents the number of repetitions.  I'm going to assume that these come in sets of three: Asat, Amax, Rd.  If the end of the dataset does not yield three groups, then kick out an error.
    count = 1
    leaf_count = 0
    leaf_measurements=[1,]
    while count<len(data):
        if data['code'][count]==data['code'][count-1]:
            leaf_measurements[leaf_count]+=1
        else:
            # new leaf
            leaf_count+=1
            leaf_measurements.append(1)

        count+=1
    leaf_count+=1

    # Next up, get Asat and Rd, the Michaelis-Menten constant and CO2 compensation point for the leaf
    Asat=np.zeros(leaf_count)
    Amax=np.zeros(leaf_count)
    Rd=np.zeros(leaf_count)
    Km=np.zeros(leaf_count)
    c2co=np.zeros(leaf_count)
    Jmax=np.zeros(leaf_count)
    Vcmax=np.zeros(leaf_count)
    Ci=np.zeros(leaf_count)
    VPD=np.zeros(leaf_count)
    LeafT=np.zeros(leaf_count)
    gs=np.zeros(leaf_count)
    leaf_ID = []
    BALIplot = []
    Tree_tag = np.zeros(leaf_count)
    leaf_height = np.zeros(leaf_count)
    tree_height = np.zeros(leaf_count)
    leaf_thickness = np.zeros(leaf_count)
    leaf_area = np.zeros(leaf_count)
    SLA = np.zeros(leaf_count)
    shade_tag = np.zeros(leaf_count)
    branch_tag = []
    forest_type = []

    start_point = 0
    err_count=0
    for i in range(0,leaf_count):

        i1 = start_point+leaf_measurements[i]
        Leaf_data = data[start_point:i1]
        leaf_ID.append(Leaf_data['code'][0])
        #Tree_tag.append(Leaf_data['code'][0].split("-")[1])
        #Tree_tag[i] = Leaf_data['tree'][0]

        # 1) First, split up Asat, Amax and Rd measurements
        Rd_data = Leaf_data[Leaf_data['PARi']==0]
        temp_data = Leaf_data[Leaf_data['PARi']>0]

        Asat_data = temp_data[temp_data['CO2R']<500]
        Amax_data = temp_data[temp_data['CO2R']>500]

        temp_data = None

        # post processing
        # 2) Rd if photosynthesis is positive when insolation(PARi) is zero, delete
        Rd_data=Rd_data[Rd_data['Photo']<=0]
        # 3) Filter out bad photosynthetic rates for Amax and Asat
        Amax_data=Amax_data[Amax_data['Photo']>=0]
        Asat_data=Asat_data[Asat_data['Photo']>=0]
        # 4) Stable parameters threshold (0.7)
        Rd_data=Rd_data[Rd_data['StableF']>0.7]
        Amax_data=Amax_data[Amax_data['StableF']>0.7]
        Asat_data=Asat_data[Asat_data['StableF']>0.7]
        """
        # 5) Conductance threshold for Amax and Asat
        Amax_data=Amax_data[Amax_data['Cond']<=0.04]
        Asat_data=Asat_data[Asat_data['Cond']<=0.04]

        # 6) Filter bad Asat Ci
        Asat_data=Asat_data[Asat_data['Ci']<=300]
        Asat_data=Asat_data[Asat_data['Ci']>=150]
        """
        Asat_data = Asat_data[Asat_data['Ci']<=Asat_data['CO2R']]
        Amax_data = Amax_data[Amax_data['Ci']<=Amax_data['CO2R']]

        Asat_data = Asat_data[Asat_data['Ci']>0]
        Amax_data = Amax_data[Amax_data['Ci']>0]
        Rd_data = Rd_data[Rd_data['Ci']>0]

        # 7) Check H2O - should be between 15 and 35
        Rd_data=Rd_data[Rd_data['H2OR']>=15]
        Amax_data=Amax_data[Amax_data['H2OR']>=15]
        Asat_data=Asat_data[Asat_data['H2OR']>=15]

        Rd_data=Rd_data[Rd_data['H2OR']<=35]
        Amax_data=Amax_data[Amax_data['H2OR']<=35]
        Asat_data=Asat_data[Asat_data['H2OR']<=35]

        # 8) For Asat and Amax, make sure that PARi is sufficiently high
        Asat_data=Asat_data[Asat_data['PARi']>=1900]
        Amax_data=Amax_data[Amax_data['PARi']>=1900]

        #---------------------------------------------------------------------
        # Now post-processing is complete calculate c2co and Km
        # following Crous et al., 2013
        # 1) carbon compensation point
        c2co_Asat=38.892*np.exp(20437*(Asat_data['Tleaf']+273.15-298.15)/(298.15*R*((Asat_data['Tleaf'])+273.15)))
        # 2) Michaelis constant for CO2
        Kc_Asat=404.9*np.exp(79403*(Asat_data['Tleaf']+273.15-298.15)/(298.15*R*((Asat_data['Tleaf'])+273.15)))
        # 3) Michaelis constant for O2
        Ko_Asat=278.4*np.exp(36380*(Asat_data['Tleaf']+273.15-298.15)/(298.15*R*((Asat_data['Tleaf'])+273.15)))
        # 4) Michaelis-Menten constant
        Km_Asat=Kc_Asat*(1+Oi/Ko_Asat)

        # Now get the average values for the leaf
        if np.sum(np.isfinite(Asat_data['Photo']))>0 and np.sum(np.isfinite(Rd_data['Photo']))>0 and np.sum(np.isfinite(Amax_data['Photo']))>0:
            Asat[i]=np.mean(Asat_data['Photo'][np.isfinite(Asat_data['Photo'])])
            Amax[i]=np.mean(Amax_data['Photo'][np.isfinite(Amax_data['Photo'])])
            Rd[i]=-np.mean(Rd_data['Photo'][np.isfinite(Rd_data['Photo'])])
            Km[i]=np.mean(Km_Asat[np.isfinite(Km_Asat)])
            c2co[i]=np.mean(c2co_Asat[np.isfinite(c2co_Asat)])
            Ci[i]=np.mean(Asat_data['Ci'][np.isfinite(Asat_data['Ci'])])
            Vcmax[i] = (Asat[i]+Rd[i])*((Ci[i]+Km[i])/(Ci[i]-c2co[i])) # De Kauwe et al. 2016
            Jmax[i] = 1.197*0.847*Vcmax[i] # Walker et al. 2014
            LeafT[i] = np.mean(Asat_data['Tleaf'][np.isfinite(Asat_data['Tleaf'])])
            gs[i] = np.mean(Asat_data['Cond'][np.isfinite(Asat_data['Cond'])])
            VPD[i] = np.mean(Asat_data['VpdL'][np.isfinite(Asat_data['VpdL'])])
            if np.isnan(Asat[i]):
                print("error")
                err_count+=1
                print(start_point, i1, len(Asat_data), err_count)
        else:
            Asat[i]=np.nan
            Amax[i]=np.nan
            Rd[i]=np.nan
            Km[i]=np.nan
            c2co[i]=np.nan
            Ci[i]=np.nan
            Vcmax[i] = np.nan
            Jmax[i] = np.nan
            LeafT[i] = np.nan
            gs[i] = np.nan
            VPD[i] = np.nan

        # now get the LMA for this leaf
        if np.sum(leaf_data['code']==leaf_ID[i])==1:
            SLA[i]=leaf_data['SLA'][leaf_data['code']==leaf_ID[i]]
            leaf_thickness[i]=leaf_data['leaf_thickness'][leaf_data['code']==leaf_ID[i]]
            leaf_area[i]=leaf_data['leafA_mm2'][leaf_data['code']==leaf_ID[i]]
            # now get leaf height
            plot = leaf_data['plot'][leaf_data['code']==leaf_ID[i]]
            branch = leaf_data['branch'][leaf_data['code']==leaf_ID[i]]
            tree = leaf_data['tree'][leaf_data['code']==leaf_ID[i]]
            Tree_tag[i] = tree[0]
            branch_ID = plot[0]+'-T'+tree[0]+'-'+branch[0]
            branch_tag.append(branch_ID)
            BALIplot.append(plot[0])
            if np.sum(branch_data['code']==branch_ID)==1:
                leaf_height[i] = branch_data['branch_height'][branch_data['code']==branch_ID]
                tree_height[i] = branch_data['tree_height'][branch_data['code']==branch_ID]
                forest_type.append(branch_data['forest_type'][branch_data['code']==branch_ID][0])
            elif np.sum(branch_data['code']==branch_ID)==0:
                leaf_height[i] = np.nan
                tree_height[i] = np.nan
                forest_type.append('_')
                print(branch_ID, "\t missing branch")
            else:
                leaf_height[i] = np.nan
                tree_height[i] = np.nan
                forest_type.append('_')
                print(branch_ID, "\t multiple branches with this ID")

            if branch[0][-1] == 'H':
                shade_tag[i] = 1
            elif branch[0][-1] == 'S':
                shade_tag[i] = 0
            else:
                shade_tag[i] = np.nan


        elif np.sum(leaf_data['code']==leaf_ID[i])==0:
            SLA[i]=np.nan
            leaf_height[i] = np.nan
            shade_tag[i] = np.nan
            Tree_tag[i] = np.nan
            BALIplot.append('_')
            branch_tag.append('_')
            forest_type.append('_')
            print(leaf_ID[i], "\t missing leaf in leaf area data")
        else:
            SLA[i]=np.nan
            leaf_height[i] = np.nan
            tree_height[i] = np.nan
            shade_tag[i] = np.nan
            BALIplot.append('_')
            branch_tag.append('_')
            forest_type.append('_')
            Tree_tag[i] = np.nan
            print(leaf_ID[i], "\t more than one leaf associated with this ID")

        # update index
        start_point=i1

    # Convert SLA to LMA
    LMA = 1/SLA*1000

    plot = np.asarray(BALIplot,dtype='S8')
    leafID = np.asarray(leaf_ID,dtype='S16')
    branchID = np.asarray(branch_tag)
    Ftype = np.asarray(forest_type)
    return plot, Ftype, Tree_tag, leafID, branchID, shade_tag, tree_height, leaf_height, leaf_thickness, leaf_area, SLA, LMA, Ci, VPD, LeafT, gs, Km, c2co, Asat, Amax, Rd, Vcmax, Jmax

# Collate traits data at the branch level
# Version refers to the procedure used to clean LiCOR photosynthesis data
# Version 0 (default) is my own cleaning scheme
# Version 1 is a very similar cleaning scheme from Sabine
# Version 2 is a conservative cleaning scheme based on liCOR instructions
def collate_branch_level_traits(chem_file,photo_file,branch_file,leaf_file,spp_file,version=0):
    plot_photo, leaf_ID_photo, branch_ID_photo, shade_tag_photo, Ci, VPD, LeafT, gs, Km, c2co, Asat, Amax, Rd, Vcmax, Jmax= load_photosynthesis_data(photo_file,version)

    forest_type_branch, branch_ID, tree_height, branch_height = load_branch_info(branch_file)

    forest_type_leaf, leaf_ID_leaf, branch_ID_leaf, shade_tag_leaf, leaf_thickness, leaf_area, SLA, LMA = load_leaf_traits(leaf_file)

    forest_type_chem, branch_ID_chem, N_perc, C_perc, CN = load_leaf_chemistry(chem_file)

    branch_ID_spp, spp, genus = load_species_list(spp_file)

    unique_branches = np.unique(np.concatenate((branch_ID,branch_ID_chem,branch_ID_photo, branch_ID_leaf,branch_ID_spp)))
    N = unique_branches.size

    branch_SLA = np.zeros(N)*np.nan
    branch_LMA = np.zeros(N)*np.nan
    branch_LeafArea = np.zeros(N)*np.nan
    branch_VPD = np.zeros(N)*np.nan
    branch_Rd = np.zeros(N)*np.nan
    branch_Vcmax = np.zeros(N)*np.nan
    branch_Jmax = np.zeros(N)*np.nan
    branch_LeafThickness = np.zeros(N)*np.nan
    branch_ShadeTag = np.zeros(N)*np.nan
    branch_N = np.zeros(N)*np.nan
    branch_Na = np.zeros(N)*np.nan
    branch_C = np.zeros(N)*np.nan
    branch_CNratio = np.zeros(N)*np.nan
    branch_LeafHeight = np.zeros(N)*np.nan
    branch_gs = np.zeros(N)*np.nan
    branch_spp = []
    branch_Ftype = []
    for i in range(0,N):

        br = unique_branches[i]
        ftype_flag = 0
        leafdata_flag = 0
        # get branch info
        if np.sum(branch_ID==br)>0:
            branch_LeafHeight[i]=np.mean(branch_height[branch_ID==br])
            branch_Ftype.append(forest_type_branch[branch_ID==br][0])
            ftype_flag = 1
        else:
            branch_Ftype.append('_')

        # get leaf info
        if np.sum(branch_ID_leaf==br)>0:
            branch_SLA[i]=np.mean(SLA[branch_ID_leaf==br])
            branch_LMA[i]=np.mean(LMA[branch_ID_leaf==br])
            branch_LeafArea[i]=np.mean(leaf_area[branch_ID_leaf==br])
            branch_LeafThickness[i]=np.mean(leaf_thickness[branch_ID_leaf==br])
            branch_ShadeTag[i]=np.mean(shade_tag_leaf[branch_ID_leaf==br])
            leafdata_flag = 1
            if ftype_flag==0:
                branch_Ftype[i]=forest_type_leaf[branch_ID_leaf==br][0]
                ftype_flag = 1

        # get photo info
        if np.sum(branch_ID_photo==br)>0:
            branch_VPD[i]=np.mean(VPD[branch_ID_photo==br])
            branch_Rd[i]=np.mean(Rd[branch_ID_photo==br])
            branch_gs[i]=np.mean(gs[branch_ID_photo==br])
            branch_Vcmax[i]=np.mean(Vcmax[branch_ID_photo==br])
            branch_Jmax[i]=np.mean(Jmax[branch_ID_photo==br])
            if np.isnan(branch_ShadeTag[i]):
                branch_ShadeTag[i]=np.mean(shade_tag_photo[branch_ID_photo==br])

        # get chem info
        if np.sum(branch_ID_chem==br)>0:
            branch_N[i]=np.mean(N_perc[branch_ID_chem==br])
            branch_C[i]=np.mean(C_perc[branch_ID_chem==br])
            branch_CNratio[i]=np.mean(CN[branch_ID_chem==br])
            if ftype_flag==0:
                branch_Ftype[i]=forest_type_chem[branch_ID_chem==br][0]
                ftype_flag = 1
            if leafdata_flag == 1:
                branch_Na[i] = branch_N[i]/100.*branch_LMA[i]

        if np.sum(branch_ID_spp==br)==1:
            branch_spp.append(spp[branch_ID_spp==br][0])
        elif np.sum(branch_ID_spp==br)>1:
            branch_spp.append(spp[branch_ID_spp==br][0])
        else:
            branch_spp.append('_')

    species = np.asarray(branch_spp)
    ForestType = np.asarray(branch_Ftype)

    return unique_branches, species, genus, branch_N, branch_Na, branch_C, branch_CNratio, branch_SLA, branch_LMA, branch_LeafArea, branch_LeafThickness, branch_LeafHeight, branch_VPD, branch_Rd, branch_Vcmax, branch_Jmax, branch_ShadeTag, branch_gs, ForestType



def read_ICP_census_data(census_file):

    datatype = {'names': ('ForestType', 'Plot', 'Subplot', 'CensusDate1', 'Observers1',
                'TagNumber1', 'D_POM1', 'H_POM1', 'Height1', 'RAINFOR_flag1', 'Alive_flag1',
                'C_stem1','C_coarse_root1', 'Comment1', 'CommentData1', 'CensusDate2',
                'Observers2', 'TagNumber2', 'D_POM2', 'H_POM2', 'Height2', 'RAINFOR_flag2',
                'Alive_flag2', 'C_stem2','C_coarse_root2', 'Comment2', 'CommentData2',
                'CensusDate3', 'Observers3', 'TagNumber3', 'D_POM3', 'H_POM3', 'Height3',
                'RAINFOR_flag3', 'Alive_flag3', 'C_stem3','C_coarse_root3', 'Comment3',
                'CommentData3', 'SubplotX', 'SubplotY', 'SpeciesID', 'WoodDensity','Source',
                'Quality'),
                'formats': ('U16','U16','int_','U10','U64','f','f','f','f',
                'U8','i8','f','f','U16','U64','U10','U64','f','f','f','f',
                'U8','i8','f','f','U16','U64','U10','U64','f','f','f','f',
                'U8','i8','f','f','U16','U64','int_','int_','U32','f','U100','U32')}
    census_data = np.genfromtxt(census_file, skip_header = 1, delimiter = ',',dtype=datatype)

    nrows = census_data['Plot'].size
    BALI_plot = []
    Subplot = np.zeros(nrows)
    CensusDates = np.zeros((nrows,3),dtype = 'datetime64[D]')
    TreeTag = np.zeros(nrows)
    AltTag = np.zeros(nrows)
    DPOM = np.zeros((nrows,3))
    HPOM = np.zeros((nrows,3))
    Height = np.zeros((nrows,3))
    C_stem = np.zeros((nrows,3))
    C_coarse_root = np.zeros((nrows,3))
    RAINFOR_flag = []
    Alive_flag = np.ones((nrows,3)) # assume alive unless specified
    Species = []
    SubplotCoords = np.ones((nrows,2))
    WoodDensity = np.zeros(nrows)
    AltTag[:] = np.nan

    for i in range(0,nrows):

        BALI_plot.append(census_data['Plot'][i])
        Subplot[i] = census_data['Subplot'][i]
        Species.append(census_data['SpeciesID'][i])
        SubplotCoords[i,0] = census_data['SubplotX'][i]
        SubplotCoords[i,1] = census_data['SubplotY'][i]
        WoodDensity[i] = census_data['WoodDensity'][i]
        RAINFOR_flag.append([])

        # Census #1
        prior_census = 0
        if np.isfinite(census_data['TagNumber1'][i]):
            prior_census = 1
            if len(census_data['CensusDate1'][i])<3:
                CensusDates[i,0] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
            elif census_data['CensusDate1'][i][2]=="/":
                day,month,year = census_data['CensusDate1'][i].split("/")
                CensusDates[i,0] = np.datetime64(year+'-'+month+'-'+day)
            elif census_data['CensusDate1'][i][2]=="-":
                month_string = np.asarray(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
                month_num = np.asarray(['01','02','03','04','05','06','07','08','09','10','11','12'])
                day,month,year = census_data['CensusDate1'][i].split("-")
                year_str = str(2000+int(year))
                CensusDates[i,0] = np.datetime64(year_str+'-'+month_num[month_string==month][0]+'-'+day)
            else:
                print("row", i, "\t incorrect date format in CensusDate1")
                CensusDates[i,0] =np.datetime64(0,'D')# np.datetime64('1990-01-01')

            TreeTag[i] = census_data['TagNumber1'][i]
            DPOM[i,0] = census_data['D_POM1'][i]
            HPOM[i,0] = census_data['H_POM1'][i]
            Height[i,0] = census_data['Height1'][i]
            C_stem[i,0] = census_data['C_stem1'][i]
            C_coarse_root[i,0] = census_data['C_coarse_root1'][i]
            RAINFOR_flag[i].append(census_data['RAINFOR_flag1'][i])
            Alive_flag[i,0] = census_data['Alive_flag1'][i]
        else:
            CensusDates[i,0] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
            TreeTag[i] = np.nan
            DPOM[i,0] = np.nan
            HPOM[i,0] = np.nan
            Height[i,0] = np.nan
            C_stem[i,0] = 0#np.nan
            C_coarse_root[i,0] = 0#np.nan
            RAINFOR_flag[i].append('NaN')
            Alive_flag[i,0] = np.nan

        # Subsequent censuses
        for j in range(1,3):
            if prior_census==0: # New recruits
                if np.isfinite(census_data['TagNumber'+str(j+1)][i]):
                    prior_census = 1
                    if len(census_data['CensusDate'+str(j+1)][i])<3:
                        CensusDates[i,j] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
                    elif census_data['CensusDate'+str(j+1)][i][2]=="/":
                        day,month,year = census_data['CensusDate'+str(j+1)][i].split("/")
                        CensusDates[i,j] = np.datetime64(year+'-'+month+'-'+day)
                    elif census_data['CensusDate'+str(j+1)][i][2]=="-":
                        month_string = np.asarray(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
                        month_num = np.asarray(['01','02','03','04','05','06','07','08','09','10','11','12'])
                        day,month,year = census_data['CensusDate'+str(j+1)][i].split("-")
                        year_str = str(2000+int(year))
                        CensusDates[i,j] = np.datetime64(year_str+'-'+month_num[month_string==month][0]+'-'+day)
                    else:
                        print("row", i, "\t incorrect date format in CensusDate"+str(j+1))
                        CensusDates[i,j] =np.datetime64(0,'D')# np.datetime64('1990-01-01')

                    #day,month,year = census_data['CensusDate'+str(j+1)][i].split("/")
                    #CensusDates[i,j] = np.datetime64(year+'-'+month+'-'+day)
                    TreeTag[i] = census_data['TagNumber'+str(j+1)][i]
                    DPOM[i,j] = census_data['D_POM'+str(j+1)][i]
                    HPOM[i,j] = census_data['H_POM'+str(j+1)][i]
                    Height[i,j] = census_data['Height'+str(j+1)][i]
                    C_stem[i,j] = census_data['C_stem'+str(j+1)][i]
                    C_coarse_root[i,j] = census_data['C_coarse_root'+str(j+1)][i]
                    RAINFOR_flag[i].append(census_data['RAINFOR_flag'+str(j+1)][i])
                    Alive_flag[i,j] = census_data['Alive_flag'+str(j+1)][i]
                else:
                    CensusDates[i,j] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
                    TreeTag[i] = np.nan
                    DPOM[i,j] = np.nan
                    HPOM[i,j] = np.nan
                    Height[i,j] = np.nan
                    C_stem[i,j] = 0
                    C_coarse_root[i,j] = 0
                    RAINFOR_flag[i].append('NaN')
                    Alive_flag[i,j] = np.nan

            else: # trees present in prior census
                if np.isfinite(census_data['TagNumber'+str(j+1)][i]):

                    if len(census_data['CensusDate'+str(j+1)][i])<3:
                        CensusDates[i,j] = np.datetime64(0,'D')#np.datetime64('1990-01-01')
                    elif census_data['CensusDate'+str(j+1)][i][2]=="/":
                        day,month,year = census_data['CensusDate'+str(j+1)][i].split("/")
                        CensusDates[i,j] = np.datetime64(year+'-'+month+'-'+day)
                    elif census_data['CensusDate'+str(j+1)][i][2]=="-":
                        month_string = np.asarray(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
                        month_num = np.asarray(['01','02','03','04','05','06','07','08','09','10','11','12'])
                        day,month,year = census_data['CensusDate'+str(j+1)][i].split("-")
                        year_str = str(2000+int(year))
                        CensusDates[i,j] = np.datetime64(year_str+'-'+month_num[month_string==month][0]+'-'+day)
                    else:
                        print("row", i, "\t incorrect date format in CensusDate"+str(j+1))
                        CensusDates[i,j] = np.datetime64(0,'D')#np.datetime64('1990-01-01')

                    DPOM[i,j] = census_data['D_POM'+str(j+1)][i]
                    HPOM[i,j] = census_data['H_POM'+str(j+1)][i]
                    Height[i,j] = census_data['Height'+str(j+1)][i]
                    C_stem[i,j] = census_data['C_stem'+str(j+1)][i]
                    C_coarse_root[i,j] = census_data['C_coarse_root'+str(j+1)][i]
                    RAINFOR_flag[i].append(census_data['RAINFOR_flag'+str(j+1)][i])
                    Alive_flag[i,j] = census_data['Alive_flag'+str(j+1)][i]
                    if census_data['TagNumber'+str(j+1)][i]!=TreeTag[i]:
                        AltTag[i]= census_data['TagNumber'+str(j+1)][i]
                else:
                    CensusDates[i,j] =np.datetime64(0,'D')# np.datetime64('1990-01-01')
                    #TreeTag[i] = np.nan
                    DPOM[i,j] = np.nan
                    HPOM[i,j] = np.nan
                    Height[i,j] = np.nan
                    C_stem[i,j] = 0#np.nan
                    C_coarse_root[i,j] = 0#np.nan
                    RAINFOR_flag[i].append('NaN')
                    Alive_flag[i,j] = np.nan

    plot_names=np.asarray(BALI_plot)
    RAINFOR = np.asarray(RAINFOR_flag)
    spp = np.asarray(Species)
    C_stem[np.isnan(C_stem)]=0
    C_coarse_root[np.isnan(C_coarse_root)]=0
    return plot_names, Subplot, CensusDates, TreeTag, AltTag, DPOM, HPOM, Height, C_stem, C_coarse_root, RAINFOR, Alive_flag, spp, SubplotCoords, WoodDensity

def collate_tree_level_data(census_file,traits_file,branch_file,leaf_file):

    BALI_plot_traits, forest_type, TreeTag_traits, leaf_ID, branchID, shade_tag, tree_height, leaf_height, leaf_thickness, leaf_area, SLA, LMA, Ci, VPD, LeafT, gs, Km, c2co, Asat, Amax, Rd, Vcmax, Jmax = load_leaf_traits(traits_file, branch_file, leaf_file)
    BALI_plot_traits[BALI_plot_traits=='BEL']='Belian'
    BALI_plot_traits[BALI_plot_traits=='SER']='Seraya'
    BALI_plot_traits[BALI_plot_traits=='DAS1']='Danum1'
    BALI_plot_traits[BALI_plot_traits=='DAF2']='Danum2'

    BALI_plot_census, Subplot, CensusDates, TreeTag_census, AltTag, DPOM, HPOM, Height, C_stem, C_coarse_root, RAINFOR_flag, Alive_flag, Species, SubplotCoords, WoodDensity = read_ICP_census_data(census_file)

    # Loop through the trees in the census and attach the relevant traits data
    N_trees = BALI_plot_census.size
    leaf_template = np.zeros((N_trees,10))*np.nan
    Plot = np.zeros(N_trees,dtype='S8')
    shade_tag_tree = leaf_template.copy()
    leaf_height_tree = leaf_template.copy()
    canopy_position_tree = leaf_template.copy()
    canopy_depth = leaf_template.copy()
    leaf_thickness_tree = leaf_template.copy()
    SLA_tree = leaf_template.copy()
    LMA_tree = leaf_template.copy()
    Ci_tree = leaf_template.copy()
    VPD_tree = leaf_template.copy()
    LeafT_tree = leaf_template.copy()
    gs_tree = leaf_template.copy()
    Km_tree = leaf_template.copy()
    c2co_tree = leaf_template.copy()
    Asat_tree = leaf_template.copy()
    Amax_tree = leaf_template.copy()
    Rd_tree = leaf_template.copy()
    Vcmax_tree = leaf_template.copy()
    Jmax_tree = leaf_template.copy()
    traits_flag = np.zeros(N_trees)

    for tt in range(0,N_trees):
        tree = TreeTag_census[tt]
        alt_tree = AltTag[tt]
        con1 = np.all([TreeTag_traits==tree,BALI_plot_traits==BALI_plot_census[tt]],axis=0)
        con2 = np.all([TreeTag_traits==alt_tree,BALI_plot_traits==BALI_plot_census[tt]],axis=0)
        #indices =np.any(np.asarray([TreeTag_traits==tree,TreeTag_traits==alt_tree]),axis=0)
        indices =np.any([con1,con2],axis=0)
        Plot[tt]=BALI_plot_census[tt]
        N_leaves = np.sum(indices)
        if N_leaves > 0:
            traits_flag[tt]=1
            print(tree, alt_tree, N_leaves, Plot[tt])
        if np.sum(N_leaves>10):
            print( "more leaves than expected for this tree \n-suggest increasing number of permitted leaves in code")

        for ll in range(0,N_leaves):
            shade_tag_tree[tt,ll] = shade_tag[indices][ll]
            leaf_height_tree[tt,ll] = leaf_height[indices][ll]
            leaf_thickness_tree[tt,ll] = leaf_thickness[indices][ll]
            SLA_tree[tt,ll] = SLA[indices][ll]
            LMA_tree[tt,ll] = LMA[indices][ll]
            Ci_tree[tt,ll] = Ci[indices][ll]
            VPD_tree[tt,ll] = VPD[indices][ll]
            LeafT_tree[tt,ll] = LeafT[indices][ll]
            gs_tree[tt,ll] = gs[indices][ll]
            Km_tree[tt,ll] = Km[indices][ll]
            c2co_tree[tt,ll] = c2co[indices][ll]
            Asat_tree[tt,ll] = Asat[indices][ll]
            Amax_tree[tt,ll] = Amax[indices][ll]
            Rd_tree[tt,ll] = Rd[indices][ll]
            Vcmax_tree[tt,ll] = Vcmax[indices][ll]
            Jmax_tree[tt,ll] = Jmax[indices][ll]

            # finally loop through the census data and find earliest tree height recorded - assume this represents tree height at time of traits data collection
            census_number = 0
            for cc in range(0,3):
                if np.isfinite(Height[tt,cc]):
                    census_number += 1
                    if census_number == 1:
                        canopy_position_tree[tt,ll] = leaf_height[indices][ll]/Height[tt,cc]
                        canopy_depth[tt,ll] = Height[tt,cc]-leaf_height[indices][ll]
            if census_number == 0:
                if np.isfinite(tree_height[indices][ll]): # this is the case whereby the tree height has not been measured in the census, but is present in the traits data.  Not sure when this would ever happen, but just in case
                    canopy_position_tree[tt,ll] = leaf_height[indices][ll]/tree_height[indices][ll]
                    canopy_depth[tt,ll] = tree_height[indices][ll]-leaf_height[indices][ll]


    # Now loop through all the leaf trait trees in case there are trees here that are not currently present in the tree census, for example if data is present from a plot in one dataset, but not in the other.  Ultimately this should be limited to plots for which we don't have census data, but this may not be the case with a small number of additional trees.
    print("==================================================")
    print("now finding traits trees not in the census")
    MissingTreeTags = []
    MissingTreePlots = []

    # get the number of plots for which we have traits data
    TraitPlots = np.unique(BALI_plot_traits)
    print(TraitPlots)
    for pp in range(0,TraitPlots.size):
        # get census data for plot
        TreeTag_census_pp = TreeTag_census[BALI_plot_census==TraitPlots[pp]]
        AltTag_pp= AltTag[BALI_plot_census==TraitPlots[pp]]
        # Now loop through this looking for tags which are missing from census
        TraitTrees = np.unique(TreeTag_traits[BALI_plot_traits==TraitPlots[pp]])
        N_traits_pp = TraitTrees.size
        for i in range(0,N_traits_pp):
            # Check if there is a census record for every tree in traits record
            if not np.any([TreeTag_census_pp==TraitTrees[i],AltTag_pp==TraitTrees[i]]):
                MissingTreeTags.append(TraitTrees[i])
                MissingTreePlots.append(TraitPlots[pp])
    missing_trees=np.asarray(MissingTreeTags)
    missing_tree_plots = np.asarray(MissingTreePlots)
    N_missing = missing_trees.size

    # now loop through the missing trees in the census dataset
    leaf_template = np.zeros((N_missing,10))*np.nan
    shade_tag_tree_new = leaf_template.copy()
    leaf_height_tree_new = leaf_template.copy()
    canopy_position_tree_new = leaf_template.copy()
    canopy_depth_new = leaf_template.copy()
    leaf_thickness_tree_new = leaf_template.copy()
    SLA_tree_new = leaf_template.copy()
    LMA_tree_new = leaf_template.copy()
    Ci_tree_new = leaf_template.copy()
    VPD_tree_new = leaf_template.copy()
    LeafT_tree_new = leaf_template.copy()
    gs_tree_new = leaf_template.copy()
    Km_tree_new = leaf_template.copy()
    c2co_tree_new = leaf_template.copy()
    Asat_tree_new = leaf_template.copy()
    Amax_tree_new = leaf_template.copy()
    Rd_tree_new = leaf_template.copy()
    Vcmax_tree_new = leaf_template.copy()
    Jmax_tree_new = leaf_template.copy()
    # Missing variables (missing census data)
    AltTag_new = np.zeros(N_missing)*np.nan
    Plot_new = np.zeros(N_missing,dtype = 'S8')
    Subplot_new = np.zeros(N_missing)*np.nan
    CensusDates_new = np.zeros((N_missing,3),dtype = 'datetime64[D]')
    DPOM_new = np.zeros((N_missing,3))*np.nan
    HPOM_new = np.zeros((N_missing,3))*np.nan
    Height_new = np.zeros((N_missing,3))*np.nan
    C_stem_new = np.zeros((N_missing,3))*np.nan
    C_coarse_root_new = np.zeros((N_missing,3))*np.nan
    RAINFOR_flag_new = np.zeros((N_missing,3),dtype='S8')
    Alive_flag_new = np.zeros((N_missing,3))*np.nan
    Species_new = np.zeros(N_missing,dtype='S32')
    SubplotCoords_new =  np.zeros((N_missing,2))*np.nan
    WoodDensity_new = np.zeros(N_missing)*np.nan
    traits_flag_new = np.ones(N_missing)
    for tt in range(0,N_missing):
        indices = np.all([TreeTag_traits==missing_trees[tt],BALI_plot_traits==missing_tree_plots[tt]],axis=0)
        Plot_new[tt]=missing_tree_plots[tt]

        N_leaves = np.sum(indices)
        print(Plot_new[tt], missing_trees[tt], N_leaves, np.sum(TreeTag_traits==missing_trees[tt]))
        if np.sum(N_leaves>10):
            print("more leaves than expected for this tree \n-suggest increasing number of permitted leaves in code")

        for ll in range(0,N_leaves):
            shade_tag_tree_new[tt,ll] = shade_tag[indices][ll]
            leaf_height_tree_new[tt,ll] = leaf_height[indices][ll]
            leaf_thickness_tree_new[tt,ll] = leaf_thickness[indices][ll]
            SLA_tree_new[tt,ll] = SLA[indices][ll]
            LMA_tree_new[tt,ll] = LMA[indices][ll]
            Ci_tree_new[tt,ll] = Ci[indices][ll]
            VPD_tree_new[tt,ll] = VPD[indices][ll]
            LeafT_tree_new[tt,ll] = LeafT[indices][ll]
            gs_tree_new[tt,ll] = gs[indices][ll]
            Km_tree_new[tt,ll] = Km[indices][ll]
            c2co_tree_new[tt,ll] = c2co[indices][ll]
            Asat_tree_new[tt,ll] = Asat[indices][ll]
            Amax_tree_new[tt,ll] = Amax[indices][ll]
            Rd_tree_new[tt,ll] = Rd[indices][ll]
            Vcmax_tree_new[tt,ll] = Vcmax[indices][ll]
            Jmax_tree_new[tt,ll] = Jmax[indices][ll]
            if np.isfinite(tree_height[indices][ll]):
                canopy_position_tree_new[tt,ll] = leaf_height[indices][ll]/tree_height[indices][ll]
                canopy_depth[tt,ll] = tree_height[indices][ll]-leaf_height[indices][ll]
                Height_new[tt] = tree_height[indices][ll]
            else:
                canopy_position_tree_new[tt,ll] = np.nan
                canopy_depth[tt,ll] = np.nan

    # Now wrap everything into a dictionary
    TreeDict = {}
    TreeDict['plot'] = np.concatenate((Plot,Plot_new),axis=0)
    TreeDict['subplot'] = np.concatenate((Subplot, Subplot_new),axis=0)
    TreeDict['CensusDates']=np.concatenate((CensusDates,CensusDates_new),axis=0)
    TreeDict['TreeTag']=np.concatenate((TreeTag_census,missing_trees),axis=0)
    TreeDict['AltTag']=np.concatenate((AltTag,AltTag_new),axis=0)
    TreeDict['DPOM']=np.concatenate((DPOM,DPOM_new),axis=0)
    TreeDict['HPOM']=np.concatenate((HPOM,HPOM_new),axis=0)
    TreeDict['TreeHeight']=np.concatenate((Height,Height_new),axis=0)
    TreeDict['C_stem']=np.concatenate((C_stem,C_stem_new),axis=0)
    TreeDict['C_coarse_roots']=np.concatenate((C_coarse_root,C_coarse_root_new),axis=0)
    TreeDict['C_wood']=np.concatenate((C_stem,C_stem_new),axis=0)+np.concatenate((C_coarse_root,C_coarse_root_new),axis=0)
    TreeDict['RAINFOR_flag'] = np.concatenate((RAINFOR_flag, RAINFOR_flag_new),axis=0)
    TreeDict['Alive_flag']=np.concatenate((Alive_flag,Alive_flag_new),axis=0)
    TreeDict['Species']=np.concatenate((Species,Species_new),axis=0)
    TreeDict['SubplotCoords']=np.concatenate((SubplotCoords,SubplotCoords_new),axis=0)
    TreeDict['WoodDensity']= np.concatenate((WoodDensity,WoodDensity_new),axis=0)
    TreeDict['LeafHeight']=np.concatenate((leaf_height_tree,leaf_height_tree_new),axis=0)
    TreeDict['LeafCanopyPosition'] =np.concatenate(( canopy_position_tree, canopy_position_tree_new),axis=0)
    TreeDict['LeafCanopyDepth'] =np.concatenate(( canopy_depth, canopy_depth_new),axis=0)
    TreeDict['ShadeTag'] =np.concatenate( (shade_tag_tree,shade_tag_tree_new),axis=0)
    TreeDict['LeafThickness']=np.concatenate((leaf_thickness_tree,leaf_thickness_tree_new),axis=0)
    TreeDict['SLA']=np.concatenate((SLA_tree,SLA_tree_new),axis=0)
    TreeDict['LMA']=np.concatenate((LMA_tree,LMA_tree_new),axis=0)
    TreeDict['Ci']=np.concatenate((Ci_tree,Ci_tree_new),axis=0)
    TreeDict['VPD']=np.concatenate((VPD_tree,VPD_tree_new),axis=0)
    TreeDict['LeafT']=np.concatenate((LeafT_tree,LeafT_tree_new),axis=0)
    TreeDict['gs']=np.concatenate((gs_tree,gs_tree_new),axis=0)
    TreeDict['Km']=np.concatenate((Km_tree,Km_tree_new),axis=0)
    TreeDict['c2co']=np.concatenate((c2co_tree,c2co_tree_new),axis=0)
    TreeDict['Asat']=np.concatenate((Asat_tree,Asat_tree_new),axis=0)
    TreeDict['Amax']=np.concatenate((Amax_tree,Amax_tree_new),axis=0)
    TreeDict['Rd']=np.concatenate((Rd_tree,Rd_tree_new),axis=0)
    TreeDict['Vcmax']=np.concatenate((Vcmax_tree,Vcmax_tree_new),axis=0)
    TreeDict['Jmax']=np.concatenate((Jmax_tree,Jmax_tree_new),axis=0)
    TreeDict['TraitAvailability']=np.concatenate((traits_flag,traits_flag_new),axis=0)


    # Make some additional calculations
    Ntrees,Ncens=TreeDict['C_wood'].shape
    C_wood_increment=TreeDict['C_wood'][:,1:]-TreeDict['C_wood'][:,:-1]
    C_wood_increment[C_wood_increment<0]=0 # reset wood increment to zero if negative
    TreeDict['C_wood_increment']=C_wood_increment

    plots = np.unique(TreeDict['plot'])
    N_plots = plots.size
    Dates_temp = TreeDict['CensusDates'].copy()
    for pp in range(0, N_plots):
        ii = TreeDict['plot']==[plots[pp]]
        for cc in range(0,Ncens):
            plot_dates = TreeDict['CensusDates'][ii,cc]
            plot_dates[plot_dates==np.datetime64(0,'D')]=np.max(plot_dates)
            Dates_temp[ii,cc] = plot_dates.copy()
    growth_interval = Dates_temp[:,1:]-Dates_temp[:,:-1]
    growth_interval[growth_interval<0]=0

    TreeDict['growth_interval'] = growth_interval

    return TreeDict

def collate_plot_level_census_data(census_file):
    BALI_plot, Subplot, CensusDates, TreeTag_census, AltTag, DPOM, HPOM, Height, C_stem, C_coarse_root, RAINFOR_flag, Alive_flag, Species, SubplotCoords, WoodDensity = read_ICP_census_data(census_file)

    plot_names = np.unique(BALI_plot)
    N_plots = plot_names.size

    CensusDict={}
    # Interesting properties: CanopyHeight, C_stem, C_coarse_roots, CensusDate
    for i in range(0,N_plots):
        plot_indices = BALI_plot==plot_names[i]
        #Set up arrays to save
        subplot_ids = np.unique(Subplot)
        n_subplots = subplot_ids.size
        dates = np.zeros((n_subplots,3),dtype = 'datetime64[D]')
        CanHt = np.zeros((n_subplots,3))
        Cstem = np.zeros((n_subplots,3))
        Croot = np.zeros((n_subplots,3))
        BasalArea = np.zeros((n_subplots,3))
        # note that for growth, mortality and recruitment, first year of census will be nan values because no there are no previous surveys!
        Growth = np.zeros((n_subplots,3))*np.nan
        Mortality = np.zeros((n_subplots,3))*np.nan
        Recruitment = np.zeros((n_subplots,3))*np.nan

        for s in range(0,n_subplots):
            subplot_indices = plot_indices * Subplot==subplot_ids[s]
            for y in range(0,3):
                datetemp = CensusDates[subplot_indices,y]
                dates[s,y]= np.max(datetemp)

                Cstemtemp = C_stem[subplot_indices,y]
                if np.isfinite(Cstemtemp).sum()>0:
                    Cstem[s,y]= np.sum(Cstemtemp[np.isfinite(Cstemtemp)])
                else:
                    Cstem[s,y]=np.nan

                Croottemp = C_coarse_root[subplot_indices,y]
                if np.isfinite(Croottemp).sum()>0:
                    Croot[s,y]= np.sum(Croottemp[np.isfinite(Croottemp)])
                else:
                    Croot[s,y]= np.nan

                httemp = Height[subplot_indices,y]
                if np.isfinite(httemp).sum()>0:
                    CanHt[s,y]= np.mean(httemp[np.isfinite(httemp)])
                else:
                    CanHt[s,y]=np.nan

                DBHtemp = DPOM[subplot_indices,y]
                if np.isfinite(DBHtemp).sum()>0:
                    BasalArea[s,y]= np.pi*np.sum((DBHtemp[np.isfinite(DBHtemp)]/2)**2)
                else:
                    BasalArea[s,y]=np.nan

        # need to do catch for where there are no trees in subplot that has been surveyed!
        for s in range(0,n_subplots):
            for y in range(0,3):
                if dates[s,y]==np.datetime64('1970-01-01','D'):
                    if np.max(dates[:,y])>np.datetime64('1970-01-01','D'):
                        dates[s,y]=np.max(dates[:,y])
                        Croot[s,y]=0.
                        Cstem[s,y]=0.
        # now lets do the growth, mortality and recruitment
        for s in range(0,n_subplots):
            subplot_indices = plot_indices * Subplot==subplot_ids[s]
            Cwood_temp = C_stem[subplot_indices]+C_coarse_root[subplot_indices]
            for y in range(1,3):
                """
                growth_indices = np.all((np.isfinite(Cwood_temp[:,y-1]),np.isfinite(Cwood_temp[:,y])),axis=0)
                recruit_indices = np.all((np.isfinite(Cwood_temp[:,y]),~np.isfinite(Cwood_temp[:,y-1])),axis=0)
                mortality_indices = np.all((np.isfinite(Cwood_temp[:,y-1]),~np.isfinite(Cwood_temp[:,y])),axis=0)
                """
                growth_indices = np.all((Cwood_temp[:,y-1]>0,Cwood_temp[:,y]>0),axis=0)
                recruit_indices = np.all((Cwood_temp[:,y]>0,Cwood_temp[:,y-1]==0),axis=0)
                mortality_indices = np.all((Cwood_temp[:,y-1]>0,Cwood_temp[:,y]==0),axis=0)
                if np.isfinite(Cwood_temp).sum()>0:
                    Growth[s,y] = np.sum(Cwood_temp[:,y][growth_indices]-Cwood_temp[:,y-1][growth_indices])
                    Recruitment[s,y] = np.sum(Cwood_temp[:,y][recruit_indices])
                    Mortality[s,y] = np.sum(Cwood_temp[:,y-1][mortality_indices])
                else:
                     Growth[s,y] = np.nan
                     Recruitment[s,y] = np.nan
                     Mortality[s,y] = np.nan

        PlotDict = {}
        PlotDict['n_subplots']=n_subplots
        PlotDict['CanopyHeight']=CanHt
        PlotDict['C_stem']=Cstem
        PlotDict['C_coarse_roots']=Croot
        PlotDict['C_wood']=Croot+Cstem
        PlotDict['CensusDate']=dates
        PlotDict['BasalArea']=BasalArea
        PlotDict['Growth']= Growth
        PlotDict['Recruitment']= Recruitment
        PlotDict['Mortality']=Mortality
        CensusDict[plot_names[i]]=PlotDict
    return CensusDict

# Read litterfall from GEM plot database file (converted to csv).  Mass collections (denoted with prefix m) are in g accumulated.  Otherwise given as flux in Mg C ha-1 yr-1
def read_litterfall_data(litter_file):

    datatype = {'names': ('ForestType', 'Plot', 'CollectionDate', 'PreviousCollectionDate',
                'AccumulationDays', 'Trap', 'TrapSize', 'mLeaves', 'mTwigs', 'mFruit',
                'mFlowers', 'mSeeds', 'mMisc','Comments','Leaves', 'Twigs', 'Fruit',
                'Flowers', 'Seeds', 'Misc', 'Reproductive','Total'),
                'formats': ('S16','S16','S10','S10','f16','i8','f16','f16','f16','f16',
                'f16','f16','f16','S64','f16','f16','f16','f16','f16','f16','f16','f16')}
    litter_data = np.genfromtxt(litter_file, skip_header = 1, delimiter = ',',dtype=datatype)

    plot_names = np.unique(litter_data['Plot'])
    N_plots = plot_names.size

    LitterDict = {}

    for i in range(0,N_plots):
        plot_data = litter_data[litter_data['Plot']==plot_names[i]]
        N_sub = np.max(np.unique(plot_data['Trap']))
        PlotDict = {}
        CollectionDates = []
        PreviousDates = []
        AccumulationDays = []
        TrapSize = []
        mLeaves = []
        mTwigs = []
        mFruit = []
        mFlowers = []
        mSeeds = []
        mMisc = []
        rLeaves = []
        rTwigs = []
        rFruit = []
        rFlowers = []
        rSeeds = []
        rRepro = []
        rTotal = []
        rMisc = []

        for j in range(0,N_sub):
            subplot_data = plot_data[plot_data['Trap']==j+1]
            N_collections = subplot_data.size

            cDates = np.zeros(N_collections, dtype = 'datetime64[D]')
            pDates = np.zeros(N_collections, dtype = 'datetime64[D]')

            # mass collected by component - note that masses need converting from g per trap to g m-2
            mLeaves.append(subplot_data['mLeaves']/subplot_data['TrapSize'])
            mTwigs.append(subplot_data['mTwigs']/subplot_data['TrapSize'])
            mFruit.append(subplot_data['mFruit']/subplot_data['TrapSize'])
            mFlowers.append(subplot_data['mFlowers']/subplot_data['TrapSize'])
            mSeeds.append(subplot_data['mSeeds']/subplot_data['TrapSize'])
            mMisc.append(subplot_data['mMisc']/subplot_data['TrapSize'])
            # Rates of litter fall by component
            rLeaves.append(subplot_data['Leaves'])
            rTwigs.append(subplot_data['Twigs'])
            rFruit.append(subplot_data['Fruit'])
            rFlowers.append(subplot_data['Flowers'])
            rSeeds.append(subplot_data['Seeds'])
            rMisc.append(subplot_data['Misc'])
            rRepro.append(subplot_data['Reproductive'])
            rTotal.append(subplot_data['Total'])
            AccumulationDays.append(subplot_data['AccumulationDays'])
            for k in range(0, N_collections):
                if len(subplot_data['CollectionDate'][k])>0:
                    day,month,year = subplot_data['CollectionDate'][k].split("/")
                    cDates[k]= np.datetime64(year+'-'+month+'-'+day)
                else:
                    cDates[k]= np.datetime64(0,'D')#

                if len(subplot_data['PreviousCollectionDate'][k])>0:
                    day,month,year = subplot_data['PreviousCollectionDate'][k].split("/")
                    pDates[k]= np.datetime64(year+'-'+month+'-'+day)
                else:
                    pDates[k]= np.datetime64(0,'D')
                #accDays[k]= cDates[k]-pDates[k]

            #AccumulationDays.append(accDays)
            PreviousDates.append(pDates)
            CollectionDates.append(cDates)
            TrapSize.append(subplot_data['TrapSize'][0])

        # Load plot data into plot dictionary
        PlotDict['N_Subplots']=N_sub
        PlotDict['CollectionDate']=np.asarray(CollectionDates)
        PlotDict['PreviousCollectionDate']=np.asarray(PreviousDates)
        PlotDict['AccumulationDays']=np.asarray(AccumulationDays)
        PlotDict['mLeaves']=np.asarray(mLeaves)
        PlotDict['mTwigs']=np.asarray(mTwigs)
        PlotDict['mFruit']=np.asarray(mFruit)
        PlotDict['mFlowers']=np.asarray(mFlowers)
        PlotDict['mSeeds']=np.asarray(mSeeds)
        PlotDict['mMisc']=np.asarray(mMisc)
        PlotDict['mTotal']=PlotDict['mLeaves']+PlotDict['mTwigs']+PlotDict['mFruit']+PlotDict['mFlowers']+PlotDict['mSeeds']+PlotDict['mMisc']
        PlotDict['rLeaves']=np.asarray(rLeaves)
        PlotDict['rTwigs']=np.asarray(rTwigs)
        PlotDict['rFruit']=np.asarray(rFruit)
        PlotDict['rFlowers']=np.asarray(rFlowers)
        PlotDict['rSeeds']=np.asarray(rSeeds)
        PlotDict['rMisc']=np.asarray(rMisc)
        PlotDict['rReproductive']=np.asarray(rRepro)
        PlotDict['rTotal']=np.asarray(rTotal)
        PlotDict['TrapSize']=np.asarray(TrapSize)

        LitterDict[plot_names[i]]=PlotDict

    return LitterDict


def read_soil_respiration_data(soil_resp_file):
    datatype = {'names': ('ForestType', 'Plot', 'Date', 'Observer', 'Collar', 'SoilMoisture', 'SoilT', 'AirT', 'RespirationFlux', 'Remarks', 'QualityFlag'), 'formats': ('S16','S16','S10','S32','i8','f16','f16','f16','f16','S64','i8')}
    resp_data = np.genfromtxt(soil_resp_file, skip_header = 1, delimiter = ',',dtype=datatype)

    plot_names = np.unique(resp_data['Plot'])
    N_plots = plot_names.size

    SoilRespDict = {}

    for i in range(0,N_plots):
        plot_data = resp_data[resp_data['Plot']==plot_names[i]]
        N_collars = np.max(np.unique(resp_data['Collar']))
        PlotDict = {}
        Dates = []
        SoilMoisture = []
        SoilT = []
        AirT = []
        Flux = []
        Flag = []

        for j in range(0,N_collars):
            collar_data = plot_data[plot_data['Collar']==j+1]
            N_collections = collar_data.size
            cDates = np.zeros(N_collections, dtype = 'datetime64[D]')

            SoilMoisture.append(collar_data['SoilMoisture'])
            SoilT.append(collar_data['SoilT'])
            AirT.append(collar_data['AirT'])
            Flux.append(collar_data['RespirationFlux'])
            Flag.append(collar_data['QualityFlag'])

            for k in range(0, N_collections):
                day,month,year = collar_data['Date'][k].split("/")
                cDates[k]= np.datetime64(year+'-'+month+'-'+day)

            Dates.append(cDates)

        # Load plot data into plot dictionary
        PlotDict['N_Collars']=N_collars
        PlotDict['CollectionDate']=np.asarray(Dates)
        PlotDict['SoilMoisture']=np.asarray(SoilMoisture)
        PlotDict['SoilT']=np.asarray(SoilT)
        PlotDict['AirT']=np.asarray(AirT)
        PlotDict['RespirationFlux']=np.asarray(Flux)
        PlotDict['QualityFlag']=np.asarray(Flag)

        SoilRespDict[plot_names[i]]=PlotDict

    return SoilRespDict

def read_soil_stocks_and_npp(roots_file):
    datatype = {'names': ('ForestType', 'Plot', 'Core', 'DataType', 'CollectionDate', 'PreviousCollectionDate', 'AccumulationDays', 'CoarseRoots', 'FineRoots', 'Remarks'), 'formats': ('S16','S16','i8','S5','S10','S10','i8','f16','f16','S64')}
    roots_data = np.genfromtxt(roots_file, skip_header = 1, delimiter = ',',dtype=datatype)

    plot_names = np.unique(roots_data['Plot'])
    N_plots = plot_names.size

    RootStocksDict = {}
    RootNPPDict = {}

    for i in range(0,N_plots):
        plot_data = roots_data[roots_data['Plot']==plot_names[i]]
        NPP_data = plot_data[plot_data['DataType']=='NPP']
        Stocks_data = plot_data[plot_data['DataType']=='Stock']
        N_cores = np.max(np.unique(roots_data['Core']))
        StocksDict = {}
        NPPDict = {}

        # First do stocks
        cDates = np.unique(Stocks_data['CollectionDate'])
        N_stocks = cDates.size
        StocksDates = np.zeros(N_stocks,dtype='datetime64[D]')
        for dd in range(0,N_stocks):
            day,month,year = cDates[dd].split("/")
            StocksDates[dd]= np.datetime64(year+'-'+month+'-'+day)

        StocksDates = np.sort(StocksDates)

        FineRootsStocks = np.zeros((N_cores,N_stocks))*np.nan
        CoarseRootsStocks =np.zeros((N_cores,N_stocks))*np.nan

        for cc in range(0,N_cores):
            core_data = Stocks_data[Stocks_data['Core']==cc+1]
            N_collections = core_data.size

            for kk in range(0, N_collections):
                day,month,year = core_data['CollectionDate'][kk].split("/")
                dd = (StocksDates == np.datetime64(year+'-'+month+'-'+day))
                CoarseRootsStocks[cc,dd] = core_data['CoarseRoots'][kk]
                FineRootsStocks[cc,dd] = core_data['FineRoots'][kk]

        StocksDict['N_Cores']=N_cores
        StocksDict['CollectionDate']=StocksDates
        StocksDict['CoarseRootStocks']=CoarseRootsStocks
        StocksDict['FineRootStocks']=FineRootsStocks
        StocksDict['Mask'] = np.isfinite(FineRootsStocks)


        # Next do NPP
        cDates = np.unique(NPP_data['CollectionDate'])
        N_NPP = cDates.size
        CollectionDates = np.zeros(N_NPP,dtype='datetime64[D]')
        for dd in range(0,N_NPP):
            day,month,year = cDates[dd].split("/")
            CollectionDates[dd]= np.datetime64(year+'-'+month+'-'+day)

        CollectionDates = np.sort(CollectionDates)
        PreviousDates = np.zeros((N_cores,N_NPP),dtype='datetime64[D]')
        AccumulationDays = np.zeros((N_cores,N_NPP))*np.nan
        FineRootsNPP = np.zeros((N_cores,N_NPP))*np.nan
        CoarseRootsNPP = np.zeros((N_cores,N_NPP))*np.nan

        for cc in range(0,N_cores):
            core_data = NPP_data[NPP_data['Core']==cc+1]
            N_collections = core_data.size

            for kk in range(0, N_collections):
                day,month,year = core_data['CollectionDate'][kk].split("/")
                dd = (CollectionDates == np.datetime64(year+'-'+month+'-'+day))

                if len(core_data['PreviousCollectionDate'][kk])==0:
                    PreviousDates[cc,dd]= np.datetime64(0,'D')
                else:
                    pday,pmonth,pyear = core_data['PreviousCollectionDate'][kk].split("/")
                    PreviousDates[cc,dd]= np.datetime64(pyear+'-'+pmonth+'-'+pday)
                    AccumulationDays[cc,dd]=core_data['AccumulationDays'][kk]
                    CoarseRootsNPP[cc,dd]=core_data['CoarseRoots'][kk]
                    FineRootsNPP[cc,dd]=core_data['FineRoots'][kk]

        NPPDict['N_Cores']=N_cores
        NPPDict['CollectionDate']=CollectionDates
        NPPDict['PreviousCollectionDate']=PreviousDates
        NPPDict['AccumulationDays']=AccumulationDays
        NPPDict['CoarseRootNPP']=CoarseRootsNPP
        NPPDict['FineRootNPP']=FineRootsNPP
        NPPDict['Mask'] = np.isfinite(FineRootsNPP)

        RootStocksDict[plot_names[i]]=StocksDict
        RootNPPDict[plot_names[i]]=NPPDict

    return RootStocksDict,RootNPPDict




def write_leaf_traits_to_file(census_file,traits_file,branch_file,leaf_file):

    TreeData = collate_tree_level_data(census_file,traits_file,branch_file,leaf_file)
    ii = TreeData['TraitAvailability']==1
    plot = TreeData['plot'][ii]
    tree = TreeData['TreeTag'][ii]
    shade_tag = TreeData['ShadeTag'][ii,:]
    Asat = TreeData['Asat'][ii,:]
    Amax = TreeData['Amax'][ii,:]
    Rd = TreeData['Rd'][ii,:]
    Km = TreeData['Km'][ii,:]
    c2co = TreeData['c2co'][ii,:]
    Jmax = TreeData['Jmax'][ii,:]
    Vcmax= TreeData['Vcmax'][ii,:]
    Ci = TreeData['Ci'][ii,:]
    VPD = TreeData['VPD'][ii,:]
    LeafT= TreeData['LeafT'][ii,:]
    gs = TreeData['gs'][ii,:]
    SLA = TreeData['SLA'][ii,:]
    leaf_height= TreeData['LeafHeight'][ii,:]
    tree_height= TreeData['TreeHeight'][ii,:]
    canopy_pos= TreeData['LeafCanopyPosition'][ii,:]
    canopy_depth= TreeData['LeafCanopyDepth'][ii,:]



    # Write leaf traits to csv file
    print('writing traits to file')
    LeafTraitsFile = 'BALI_leaf_traits.csv'
    out = open(LeafTraitsFile,'w')
    out.write('Plot, TreeTag, ShadeTag, Asat, Amax, Rd, Km, c2co, Jmax, Vcmax, Ci, VPD, LeafT, gs, SLA, LMA, leaf_height, tree_height, canopy_position, canopy_depth\n')
    N_trees, N_leaves = Asat.shape
    for tt in range(0,N_trees):
        if np.sum(np.isfinite(tree_height[tt,:])):
            tree_ht = np.min(tree_height[tt,:][np.isfinite(tree_height[tt,:])])
        else:
            tree_ht=np.nan
        for ll in range(0,N_leaves):
            if np.isfinite(Asat[tt,ll]):
                out.write( plot[tt]+','+ str(tree[tt])+','+str(shade_tag[tt,ll])+','+str(Asat[tt,ll])+','+str(Amax[tt,ll])+','+str(Rd[tt,ll])+','+str(Km[tt,ll])+','+str(c2co[tt,ll])+','+str(Jmax[tt,ll])+','+str(Vcmax[tt,ll])+','+str(Ci[tt,ll])+','+str(VPD[tt,ll])+','+str(LeafT[tt,ll])+','+str(gs[tt,ll])+','+str(SLA[tt,ll])+','+str(1/SLA[tt,ll])+','+str(leaf_height[tt,ll])+','+str(tree_ht)+','+str(canopy_pos[tt,ll])+','+str(canopy_depth[tt,ll])+'\n' )

    out.close()


# get turnover rate estimates for wood, foliar and fine root biomass pools
def get_woody_biomass_info_for_SPA(census_file,traits_file,branch_file,leaf_file):

    # get plot trees
    tree_data = collate_tree_level_data(census_file,traits_file,branch_file,leaf_file)

    # first woody biomass turnover
    plots = np.unique(tree_data['plot'])
    npp=np.zeros(plots.size)*np.nan
    c_init=np.zeros(plots.size)*np.nan
    turnover=np.zeros(plots.size)*np.nan
    for pp in range(0,plots.size):
        ii = tree_data['plot']==plots[pp]
        increment = np.sum(tree_data['C_wood_increment'][ii,:],axis=1) # this gives the biomass increment in units of kg C
        interval = np.sum(tree_data['growth_interval'][ii,:],axis=1) # this should give the time interval in days
        av_biomass = np.mean(np.sum(tree_data['C_wood'][ii,:],axis=0)) # this is the plot biomass for each census in units of kg C
        npp[pp] = np.sum(increment/np.asarray(interval,dtype='float'))*10**3/10**4  # npp in kg C m-2 d-1

        turnover[pp] = npp[pp]/av_biomass # note that this assumes steady state - likely to be false especially for regenerating forest
        c_init[pp] = np.sum(tree_data['C_wood'][ii,:],axis=0)[0]*1000/10000 # g m-2
        print("===============================")
        print("PLOT: ", plots[pp])
        print("C_wood_i: ", c_init[pp])
        print("biomass increase: ", np.sum(increment))
        print("NPP: ", npp[pp])
        print("turnover: ", turnover[pp])
        print("===============================")
    return plots, c_init, npp, turnover


# Load spreadsheet of LAI derived from hemispherical photographs (courtesy of Terhi Riutta at Oxford).  LAI estimated using Hemisfer.
def load_field_LAI(LAI_file):
    datatype = {'names': ('ForestType','Plot', 'Subplot', 'LAI'), 'formats': ('S32','S32'
,'i8','f16')}
    hemisfer_LAI = np.genfromtxt(LAI_file, skip_header = 1, delimiter = ',',dtype=datatype)

    return hemisfer_LAI

# get foliar biomass estimates from LAI and average LCA.
def get_foliage_info_for_SPA(traits_file,branch_file,leaf_file,lai_file,litter_file):
    # some loading of data
    litter = read_litterfall_data(litter_file)
    BALI_plot_traits, TreeTag_traits, leaf_ID, shade_tag, tree_height, leaf_height, leaf_thickness, leaf_area, SLA, LMA, Ci, VPD, LeafT, gs, Km, c2co, Asat, Amax, Rd, Vcmax, Jmax = load_leaf_traits(traits_file, branch_file, leaf_file)
    BALI_plot_traits[BALI_plot_traits=='BEL']='Belian'
    BALI_plot_traits[BALI_plot_traits=='SER']='Seraya'
    BALI_plot_traits[BALI_plot_traits=='DAS1']='DC1'
    BALI_plot_traits[BALI_plot_traits=='DAF2']='DC2'

    LAI_hemisfer = load_field_LAI(lai_file)

    # declare some stuff
    plots = litter.keys()
    NPPfol = np.zeros(len(plots))
    LAI_hemisfer = load_field_LAI(lai_file)
    LAI = np.zeros(len(plots))
    LCA = np.zeros(len(plots))
    Cfol = np.zeros(len(plots))
    turnover = np.zeros(len(plots))

    for pp in range(0,len(plots)):
        # let's start with the litter fall.  Get time average fluxes for each subplot, which will be taken to represent the foliage NPP.  This assumes a steady state canopy i.e. leaves are produced at the same rate as they fall
        lit=litter[plots[pp]]
        lit['AccumulationDays'][np.isnan(lit['mTotal'])]=np.nan
        trap_size = np.mean(lit['TrapSize'])
        mLit = np.nansum(lit['mTotal'],axis=1)
        AccDays = np.nansum(lit['AccumulationDays'],axis=1)
        NPPfol_sp = mLit/AccDays/trap_size
        NPPfol[pp] = np.mean(NPPfol_sp)

        # Now get LAI from hemisfer
        LAI[pp] = np.mean(LAI_hemisfer['LAI'][LAI_hemisfer['Plot']==plots[pp]])

        # Now get LMA
        if np.sum(BALI_plot_traits==plots[pp])>0:
            LCA[pp] = np.mean(LMA[np.all((BALI_plot_traits==plots[pp],np.isfinite(LMA)),axis=0)])
        else:
            LCA[pp] = np.mean(LMA[np.isfinite(LMA)])

        # Now get Cfol in g m-2
        Cfol[pp] = LAI[pp]*LCA[pp]

        turnover[pp] = NPPfol[pp]/Cfol[pp]/2.

        print("============================")
        print("Plot: ", plots[pp])
        print("NPPfol: ", NPPfol[pp], " g m-2 d-1")
        print("LAI: ", LAI[pp])
        print("LCA: ", LCA[pp], " g m-2")
        print("Cfol: ", Cfol[pp]," g m-2")
        print("turnover: ", turnover[pp], " d-1")

    return plots, Cfol, NPPfol, turnover

def get_root_info_for_SPA(fine_roots_file):
    stocks,npp=read_soil_stocks_and_npp(fine_roots_file)
    plots = stocks.keys()
    NPP=np.zeros(len(plots))*np.nan
    c_root=np.zeros(len(plots))*np.nan
    turnover=np.zeros(len(plots))*np.nan
    core_area = np.pi*0.06**2
    for pp in range(0,len(plots)):
        #core_npp = np.nansum(npp[plots[pp]]['FineRootNPP']*npp[plots[pp]]['AccumulationDays'])/np.nansum(npp[plots[pp]]['AccumulationDays'],axis=1)
        ii = np.isfinite(npp[plots[pp]]['FineRootNPP'])
        core_npp = np.nansum(npp[plots[pp]]['FineRootNPP'],axis=1)/np.sum(ii,axis=1)

        core_npp*=10**2/356.25 # convert to g m-2 d-1 from Mg Ha-1 yr-1
        root_npp = np.mean(core_npp) # average across the cores
        ii = np.isfinite(stocks[plots[pp]]['FineRootStocks'])
        core_stocks = np.nansum(stocks[plots[pp]]['FineRootStocks'],axis=1)/np.sum(ii,axis=1)
        #core_stocks = stocks[plots[pp]]['FineRootStocks'][:,0]
        root_stocks = np.mean(core_stocks)/core_area # convert g per core to g m-2
        c_root[pp] = root_stocks
        NPP[pp]=np.sum(root_npp)
        turnover[pp] = root_npp/root_stocks
        print("===============================")
        print("PLOT: ", plots[pp])
        print("C_root: ", root_stocks)
        print("NPP: ", np.sum(root_npp))
        print("turnover: ", turnover[pp])
        print("===============================")
    return plots, c_root, NPP, turnover

# get allocation fractions (alloc to wood, remaining alloc to roots)
def get_allocation_fractions(census_file,traits_file,branch_file,leaf_file,lai_file,litter_file,fine_roots_file):

    plots_root, c_root, npp_root, turnover_root = get_root_info_for_SPA(fine_roots_file)
    plots_fol, c_fol, npp_fol, turnover_fol = get_foliage_info_for_SPA(traits_file,branch_file,leaf_file,lai_file,litter_file)
    plots_wood, c_wood, npp_wood, turnover_wood =  get_woody_biomass_info_for_SPA(census_file,traits_file,branch_file,leaf_file)

    plots = np.asarray(plots_fol)
    plots_w = np.asarray(plots_wood)
    plots_w[plots_w=='Danum1'] = 'DC1'
    plots_w[plots_w=='Danum2'] = 'DC2'
    plots_r = np.asarray(plots_root)

    for i in range(0,len(plots_fol)):
        npp_total = npp_root[plots_r==plots[i]]+npp_wood[plots_w==plots[i]]+npp_fol[i]
        npp_minus_fol = npp_root[plots_r==plots[i]]+npp_wood[plots_w==plots[i]]
        f_fol = npp_fol[i]/npp_total
        f_roo = npp_root[plots_r==plots[i]]/npp_minus_fol
        print("===============================")
        print("PLOT: ", plots_fol[i])
        print("f_fol: ", f_fol)
        print("f_root: ", f_roo)
        print("===============================")

    return 0


def load_LAI_time_series(LAI_file):
    datatype = {'names': ('ForestType', 'Plot', 'Date', 'Subplot', 'SkyConditions', 'Exposure', 'LAI', 'Remarks', 'Qflag', 'Qreason'), 'formats': ('S16','S16','S10','i8','S16','i8','f8','S64','i8','S64')}
    LAI_data = np.genfromtxt(LAI_file, skip_header = 1, delimiter = ',',dtype=datatype)

    plot_names = np.unique(LAI_data['Plot'])
    N_plots = plot_names.size
    plot_dict = {}
    for pp in range(0,N_plots):
        LAI_dict = {}
        plot_data = LAI_data[LAI_data['Plot']==plot_names[pp]]
        N_subplots = np.unique(plot_data['Subplot']).size
        # set up some arrays for data output
        date_str = np.unique(plot_data['Date'])
        N_dates = date_str.size
        dates = np.zeros(N_dates,dtype = 'datetime64[D]')
        for dd in range(0,N_dates):
            day,month,year = date_str[dd].split("/")
            dates[dd] = np.datetime64(year+'-'+month+'-'+day)
        dates=np.sort(dates)

        # Need to loop through the dates column in the input array, and produce equivalent but in np.datetime64 so that cross referencing is easier
        date_ref = np.zeros(plot_data['Date'].size,dtype = 'datetime64[D]')
        for dd in range(0,plot_data['Date'].size):
            day,month,year = plot_data['Date'][dd].split("/")
            date_ref[dd] =  np.datetime64(year+'-'+month+'-'+day)

        LAI = np.zeros((N_subplots,N_dates))
        for dd in range(0,N_dates):
            for ss in range(0,N_subplots):
                index=np.all((date_ref==dates[dd],plot_data['Subplot']==ss+1,plot_data['Exposure']==2,plot_data['Qflag']==1),axis=0)
                if np.sum(index)>0:
                    LAI[ss,dd] = plot_data['LAI'][index][0]
                else:
                    LAI[ss,dd] = np.nan

        LAI_dict['date']=dates.copy()
        LAI_dict['LAI']=LAI.copy()
        plot_dict[plot_names[pp]]=LAI_dict
    return plot_dict



if __name__ == "__main__":
    census_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_TreeCensus.csv'
    traits_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_traits_data/Photosynthesis_OG.csv'
    branch_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_traits_data/trees_OG.csv'
    leaf_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_traits_data/Leaf_area_OG.csv'
    fine_roots_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_FineRoots_Stock_NPP_RawData.csv'
    soil_respiration_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_TotalSoilRespiration_RawData.csv'
    litter_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_Litterfall_RawData.csv'
    lai_file='/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos.csv'
