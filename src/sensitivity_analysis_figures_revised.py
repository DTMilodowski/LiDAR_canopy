import numpy as np
import sensitivity_analysis_plots as sp

output_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/FiguresRevised/201910/'


# Load in the sensitivity analysis results
max_height=80
heights = np.arange(0.,max_height)+1

resolution_wt_Belian = np.load('MH_resolution_wt_sensitivity_adaptive_Belian.npy')[()]
resolution_wt_BNorth = np.load('MH_resolution_wt_sensitivity_adaptive_BNorth.npy')[()]
resolution_wt_E = np.load('MH_resolution_wt_sensitivity_adaptive_E.npy')[()]

resolution_MH_Belian = np.load('MH_resolution_sensitivity_adaptive_Belian.npy')[()]
resolution_MH_BNorth = np.load('MH_resolution_sensitivity_adaptive_BNorth.npy')[()]
resolution_MH_E = np.load('MH_resolution_sensitivity_adaptive_E.npy')[()]

resolution_rad_Belian = np.load('rad2_resolution_sensitivity_adaptive_Belian_sensitivity2.npy')[()]
resolution_rad_BNorth = np.load('rad2_resolution_sensitivity_adaptive_BNorth_sensitivity2.npy')[()]
resolution_rad_E = np.load('rad2_resolution_sensitivity_adaptive_E_sensitivity2.npy')[()]

penetration_lim_Belian = np.load('penetration_limit_resolution_adaptive_Belian.npy')[()]
penetration_lim_BNorth = np.load('penetration_limit_resolution_adaptive_BNorth.npy')[()]
penetration_lim_E = np.load('penetration_limit_resolution_adaptive_E.npy')[()]

density_MH_Belian = np.load('MH_density_sensitivity_Belian.npy')[()]
density_MH_BNorth = np.load('MH_density_sensitivity_BNorth.npy')[()]
density_MH_E = np.load('MH_density_sensitivity_E.npy')[()]

density_rad_Belian = np.load('rad2_density_sensitivity_Belian_sensitivity2.npy')[()]
density_rad_BNorth = np.load('rad2_density_sensitivity_BNorth_sensitivity2.npy')[()]
density_rad_E = np.load('rad2_density_sensitivity_E_sensitivity2.npy')[()]

density_wt_Belian = np.load('MH_wt_density_sensitivity_Belian.npy')[()]
density_wt_BNorth = np.load('MH_wt_density_sensitivity_BNorth.npy')[()]
density_wt_E = np.load('MH_wt_density_sensitivity_E.npy')[()]

#-------------------------------
# RESULTS - SENSITIVITY ANALYSIS
#-------------------------------
"""
# Figure 5 - Sensitivity of PAI estimates to pulse density and resolution
"""
figure_name = output_dir + 'Fig5_PAI_sensitivity.png'
figure_number = 5
sp.plot_PAI_sensitivity(figure_number,figure_name,
                        resolution_MH_Belian,resolution_MH_BNorth,resolution_MH_E,
                        resolution_rad_Belian,resolution_rad_BNorth,resolution_rad_E,
                        resolution_wt_Belian,resolution_wt_BNorth,resolution_wt_E,
                        density_MH_Belian,density_MH_BNorth,density_MH_E,
                        density_rad_Belian,density_rad_BNorth,density_rad_E,
                        density_wt_Belian,density_wt_BNorth,density_wt_E)

"""
Figure 6 - Sensitivity analysis of unsampled voxels
"""
figure_number = 6
figure_name = output_dir + "fig6_unsampled_voxels.png"
sp.plot_penetration_limits(figure_number,figure_name,heights,penetration_lim_Belian,
                            penetration_lim_E,penetration_lim_BNorth)

"""
#-------------------------------
# SUPPLEMENT
#-------------------------------
# Figure S7 - Sensitivity analysis of vertical profiles to point density
"""
figure_number = 117
figure_name = output_dir + "figS7_profile_sensitivity_to_pulse_density.png"
sp.plot_profile_sensitivity_density(figure_number,figure_name,heights,density_MH_Belian,
                            density_MH_BNorth,density_MH_E,density_rad_Belian,
                            density_rad_BNorth,density_rad_E,density_wt_Belian,
                            density_wt_BNorth,density_wt_E)

"""
# Figure S8 - sensitivity analysis, confidence interval sensitivity to density
"""
figure_number = 118
figure_name = output_dir + "figS8_profile_sensitivity_to_point_density_individual_CI.png"
sp.plot_profile_sensitivity_to_point_density_individual_CI(figure_number,figure_name,heights,
                    density_MH_Belian,density_rad_Belian,density_wt_Belian)

"""
 Figure S9 - Sensitivity analysis of vertical profiles to spatial resolution
 Comparison of OG vs Moderately Logged vs. Heavily Logged
"""
figure_number = 119
figure_name = output_dir + "figS9_profile_sensitivity_to_resolution_adaptive_sensitivity.png"
sp.plot_profile_sensitivity_resolution_full(figure_number,figure_name,heights,resolution_MH_Belian,
                            resolution_MH_BNorth,resolution_MH_E,resolution_rad_Belian,
                            resolution_rad_BNorth,resolution_rad_E,resolution_wt_Belian,
                            resolution_wt_BNorth,resolution_wt_E)

"""
# Figure S10 - sensitivity analysis, confidence interval sensitivity to resolution
"""
figure_number = 1110
figure_name = output_dir + "figS10_profile_sensitivity_to_resolution_individual_CI.png"
sp.plot_profile_sensitivity_to_resolution_individual_CI(figure_number,figure_name,heights,
                    resolution_MH_Belian,resolution_rad_Belian,resolution_wt_Belian)
