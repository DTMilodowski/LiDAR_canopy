import numpy as np
import sensitivity_analysis_plots as sp

output_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/FiguresRevised/'


# Load in the sensitivity analysis results
max_height=80
heights = np.arange(0.,max_height)+1

resolution_MH_Belian = np.load('MH_resolution_sensitivity_adaptive_Belian.npy')[()]
resolution_MH_BNorth = np.load('MH_resolution_sensitivity_adaptive_BNorth.npy')[()]
resolution_MH_E = np.load('MH_resolution_sensitivity_adaptive_E.npy')[()]

resolution_rad_Belian = np.load('rad2_resolution_sensitivity_adaptive_Belian.npy')[()]
resolution_rad_BNorth = np.load('rad2_resolution_sensitivity_adaptive_BNorth.npy')[()]
resolution_rad_E = np.load('rad2_resolution_sensitivity_adaptive_E.npy')[()]

penetration_lim_Belian = np.load('penetration_limit_resolution_adaptive_Belian.npy')[()]
penetration_lim_BNorth = np.load('penetration_limit_resolution_adaptive_BNorth.npy')[()]
penetration_lim_E = np.load('penetration_limit_resolution_adaptive_E.npy')[()]

density_MH_Belian = np.load('MH_density_sensitivity_Belian.npy')[()]
density_MH_BNorth = np.load('MH_density_sensitivity_BNorth.npy')[()]
density_MH_E = np.load('MH_density_sensitivity_E.npy')[()]

density_rad_Belian = np.load('rad2_density_sensitivity_Belian.npy')[()]
density_rad_BNorth = np.load('rad2_density_sensitivity_BNorth.npy')[()]
density_rad_E = np.load('rad2_density_sensitivity_E.npy')[()]


#-------------------------------
# RESULTS - SENSITIVITY ANALYSIS
#-------------------------------
"""
 Figure 8 - Sensitivity analysis of vertical profiles to spatial resolution
 Comparison of OG vs Moderately Logged vs. Heavily Logged
"""
figure_number = 8
figure_name = output_dir + "fig8_profile_sensitivity_to_resolution_adaptive.png"
sp.plot_profile_sensitivity_resolution(figure_number,figure_name,heights,resolution_MH_Belian,
                            resolution_MH_BNorth,resolution_MH_E,resolution_rad_Belian,
                            resolution_rad_BNorth,resolution_rad_E)

"""
Figure 9 - Sensitivity analysis of unsampled voxels
"""
figure_number = 9
figure_name = output_dir + "figs9_unsampled_voxels.png"
sp.plot_penetration_limits(figure_number,figure_name,heights,penetration_lim_Belian,
                            penetration_lim_E,penetration_lim_BNorth)

"""
# Figure 10 - Sensitivity analysis of vertical profiles to point density
"""
figure_number = 10
figure_name = output_dir + "fig10_profile_sensitivity_to_pulse_density.png"
sp.plot_profile_sensitivity_density(figure_number,figure_name,heights,density_MH_Belian,
                            density_MH_BNorth,density_MH_E,density_rad_Belian,
                            density_rad_BNorth,density_rad_E)

"""
# Figure 11 - Summary plots from sensitivity analysis for PAI vs. resolution
"""
figure_number = 11
figure_name = output_dir + "fig11_PAI_sensitivity_resolution.png"
sp.plot_PAI_sensitivity_resolution(figure_number,figure_name,resolution_MH_Belian,
                        resolution_MH_BNorth,resolution_MH_E,resolution_rad_Belian,
                        resolution_rad_BNorth,resolution_rad_E)
"""
# Figure 12 - Summary plots from sensitivity analysis for PAI vs. pulse density
"""
figure_number = 12
figure_name = output_dir + "fig12_PAI_sensitivity_density.png"
sp.plot_PAI_sensitivity_density(figure_number,figure_name,density_MH_Belian,density_MH_BNorth,
                        density_MH_E,density_rad_Belian,
                        density_rad_BNorth,density_rad_E)
#-------------------------------
# SUPPLEMENT
# RESULTS
#-------------------------------
"""
# Figure S5 - sensitivity analysis, confidence interval sensitivity to resolution
"""
figure_number = 115
figure_name = output_dir + "figS5_profile_sensitivity_to_resolution_individual_CI.png"
sp.plot_profile_sensitivity_to_resolution_individual_CI(figure_number,figure_name,heights,
                    PAD_profiles_MH_Belian,PAD_profiles_rad_Belian)
"""
# Figure S6 - sensitivity analysis, confidence interval sensitivity to density
"""
figure_number = 116
figure_name = output_dir + "figS6_profile_sensitivity_to_point_density_individual_CI.png"
sp.plot_profile_sensitivity_to_point_density_individual_CI(figure_number,figure_name,heights,
                    PAD_profiles_MH_Belian,PAD_profiles_rad_Belian)
