import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import bottom_up_LAD_profiles as field


filename = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'

a, b, CF, r_sq, p, H, D = field.retrieve_crown_allometry_file(filename)

field_filename = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_carbonplots_FieldMapcensus2016.csv'

field_data = field.load_crown_survey_data(field_filename)
a_ht, b_ht, CF_ht, a_A, b_A, CF_A = field.calculate_allometric_equations_from_survey(field_data)

canopy_layers = np.arange(1,81)

plot = 'Belian'
plot_area = 10.**4
mask = field_data['plot']==plot
Ht,Area,Depth = field.calculate_crown_dimensions(field_data['DBH_field'][mask],field_data['Height_field'][mask],field_data['CrownArea'][mask], a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)

axis_a, axis_b, axis_c, z0 = field.construct_ellipses_for_subplot(Ht, Area, Depth)
LAD, CanopyV = field.calculate_LAD_profiles_ellipsoid(canopy_layers, axis_a, axis_b, axis_c, z0, plot_area)



plt.figure(1, facecolor='White',figsize=[9,3])
ax7 = plt.subplot2grid((1,3),(0,0))
ax7.set_ylabel('Crown Depth / m')
ax7.set_xlabel('Height / m')
ax7.annotate('a - Height-Crown Depth', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax7.plot(H,D,'.',color='blue',alpha=0.1)
H_mod = np.linspace(0.,H.max(),1000)
D_mod = CF*a*H_mod**b
ax7.plot(H_mod,D_mod,'-',color='black')
eq = '$D=%.3fH^{%.3f}$' % (CF*a, b)
ax7.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)

ax8 = plt.subplot2grid((1,3),(0,1))
ax8.annotate('b - DBH-Crown Area', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax8.set_xlabel('DBH / cm')
ax8.set_ylabel('Crown Area / m$^2$')
ax8.plot(field_data['DBH_field'],field_data['CrownArea'],'.',color='red',alpha=0.1)
DBH_mod = np.linspace(0.,np.nanmax(field_data['DBH_field']),1000)
CA_mod = CF_A*a_A*DBH_mod**b_A
ax8.plot(DBH_mod,CA_mod,'-',color='black')
eq = '$A=%.3fDBH^{%.3f}$' % (CF_A*a_A, b_A)
ax8.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)

ax9 = plt.subplot2grid((1,3),(0,2),sharex=ax8)
ax9.annotate('c - DBH-Height', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax9.set_xlabel('DBH / cm')
ax9.set_ylabel('Height / m')
ax9.plot(field_data['DBH_field'],field_data['Height_field'],'.',color='red',alpha=0.1)
Ht_mod = CF_ht*a_ht*DBH_mod**b_ht
ax9.plot(DBH_mod,Ht_mod,'-',color='black')
eq = '$H=%.3fDBH^{%.3f}$' % (CF_ht*a_ht, b_ht)
ax9.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)


plt.tight_layout()
plt.savefig('SAFE_allometries.png')
plt.show()
