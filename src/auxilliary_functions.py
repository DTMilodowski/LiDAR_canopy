import numpy as np

# Load spreadsheet of LAI derived from hemispherical photographs (courtesy of Terhi Riutta at Oxford).  LAI estimated using Hemisfer.
def load_field_LAI(LAI_file):
    datatype = {'names': ('ForestType','Plot', 'Subplot', 'LAI'), 'formats': ('S32','S32','i8','f16')}
    hemisfer_LAI = np.genfromtxt(LAI_file, skiprows = 1, delimiter = ',',dtype=datatype)

    return hemisfer_LAI
