from dipole_lib import importDipoleData, exportDipoleData
import numpy as np

s_num = 0
t_num = 100
folderName = "/home/jirka/JESS_2022/TRAJ"
oFolderName = "/home/jirka/JESS_2022/TRAJ"
for i in range(s_num, t_num):
    print(i)
    # import dipoles from permanent charges
    dipFileName = "%s/dip_perm_aligned_%i.csv" % (folderName, i)
    dip_data, V = importDipoleData(dipFileName)
    dip_array = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
    time_array = np.array(dip_data[[' time']])
    # import total dipoles of Amoeba ff
    amoebaDipFn = "%s/rot_tot_%i.dip" % (folderName, i)
    tot_dip = np.loadtxt(amoebaDipFn)
    # sum them
    sum_dip = dip_array + tot_dip
    sum_norm = np.linalg.norm(dip_array, axis=1, keepdims=True)
    out_dip_data = np.hstack((time_array, dip_array, sum_norm))
    outFileName = "%s/SUMDIP_%i.csv" % (oFolderName, i)
    exportDipoleData(out_dip_data, V, outFileName)
