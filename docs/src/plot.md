# Plot
I'm used to using Matplotlib for visualization, HDF5 files can be read by h5py package.


## S0\_read\_JLD2.py
```python
import numpy as np
import h5py
import os
import glob

#%% 1. set parameter
file = "../../data/BJ.081_BJ.084__2020_04_11_00_00_00T2021_04_13_00_00_00__all.jld2"
chan = "NN"
dt = 0.005 #?

#%% 2. read h5
# open file
f = h5py.File(file,'r')
# read data
data = f[chan]["data"][0] 
# read parameters
azi = f[chan]["azi"][()]
baz = f[chan]["baz"][()]
maxlag = f[chan]["maxlag"][()]
cc_len = f[chan]["cc_len"][()]
cc_step = f[chan]["cc_step"][()]
corr_type = f[chan]["corr_type"][()]
comp = f[chan]["comp"][()]
dist = f[chan]["dist"][()]   # dist = f[chan]["dist"].value
lat = f[chan]["lat"][()]
lon = f[chan]["lon"][()]
N_glob = f[chan]["N_glob"][()]
N_read = f[chan]["N_read"][()]
N_good = f[chan]["N_good"][()]
name = f[chan]["name"][()][0].decode('utf-8')
# close h5-file
f.close()
```


## S1\_plot\_waveform.py
```python
import obspy
import numpy as np
import h5py
import os
import glob
import matplotlib
import matplotlib.pyplot as plt

#%% 1. set parameter
path = "/Users/yf/3.Project/2-1.Velocity_Change/STACK_das3/stack_mean"
net_sta = "BJ.106"
chan = "NN"
dt = 0.005

#%% 2. read h5
files = sorted( glob.glob(os.path.join(path,net_sta+"/*/*")))
data_list = []
dist_list = []
name_list = []
for i in range(0,len(files)):
    f = h5py.File(files[i],'r')
    data = f[chan]["data"][0] 
    dist = f[chan]["dist"][()]
    lat = f[chan]["lat"][()]
    lon = f[chan]["lon"][()]
    N_glob = f[chan]["N_glob"][()]
    N_read = f[chan]["N_read"][()]
    N_good = f[chan]["N_good"][()]
    name = f[chan]["name"][()][0].decode('utf-8')
    data_list.append(data)
    name_list.append(name)
    dist_list.append(dist)
    f.close()

#%% 3. plot
tt = np.arange(-int(5),int(5),dt)
dist=0
for i in range(0,len(files)):
    dist += 1
    data = data_list[i]/np.max(data_list[i],axis=0)
    plt.plot(tt,data[3000:5000]+dist,'k',linewidth=0.8)

plt.xlabel('Time (s)')
plt.ylabel('Offset (m)')
#plt.savefig("seisnoise2.pdf", format='pdf', dpi=400)
```