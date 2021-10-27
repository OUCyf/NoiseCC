var documenterSearchIndex = {"docs":
[{"location":"FFT_CC.html#FFT_CC","page":"FFT_CC","title":"FFT_CC","text":"","category":"section"},{"location":"FFT_CC.html","page":"FFT_CC","title":"FFT_CC","text":"Use Multi-process and Multi-thread to accelerate FFT and CC from raw-data (sac/mseed) to corr-data (CORRData).","category":"page"},{"location":"FFT_CC.html","page":"FFT_CC","title":"FFT_CC","text":"############################################\n## 1.Set parameters \n############################################\ncores = 1\nthreads = 1\nflag = true\nprecompile = false\n...\n# Set your own parameters here...\n# Then run following code\n\n############################################\n# 2.run NoiseCC.FFT\n############################################\nusing Distributed\naddprocs(cores, topology = :master_worker)\n@everywhere using NoiseCC\n\nNoiseCC.FFT_CC(\n    # processes and threads parameters\n    cores,threads,flag,precompile,\n    # absolute path parameters\n    DATADIR,CCFDIR,\n    # select time-chunck folder for cross-correlate\n    select_time,TIMETXT,\n    # select station for cross\n    select_sta,STATXT,SRTXT,channel_regular,\n    # instrument parameters\n    rm_resp,respdir,\n    # some control parameters\n    input_fmt,time_norm,freq_norm,cc_method,auto_corr,cross_corr,\n    # pre-processing parameters\n    freqmin,freqmax,cc_len,cc_step,\n    # cross-correlation parameters\n    maxlag,substack,substack_method,\n    # criteria for data selection in time_norm\n    factor_clip_std,factor_mute,time_half_win,\n    # criteria for data selection in freq_norm\n    freq_half_win,\n)","category":"page"},{"location":"FFT_CC.html","page":"FFT_CC","title":"FFT_CC","text":"cores: Int64; total number of processes. cores number should be less than time-chunck number corresponded with \"select_time\".\nthreads: Int64; total number of threads. threads number should be less than N, N is the number of sta-pairs in each time-chunck.\nflag: true of false; print computing time for debugging purpose.\nprecompile: true or false; precompile or not.\nDATADIR: absolute dir where sac/mseed files are stored.\nCCFDIR: dir to store CC data.\nselect_time: true or false; whether specify time-chunck folder path.\nTIMETXT: time-chunck folder path.\nselect_sta: \"no\", \"sta\" or \"sr\"; \"sta\" cross-correlate only works for selected net-sta pairs base on STATXT; \"sr\" cross-correlate only works for specified source and receiver based on SRTXT.\nSTATXT: station infomation including: network-name and station-name.\nSRTXT: source and receiver information including: source-1 receiver-1 r-2 r-3...\\nsource-2 r-1 r-2 r-3...\nchannel_regular: channel(regular expression) of FFT-data; we will use 'glob(comp,FFTDIR)' to select data in 'FFTDIR' path.\nrm_resp: \"no\", \"xml\", \"RESP\", and \"PZs\"; select 'no' to not remove response and use \"xml\", \"RESP\", or \"PZs\" to remove response\nrespdir: path where resp files are located (required if rm_resp != \"no\")\ninput_fmt: \"sac\" or \"mseed\".\ntime_norm: \"no\", \"one_bit\", \"clip\", \"mute\", and \"smooth\"; \"no\" for no normalization in time domain.\nfreq_norm: \"no\", \"whiten\", \"coherence\", and \"deconvolution\"; \"no\" for no whitening, \"whiten\" for preserving only the phase component of the signal, \"coherence\" for smoothing absolute value of itself, and \"deconvolution\" for smoothing absolute value squared of itself.\ncc_method: \"CC\" and \"PCC\"; FFT transformation depends on different cross-correlation methods, \"CC\" for pure cross correlation, and \"PCC\" for Phase cross-correlation, see [Ventosa et al., 2019] for \"PCC\".\nauto_corr: true or false; perform auto-correlation or not.\ncross_corr: true or false; perform cross-correlation or not.\nfreqmin: Float64; for 'bandpass!' of raw-data and \"whiten\" parameter \nfreqmax: Float64; for 'bandpass!' of raw-data and \"whiten\" parameter\ncc_len: Int64; basic unit of data length for fft (sec)\ncc_step: Int64; overlapping between each cc_len (sec)\nmaxlag: Float64; lags of cross-correlation to save (sec).\nsubstack: true or false; sub-stack one time-chunck result of cross-correlation or not.\nsubstack_method: \"mean\", \"pws\", and \"robust\"; Vector{String}, you can choose one or more like [\"mean\",\"pws\",\"robust\"].\nfactor_clip_std: Float64; 'clip' (default value), clipping the data whose std exceeds the 'factorclipstd' to max/min mum value. || 中文: NoisePy与之对应的参数是 'maxoverstd', 其默认值设置为10, 时间序列超过10倍均方差std的时间序列，把振幅设置成上下线（最大值），限制在一个范围内.\nfactor_mute: Float64; 'mute' (default value), mute the data whose amp exceed 'factor_mute' times the median to 0, for suppressing high-amplitude signals. || 中文: 时间序列超过3倍于波形的平均中值(使用 hibert 变换), 把振幅设置成0. median of envelope of A to find outliers.\ntime_half_win: Int64; 'smooth' (default value), the length(sample points) of the window for running mean normalization, using smooth(A, halfwin) -> Smooth array A with half-window halfwin (defaults to 3). || 中文: 将一段时间序列分成许多小的时窗，对每一个小的时窗内的波形进行滑动平均，NoisePy对应的参数是 \"rma\", 其窗口长度的默认值 smoothN=10; 这也是bensen(2007)推荐的 timenorm 的方法.\nfreq_half_win: Int64; 'coherence' and 'deconvolution' (default value), the length(sample points) of the window for running mean normalization, using coherence!(F,20) or deconvolution!(F,20) -> Smooth array F with half-window halfwin (defaults to 20). || 中文: 频率域内的滑动平均，NoisePy对应的参数是 \"rma\", 其窗口长度的默认值 smoothN=10; NoisePy 中 timenorm and freqnorm 使用的参数是同一个，都是smooth_N","category":"page"},{"location":"references.html#Function-References","page":"References","title":"Function References","text":"","category":"section"},{"location":"references.html","page":"References","title":"References","text":"Modules = [NoiseCC]\nOrder = [:function, :type]","category":"page"},{"location":"references.html#NoiseCC.get_allfile_path_fft-Tuple{Bool, String, String, String, String, String}","page":"References","title":"NoiseCC.get_allfile_path_fft","text":"get_allfile_path(Parameters)\n\nget all files path –> files\n\nArguments\n\nParameters::Dict{String,Any}: ...\n\n\n\n\n\n","category":"method"},{"location":"FFT.html#FFT","page":"FFT","title":"FFT","text":"","category":"section"},{"location":"FFT.html","page":"FFT","title":"FFT","text":"Use Multi-process to accelerate FFT from raw-data (sac/mseed) to fft-data (FFTData).","category":"page"},{"location":"FFT.html","page":"FFT","title":"FFT","text":"############################################\n## 1.Set parameters \n############################################\ncores = 2\nflag = true\nprecompile = false\n... \n# Set your own parameters here...\n# Then run following code\n\n############################################\n# 2.run NoiseCC.FFT\n############################################\nusing Distributed\naddprocs(cores, topology = :master_worker)\n@everywhere using NoiseCC\n\nNoiseCC.FFT(\n    # processes parameters\n    cores,flag,precompile,\n    # absolute path parameters\n    FFTDIR,\n    # select sac/mseed files or not\n    select_file,\n    # when select_file = ture, 'TMPTXT' works\n    FILETXT,\n    # when select_file = false, 'DATADIR' 'select_sta' 'STATXT' and 'ncomp' work.\n    DATADIR,select_sta,STATXT,channel_regular,\n    # instrument parameters\n    rm_resp,respdir,\n    # some control parameters\n    input_fmt,time_norm,freq_norm,cc_method,\n    # pre-processing parameters\n    freqmin,freqmax,cc_len,cc_step,\n    # criteria for data selection in time_norm\n    factor_clip_std,factor_mute,time_half_win,\n    # criteria for data selection in freq_norm\n    freq_half_win,\n)","category":"page"},{"location":"FFT.html","page":"FFT","title":"FFT","text":"cores: total number of processes. cores number should be less than data-file number.\nflag: true of false; print computing time for debugging purpose.\nprecompile: true or false; precompile or not.\nFFTDIR: dir to store FFTData.\nselect_file: true of false; whether specify sac/mseed files path.\nFILETXT: sac/mseed files path.\nDATADIR: absolute dir where sac/mseed files are stored.\nselect_sta: \"no\" or \"sta\"; fft only works for selected net-sta base on STATXT.\nSTATXT: station infomation including: network-name and station-name.\nchannel_regular: channel(regular expression) of sac/mseed data; we will use 'glob(comp,DATADIR)' to select data in 'DATADIR' path.\nrm_resp: \"no\", \"xml\", \"RESP\", and \"PZs\"; select 'no' to not remove response and use \"xml\", \"RESP\", or \"PZs\" to remove response\nrespdir: path where resp files are located (required if rm_resp != \"no\")\ninput_fmt: \"sac\" or \"mseed\".\ntime_norm: \"no\", \"one_bit\", \"clip\", \"mute\", and \"smooth\"; \"no\" for no normalization in time domain.\nfreq_norm: \"no\", \"whiten\", \"coherence\", and \"deconvolution\"; \"no\" for no whitening, \"whiten\" for preserving only the phase component of the signal, \"coherence\" for smoothing absolute value of itself, and \"deconvolution\" for smoothing absolute value squared of itself.\ncc_method: \"CC\" and \"PCC\"; FFT transformation depends on different cross-correlation methods, \"CC\" for pure cross correlation, and \"PCC\" for Phase cross-correlation, see [Ventosa et al., 2019] for \"PCC\".\nfreqmin: Float64; for 'bandpass!' of raw-data and \"whiten\" parameter \nfreqmax: Float64; for 'bandpass!' of raw-data and \"whiten\" parameter\ncc_len: Int64; basic unit of data length for fft (sec)\ncc_step: Int64; overlapping between each cc_len (sec)\nfactor_clip_std: Float64; 'clip' (default value), clipping the data whose std exceeds the 'factorclipstd' to max/min mum value. || 中文: NoisePy与之对应的参数是 'maxoverstd', 其默认值设置为10, 时间序列超过10倍均方差std的时间序列，把振幅设置成上下线（最大值），限制在一个范围内.\nfactor_mute: Float64; 'mute' (default value), mute the data whose amp exceed 'factor_mute' times the median to 0, for suppressing high-amplitude signals. || 中文: 时间序列超过3倍于波形的平均中值(使用 hibert 变换), 把振幅设置成0. median of envelope of A to find outliers.\ntime_half_win: Int64; 'smooth' (default value), the length(sample points) of the window for running mean normalization, using smooth(A, halfwin) -> Smooth array A with half-window halfwin (defaults to 3). || 中文: 将一段时间序列分成许多小的时窗，对每一个小的时窗内的波形进行滑动平均，NoisePy对应的参数是 \"rma\", 其窗口长度的默认值 smoothN=10; 这也是bensen(2007)推荐的 timenorm 的方法.\nfreq_half_win: Int64; 'coherence' and 'deconvolution' (default value), the length(sample points) of the window for running mean normalization, using coherence!(F,20) or deconvolution!(F,20) -> Smooth array F with half-window halfwin (defaults to 20). || 中文: 频率域内的滑动平均，NoisePy对应的参数是 \"rma\", 其窗口长度的默认值 smoothN=10; NoisePy 中 timenorm and freqnorm 使用的参数是同一个，都是smooth_N","category":"page"},{"location":"STACK.html#STACK","page":"STACK","title":"STACK","text":"","category":"section"},{"location":"STACK.html","page":"STACK","title":"STACK","text":"Use Multi-process and Multi-thread to accelerate STACK with corr-data (CORRData).","category":"page"},{"location":"STACK.html","page":"STACK","title":"STACK","text":"############################################\n## 1.Set parameters \n############################################\ncores = 1\nthreads = 1\nflag = true\nprecompile = false\n...\n# Set your own parameters here...\n# Then run following code\n\n############################################\n# 2.run NoiseCC.STACK\n############################################\nusing Distributed\naddprocs(cores, topology = :master_worker)\n@everywhere using NoiseCC\n\nNoiseCC.STACK(\n    # process and threads parameters\n    cores,threads,flag,precompile,\n    # absolute path parameters\n    CCFDIR,STACKDIR,\n    # select time-chunck folder for stack\n    select_time,TIMETXT,\n    # select station for stack\n    select_sta,STATXT,SRTXT,channel_regular,comp,\n    # stack time interval\n    stack_all,time_interval,min_stack_chunck,\n    # stack method\n    substack_method,stack_method,auto_stack,cross_stack,\n    # filter \n    select_filter,freqmin,freqmax,corners,zerophase,\n    # criteria for corr-data selection\n    select_corr,median_high,median_low,\n)","category":"page"},{"location":"STACK.html","page":"STACK","title":"STACK","text":"cores: Int64; total number of processes. cores number should be less than \"NSRpairs\" number.\nthreads: Int64; total number of threads. threads number should be less than N piece of data when \"stackall=false\". please set to be \"1\" when \"stackall=true\"\nflag: true of false; print computing time for debugging purpose.\nprecompile: true or false; precompile or not.\nCCFDIR: dir where CC data are stored.\nSTACKDIR: dir to store STACK data.\nselect_time: true or false; whether specify time-chunck folder path.\nTIMETXT: time-chunck folder path.\nselect_sta: \"sta\" or \"sr\"; \"sta\" cross-correlate only works for selected net-sta pairs base on STATXT; \"sr\" cross-correlate only works for specified source and receiver based on SRTXT.\nSTATXT: station infomation including: network-name and station-name.\nSRTXT: source and receiver information including: source-1 receiver-1 r-2 r-3...\\nsource-2 r-1 r-2 r-3...\nchannel_regular: channel(regular expression) of CCF-data; we will use 'glob(comp,CCFDIR)' to select data in 'CCFDIR' path.\ncomp: component name corresponds to \"channelregular\"; different from NoiseCC.CC, we will use 'loadcorr(file,comp)' to read CCFData in 'CCFDIR' path, and \"comp\" is the group name of jld2 file(HDF5 in julia).\nstack_all: true or false; \"true\" for stack corr-data over all the selected time chunck.\ntime_interval: Int64 (sec); stack every \"timeinterval\" time, and it works when stackall = false || 中文: 每间隔多少秒叠加一次，要根据corr-data文件之间的时间间隔决定，例如corr-data文件是6个小时一个，那么time_interval=126060，表示每12个小时叠加一次，即需要2个corr文件叠加在一起。该参数工作时，先搜索指定台站对的全部时间的文件，找到其中最小的开始时间与最大的结束时间，利用range函数生成间隔为\"\"的叠加时窗，叠加每个时窗内的corr-data文件，叠加后的文件名称为该时窗的起始与终止时刻。\nmin_stack_chunck: Int64; each stack window needs at least \"minstackchunck\" to be able to stack, and it works for both stackall = true and stackall = false. || 中文:每个时间窗内叠加，至少需要 \"minstackchunck\" 个corr-data文件才会叠加,我们在读取corr-data文件后，会进行质量评估，剔除异常文件，导致corr-data文件数量减少。\nsubstack_method: \"mean\", \"pws\", and \"robust\"; substack method corresponding with NoiseCC.CC process, but can only choose one type.\nstack_method: \"mean\", \"pws\", and \"robust\"; stack method, but can only choose one type.\nauto_stack: true or false; perform auto-correlation or not.\ncross_stack: true or false; perform cross-correlation or not.\nselect_filter: true or false; whether filter before stacking.\nfreqmin: Float64;\nfreqmax: Float64;\ncorners: Int64; number of corners in Butterworth filter.\nzerophase: true or false; if true, apply filter twice for no phase shifting.\nselect_corr: true or false; whether to perform quality selecting based on median value of amp of corr-data.\nmedian_high: Float64; max median value (default value).\nmedian_low: Float64; min median value (default value).","category":"page"},{"location":"plot.html#Plot","page":"Plot","title":"Plot","text":"","category":"section"},{"location":"plot.html","page":"Plot","title":"Plot","text":"I'm used to using Matplotlib for visualization, HDF5 files can be read by h5py package.","category":"page"},{"location":"plot.html#S0_read_JLD2.py","page":"Plot","title":"S0_read_JLD2.py","text":"","category":"section"},{"location":"plot.html","page":"Plot","title":"Plot","text":"import numpy as np\nimport h5py\nimport os\nimport glob\n\n#%% 1. set parameter\nfile = \"../../data/BJ.081_BJ.084__2020_04_11_00_00_00T2021_04_13_00_00_00__all.jld2\"\nchan = \"NN\"\ndt = 0.005 #?\n\n#%% 2. read h5\n# open file\nf = h5py.File(file,'r')\n# read data\ndata = f[chan][\"data\"][0] \n# read parameters\nazi = f[chan][\"azi\"][()]\nbaz = f[chan][\"baz\"][()]\nmaxlag = f[chan][\"maxlag\"][()]\ncc_len = f[chan][\"cc_len\"][()]\ncc_step = f[chan][\"cc_step\"][()]\ncorr_type = f[chan][\"corr_type\"][()]\ncomp = f[chan][\"comp\"][()]\ndist = f[chan][\"dist\"][()]   # dist = f[chan][\"dist\"].value\nlat = f[chan][\"lat\"][()]\nlon = f[chan][\"lon\"][()]\nN_glob = f[chan][\"N_glob\"][()]\nN_read = f[chan][\"N_read\"][()]\nN_good = f[chan][\"N_good\"][()]\nname = f[chan][\"name\"][()][0].decode('utf-8')\n# close h5-file\nf.close()","category":"page"},{"location":"plot.html#S1_plot_waveform.py","page":"Plot","title":"S1_plot_waveform.py","text":"","category":"section"},{"location":"plot.html","page":"Plot","title":"Plot","text":"import obspy\nimport numpy as np\nimport h5py\nimport os\nimport glob\nimport matplotlib\nimport matplotlib.pyplot as plt\n\n#%% 1. set parameter\npath = \"/Users/yf/3.Project/2-1.Velocity_Change/STACK_das3/stack_mean\"\nnet_sta = \"BJ.106\"\nchan = \"NN\"\ndt = 0.005\n\n#%% 2. read h5\nfiles = sorted( glob.glob(os.path.join(path,net_sta+\"/*/*\")))\ndata_list = []\ndist_list = []\nname_list = []\nfor i in range(0,len(files)):\n    f = h5py.File(files[i],'r')\n    data = f[chan][\"data\"][0] \n    dist = f[chan][\"dist\"][()]\n    lat = f[chan][\"lat\"][()]\n    lon = f[chan][\"lon\"][()]\n    N_glob = f[chan][\"N_glob\"][()]\n    N_read = f[chan][\"N_read\"][()]\n    N_good = f[chan][\"N_good\"][()]\n    name = f[chan][\"name\"][()][0].decode('utf-8')\n    data_list.append(data)\n    name_list.append(name)\n    dist_list.append(dist)\n    f.close()\n\n#%% 3. plot\ntt = np.arange(-int(5),int(5),dt)\ndist=0\nfor i in range(0,len(files)):\n    dist += 1\n    data = data_list[i]/np.max(data_list[i],axis=0)\n    plt.plot(tt,data[3000:5000]+dist,'k',linewidth=0.8)\n\nplt.xlabel('Time (s)')\nplt.ylabel('Offset (m)')\n#plt.savefig(\"seisnoise2.pdf\", format='pdf', dpi=400)","category":"page"},{"location":"CC.html#CC","page":"CC","title":"CC","text":"","category":"section"},{"location":"CC.html","page":"CC","title":"CC","text":"Use Multi-process and Multi-thread to accelerate cross-correlation from fft-data (FFTData) to corr-data (CORRData).","category":"page"},{"location":"CC.html","page":"CC","title":"CC","text":"############################################\n## 1.Set parameters \n############################################\ncores = 2 \nthreads = 1\nflag = true\nprecompile = false\n... \n# Set your own parameters here...\n# Then run following code\n\n############################################\n# 2.run NoiseCC.CC\n############################################\nusing Distributed\naddprocs(cores, topology = :master_worker)\n@everywhere using NoiseCC\n\nNoiseCC.CC(\n    # processes and threads parameters\n    cores,threads,flag,precompile,\n    # absolute path parameters\n    FFTDIR,CCFDIR,\n    # select time-chunck folder for cross-correlate\n    select_time,TIMETXT,\n    # select station for cross\n    select_sta,STATXT,SRTXT,channel_regular,comp,\n    # some control parameters\n    cc_method,auto_corr,cross_corr,\n    # cross-correlation parameters\n    maxlag,substack,substack_method,\n)","category":"page"},{"location":"CC.html","page":"CC","title":"CC","text":"cores: Int64; total number of processes. cores number should be less than time-chunck number corresponded with \"select_time\".\nthreads: Int64; total number of threads. threads number should be less than N, N is the number of sta-pairs in each time-chunck.\nflag: true of false; print computing time for debugging purpose.\nprecompile: true or false; precompile or not.\nFFTDIR: dir where FFT data are stored.\nCCFDIR: dir to store CC data.\nselect_time: true or false; whether specify time-chunck folder path.\nTIMETXT: time-chunck folder path.\nselect_sta: \"no\", \"sta\" or \"sr\"; \"sta\" cross-correlate only works for selected net-sta pairs base on STATXT; \"sr\" cross-correlate only works for specified source and receiver based on SRTXT.\nSTATXT: station infomation including: network-name and station-name.\nSRTXT: source and receiver information including: source-1 receiver-1 r-2 r-3...\\nsource-2 r-1 r-2 r-3...\nchannel_regular: channel(regular expression) of FFT-data; we will use 'glob(comp,FFTDIR)' to select data in 'FFTDIR' path.\ncomp: component name corresponds to \"channelregular\"; we will use 'loadfft(file,comp)' to read FFTData in 'FFTDIR' path, and \"comp\" is the group name of jld2 file(HDF5 in julia).\ncc_method: \"CC\" and \"PCC\"; \"CC\" for pure cross correlation, or \"PCC\" for Phase cross-correlation, see [Ventosa et al., 2019].\nauto_corr: true or false; perform auto-correlation or not.\ncross_corr: true or false; perform cross-correlation or not.\nmaxlag: Float64; lags of cross-correlation to save (sec).\nsubstack: true or false; sub-stack one time-chunck result of cross-correlation or not.\nsubstack_method: \"mean\", \"pws\", and \"robust\"; Vector{String}, you can choose one or more like [\"mean\",\"pws\",\"robust\"].","category":"page"},{"location":"Framwork.html#Prepare-database","page":"Framwork","title":"Prepare database","text":"","category":"section"},{"location":"Framwork.html","page":"Framwork","title":"Framwork","text":"Before using NoiseCC.jl, you should prepare the database. Make sure the name of the file contains the exact network, station and component information. The following figure shows the format of other paths in the database. ","category":"page"},{"location":"Framwork.html","page":"Framwork","title":"Framwork","text":"Format: Year | Yearday | Time-Chunck | Data(sac/mseed)","category":"page"},{"location":"Framwork.html","page":"Framwork","title":"Framwork","text":"(Image: database_format)","category":"page"},{"location":"Framwork.html#Processing-flow","page":"Framwork","title":"Processing flow","text":"","category":"section"},{"location":"Framwork.html","page":"Framwork","title":"Framwork","text":"(Image: Process_flow)","category":"page"},{"location":"Framwork.html#CC-flow","page":"Framwork","title":"CC flow","text":"","category":"section"},{"location":"Framwork.html","page":"Framwork","title":"Framwork","text":"(Image: Process_flow)","category":"page"},{"location":"index.html#NoiseCC.jl:-A-framwork-for-large-ambient-noise-cross-correlation-(CC)-based-on-SeisNoise.jl","page":"Home","title":"NoiseCC.jl: A framwork for large ambient noise cross-correlation (CC) based on SeisNoise.jl","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"NoiseCC.jl use Multi-process and Multi-thread to accelerate noise cross-correlation (CC) for large dataset. We designed a generic CC framework and wrap SeisNoise module into a Julia interface for the convenience of users, which is suitable for large-scale Dense-Stations, Large-N and DAS data sets. It is convenient for the subsequent study of noise-imaging and noise-monitor.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"note: Note\nAll the core code for cross-correlation comes from SeisNoise.jl and SeisIO Documentation, and this package is an example of SeisNoise.jl application. Please read through the SeisNoise.jl and SeisIO Documentation to get familiar with seismic data processing and parameters setted for cross-correlation in Julia.","category":"page"},{"location":"index.html#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Install SeisNoise firstly. Use the Julia package manager (Press ] to enter pkg):","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia>]\n(@v1.6) pkg> add SeisNoise","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Then install NoiseCC:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"(@v1.6) pkg> add NoiseCC","category":"page"},{"location":"index.html#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Wrap SeisNoise module into a Julia interface for easy using.\nMulti-process and Multi-thread to accelerate noise cross-correlation.\nA standardized seismic data processing framework.","category":"page"},{"location":"index.html#Author-Saying","page":"Home","title":"Author Saying","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Why I wrap SeisNoise into a API in julia?","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Seismic cross-correlation calculation is the most common task in seismology. It often requires designing different program solutions based on data in different scenarios, SeisNoise was originally designed to facilitate users to write their own programs for different problems. But writing parallel programs for large data is often repetitive and unfriendly to beginners, which means learning a new language and writing your own programs.\nMy idea was to provide a common programming interface that would work with most scenarios, which would save the time to rewrite code. So, NoiseCC can be seen as an application case of SeisNoise, which meets my requirement of cross-correlation calculation in most scenarios. I want to make this part of the code available to anyone who needs it.","category":"page"}]
}
