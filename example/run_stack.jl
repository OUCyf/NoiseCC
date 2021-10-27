# Greetings! from path ~/.julia/config/startup.jl
############################################
## 1.Set parameters 
############################################
# process and threads parameters
cores = 1                                                           # Int64; total number of processes. cores number should be less than "N_SR_pairs" number.
threads = 1                                                         # Int64; total number of threads. threads number should be less than N piece of data when "stack_all=false". please set to be "1" when "stack_all=true"
flag = true                                                         # true of false; print computing time for debugging purpose.
precompile = false                                                  # true or false; precompile or not.

# absolute path parameters
CCFDIR = expanduser("./CORR_fft")                              # dir where CC data are stored.
STACKDIR = expanduser("./STACK_d")                               # dir to store STACK data.

# select time-chunck folder for stack
select_time = false                                                  # true or false; whether specify time-chunck folder path.
TIMETXT = expanduser("./time_stack.txt")                            # time-chunck folder path.

# select station for stack
select_sta = "sr"                                                   # "sta" or "sr"; "sta" cross-correlate only works for selected net-sta pairs base on STATXT; "sr" cross-correlate only works for specified source and receiver based on SRTXT.
STATXT = expanduser("./sta_stack.txt")                              # station infomation including: network-name and station-name.
SRTXT = expanduser("./sr_stack.txt")                                # source and receiver information including: source-1 receiver-1 r-2 r-3...\nsource-2 r-1 r-2 r-3...
channel_regular = "*.CHN.*.CHN.*"                                   # channel(regular expression) of CCF-data; we will use 'glob(comp,CCFDIR)' to select data in 'CCFDIR' path.
comp = "NN"                                                         # component name corresponds to "channel_regular"; different from NoiseCC.CC, we will use 'load_corr(file,comp)' to read CCFData in 'CCFDIR' path, and "comp" is the group name of jld2 file(HDF5 in julia).

# stack time interval
stack_all = true                                                    # true or false; "true" for stack corr-data over all the selected time chunck.
time_interval = 1 * 24 * 60 * 60                                    # Int64 (sec); stack every "time_interval" time, and it works when stack_all = false || 中文: 每间隔多少秒叠加一次，要根据corr-data文件之间的时间间隔决定，例如corr-data文件是6个小时一个，那么time_interval=12*60*60，表示每12个小时叠加一次，即需要2个corr文件叠加在一起。该参数工作时，先搜索指定台站对的全部时间的文件，找到其中最小的开始时间与最大的结束时间，利用range函数生成间隔为""的叠加时窗，叠加每个时窗内的corr-data文件，叠加后的文件名称为该时窗的起始与终止时刻。
min_stack_chunck = 1                                                # Int64; each stack window needs at least "min_stack_chunck" to be able to stack, and it works for both stack_all = true and stack_all = false. || 中文:每个时间窗内叠加，至少需要 "min_stack_chunck" 个corr-data文件才会叠加,我们在读取corr-data文件后，会进行质量评估，剔除异常文件，导致corr-data文件数量减少。

# stack method
substack_method = "pws"                                             # "mean", "pws", and "robust"; substack method corresponding with NoiseCC.CC process, but can only choose one type.
stack_method = "mean"                                               # "mean", "pws", and "robust"; stack method, but can only choose one type.
auto_stack = true                                                   # true or false; perform auto-correlation or not.
cross_stack = true                                                  # true or false; perform cross-correlation or not.

# filter 
select_filter = true                                                # true or false; whether filter before stacking.
freqmin = 0.5                                                       # Float64;
freqmax = 20.0                                                      # Float64;
corners = 4                                                         # Int64; number of corners in Butterworth filter.
zerophase = false                                                   # true or false; if true, apply filter twice for no phase shifting.

# criteria for corr-data selection
select_corr = true                                                  # true or false; whether to perform quality selecting based on median value of amp of corr-data.
median_high = 10.0                                                  # Float64; max median value (default value).
median_low = 0.0                                                    # Float64; min median value (default value).



############################################
# 2.run NoiseCC.FFT_CC
############################################
# a.using and import
# using Distributed
# NoiseCC_path = "/Users/yf/3.Project/2-1.Velocity_Change/SeisNoise_Code_YF/code_pmap_thread/"
# push!(LOAD_PATH, NoiseCC_path)
# @everywhere using NoiseCC # or @everywhere import NoiseCC

# b.include
using Distributed
addprocs(cores, topology = :master_worker, exeflags = "-t $(threads)")
@everywhere include(
    "/Users/yf/3.Project/2-1.Velocity_Change/SeisNoise_Code_YF/code_pmap_thread/NoiseCC.jl",
)


NoiseCC.STACK(
    cores,
    threads,
    flag,
    precompile,
    CCFDIR,
    STACKDIR,
    select_time,
    TIMETXT,
    select_sta,
    STATXT,
    SRTXT,
    channel_regular,
    comp,
    stack_all,
    time_interval,
    min_stack_chunck,
    substack_method,
    stack_method,
    auto_stack,
    cross_stack,
    select_filter,
    freqmin,
    freqmax,
    corners,
    zerophase,
    select_corr,
    median_high,
    median_low,
)
