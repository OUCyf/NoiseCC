# Greetings! from path ~/.julia/config/startup.jl

############################################
## 1.Set parameters 
############################################
# processes and threads parameters
cores = 2                                                           # Int64; total number of processes. cores number should be less than time-chunck number corresponded with "select_time".
threads = 1                                                         # Int64; total number of threads. threads number should be less than N, N is the number of sta-pairs in each time-chunck.
flag = true                                                         # true of false; print computing time for debugging purpose.
precompile = false                                                  # true or false; precompile or not.

# absolute path parameters
FFTDIR = expanduser("/Users/yf/3.Project/2-1.Velocity_Change/FFT")  # dir where FFT data are stored.
CCFDIR = expanduser("./CORR_fft")                                   # dir to store CC data.

# select time-chunck folder for cross-correlate
select_time = false                                                 # true or false; whether specify time-chunck folder path.
TIMETXT = expanduser("./time_cc.txt")                               # time-chunck folder path.

# select station for cross
select_sta = "sr"                                                   # "no", "sta" or "sr"; "sta" cross-correlate only works for selected net-sta pairs base on STATXT; "sr" cross-correlate only works for specified source and receiver based on SRTXT.
STATXT = expanduser("./sta_cc.txt")                                 # station infomation including: network-name and station-name.
SRTXT = expanduser("./sr_cc.txt")                                   # source and receiver information including: source-1 receiver-1 r-2 r-3...\nsource-2 r-1 r-2 r-3...
channel_regular = "*CHN.*"                                          # channel(regular expression) of FFT-data; we will use 'glob(comp,FFTDIR)' to select data in 'FFTDIR' path.
comp = "CHN"                                                        # component name corresponds to "channel_regular"; we will use 'load_fft(file,comp)' to read FFTData in 'FFTDIR' path, and "comp" is the group name of jld2 file(HDF5 in julia).

# some control parameters
cc_method = "CC"                                                    # "CC" and "PCC"; "CC" for pure cross correlation, or "PCC" for Phase cross-correlation, see [Ventosa et al., 2019].
auto_corr = true                                                    # true or false; perform auto-correlation or not.
cross_corr = true                                                   # true or false; perform cross-correlation or not.

# cross-correlation parameters
maxlag = 20.0                                                       # Float64; lags of cross-correlation to save (sec).
substack = true                                                     # true or false; sub-stack one time-chunck result of cross-correlation or not.
substack_method = ["mean","pws"]                                           # "mean", "pws", and "robust"; Vector{String}, you can choose one or more like ["mean","pws","robust"].



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


NoiseCC.CC(
    cores,
    threads,
    flag,
    precompile,
    FFTDIR,
    CCFDIR,
    select_time,
    TIMETXT,
    select_sta,
    STATXT,
    SRTXT,
    channel_regular,
    comp,
    cc_method,
    auto_corr,
    cross_corr,
    maxlag,
    substack,
    substack_method,
)
