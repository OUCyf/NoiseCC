############################################
## 1.Set parameters 
############################################
begin
    # processes and threads parameters
    cores = 1                                                           # Int64; total number of processes. cores number should be less than time-chunck number corresponded with "select_time".
    threads = 1                                                         # Int64; total number of threads. threads number should be less than N, N is the number of sta-pairs in each time-chunck.
    flag = true                                                         # true of false; print computing time for debugging purpose.
    precompile = false                                                   # true or false; precompile or not.

    # absolute path parameters
    DATADIR = expanduser("/Users/yf/3.Project/2-1.Velocity_Change/DATA2")                                     # absolute dir where sac/mseed files are stored.
    CCFDIR = expanduser("/Users/yf/3.Project/2-1.Velocity_Change/CORR_das_module1")                              # dir to store CC data.

    # select time-chunck folder for cross-correlate
    select_time = true                                                 # true or false; whether specify time-chunck folder path.
    TIMETXT = expanduser("/Users/yf/3.Project/2-1.Velocity_Change/time_cc.txt")                               # time-chunck folder path.

    # select station for cross
    select_sta = "sr"                                                   # "no", "sta" or "sr"; "sta" cross-correlate only works for selected net-sta pairs base on STATXT; "sr" cross-correlate only works for specified source and receiver based on SRTXT.
    STATXT = expanduser("/Users/yf/3.Project/2-1.Velocity_Change/sta_cc.txt")                                 # station infomation including: network-name and station-name.
    SRTXT = expanduser("/Users/yf/3.Project/2-1.Velocity_Change/sr_cc.txt")                                   # source and receiver information including: source-1 receiver-1 r-2 r-3...\nsource-2 r-1 r-2 r-3...
    channel_regular = "*CHN.*"                                          # channel(regular expression) of FFT-data; we will use 'glob(comp,FFTDIR)' to select data in 'FFTDIR' path.

    # instrument parameters (The functionality has not been added!! please use [rm_resp = "no"])
    rm_resp = "no"                                                      # "no", "xml", "RESP", and "PZs"; select 'no' to not remove response and use "xml", "RESP", or "PZs" to remove response
    respdir = expanduser("./XX.PZs")                                    # path where resp files are located (required if rm_resp != "no")

    # some control parameters
    input_fmt = "sac"                                                   # "sac" or "mseed".
    time_norm = "one_bit"                                               # "no", "one_bit", "clip", "mute", and "smooth"; "no" for no normalization in time domain.
    freq_norm = "whiten"                                                # "no", "whiten", "coherence", and "deconvolution"; "no" for no whitening, "whiten" for preserving only the phase component of the signal, "coherence" for smoothing absolute value of itself, and "deconvolution" for smoothing absolute value squared of itself.
    cc_method = "CC"                                                    # "CC" and "PCC"; FFT transformation depends on different cross-correlation methods, "CC" for pure cross correlation, and "PCC" for Phase cross-correlation, see [Ventosa et al., 2019] for "PCC".
    auto_corr = true                                                    # true or false; perform auto-correlation or not.
    cross_corr = true                                                   # true or false; perform cross-correlation or not.

    # pre-processing parameters
    freqmin = 0.5                                                       # Float64; for 'bandpass!' of raw-data and "whiten" parameter 
    freqmax = 20.0                                                      # Float64; for 'bandpass!' of raw-data and "whiten" parameter
    cc_len = 60 * 2                                                     # Int64; basic unit of data length for fft (sec)
    cc_step = 60 * 1                                                    # Int64; overlapping between each cc_len (sec)

    # cross-correlation parameters
    maxlag = 20.0                                                       # Float64; lags of cross-correlation to save (sec).
    substack = true                                                     # true or false; sub-stack one time-chunck result of cross-correlation or not.
    substack_method = ["mean","pws"]                                           # "mean", "pws", and "robust"; Vector{String}, you can choose one or more like ["mean","pws","robust"].

    # criteria for data selection in time_norm
    factor_clip_std = 10.0                                              # Float64; 'clip' (default value), clipping the data whose std exceeds the 'factor_clip_std' to max/min mum value. || 中文: NoisePy与之对应的参数是 'max_over_std', 其默认值设置为10, 时间序列超过10倍均方差std的时间序列，把振幅设置成上下线（最大值），限制在一个范围内.
    factor_mute = 3.0                                                   # Float64; 'mute' (default value), mute the data whose amp exceed 'factor_mute' times the median to 0, for suppressing high-amplitude signals. || 中文: 时间序列超过3倍于波形的平均中值(使用 hibert 变换), 把振幅设置成0. median of envelope of `A` to find outliers.
    time_half_win = 3                                                   # Int64; 'smooth' (default value), the length(sample points) of the window for running mean normalization, using smooth(A, half_win) -> Smooth array A with half-window half_win (defaults to 3). || 中文: 将一段时间序列分成许多小的时窗，对每一个小的时窗内的波形进行滑动平均，NoisePy对应的参数是 "rma", 其窗口长度的默认值 smooth_N=10; 这也是bensen(2007)推荐的 time_norm 的方法.

    # criteria for data selection in freq_norm
    freq_half_win = 20                                      # Int64; 'coherence' and 'deconvolution' (default value), the length(sample points) of the window for running mean normalization, using coherence!(F,20) or deconvolution!(F,20) -> Smooth array F with half-window half_win (defaults to 20). || 中文: 频率域内的滑动平均，NoisePy对应的参数是 "rma", 其窗口长度的默认值 smooth_N=10; NoisePy 中 time_norm and freq_norm 使用的参数是同一个，都是smooth_N
end



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
@everywhere include("/Users/yf/3.Project/2-1.Velocity_Change/SeisNoise_Code_YF/code_pmap_thread/NoiseCC.jl")

NoiseCC.FFT_CC(
    cores,
    threads,
    flag,
    precompile,
    DATADIR,
    CCFDIR,
    select_time,
    TIMETXT,
    select_sta,
    STATXT,
    SRTXT,
    channel_regular,
    rm_resp,
    respdir,
    input_fmt,
    time_norm,
    freq_norm,
    cc_method,
    auto_corr,
    cross_corr,
    freqmin,
    freqmax,
    cc_len,
    cc_step,
    maxlag,
    substack,
    substack_method,
    factor_clip_std,
    factor_mute,
    time_half_win,
    freq_half_win,
)
