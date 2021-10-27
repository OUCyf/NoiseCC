# Greetings! from path ~/.julia/config/startup.jl

# """
#     1. Generate tmp file: find $(pwd) -name \*.sac > ./tmp.txt
#     2. Run S0A_FFT.jl script
# """

############################################
## 1.Set parameters 
############################################
# processes parameters
cores = 2                                               # total number of processes. cores number should be less than data-file number.
flag = true                                             # true of false; print computing time for debugging purpose.
precompile = false                                      # true or false; precompile or not.

# absolute path parameters
FFTDIR = expanduser("./FFT")                            # dir to store FFTData.

# select sac/mseed files or not
select_file = false                                     # true of false; whether specify sac/mseed files path.

# when select_file = ture, 'TMPTXT' works
FILETXT = expanduser("./file.txt")                      # sac/mseed files path.

# when select_file = false, 'DATADIR' 'select_sta' 'STATXT' and 'ncomp' work.
DATADIR = expanduser("./DATA2")                         # absolute dir where sac/mseed files are stored.
select_sta = "sta"                                      # "no" or "sta"; fft only works for selected net-sta base on STATXT.
STATXT = expanduser("./station_fft.txt")                # station infomation including: network-name and station-name.
channel_regular = "*CHN*"                               # channel(regular expression) of sac/mseed data; we will use 'glob(comp,DATADIR)' to select data in 'DATADIR' path.

# instrument parameters (The functionality has not been added!! please use [rm_resp = "no"])
rm_resp = "no"                                          # "no", "xml", "RESP", and "PZs"; select 'no' to not remove response and use "xml", "RESP", or "PZs" to remove response
respdir = expanduser("./XX.PZs")                        # path where resp files are located (required if rm_resp != "no")

# some control parameters
input_fmt = "sac"                                       # "sac" or "mseed".
time_norm = "one_bit"                                   # "no", "one_bit", "clip", "mute", and "smooth"; "no" for no normalization in time domain.
freq_norm = "whiten"                                    # "no", "whiten", "coherence", and "deconvolution"; "no" for no whitening, "whiten" for preserving only the phase component of the signal, "coherence" for smoothing absolute value of itself, and "deconvolution" for smoothing absolute value squared of itself.
cc_method = "CC"                                        # "CC" and "PCC"; FFT transformation depends on different cross-correlation methods, "CC" for pure cross correlation, and "PCC" for Phase cross-correlation, see [Ventosa et al., 2019] for "PCC".

# pre-processing parameters
freqmin = 0.5                                           # Float64; for 'bandpass!' of raw-data and "whiten" parameter 
freqmax = 20.0                                          # Float64; for 'bandpass!' of raw-data and "whiten" parameter
cc_len = 60 * 2                                         # Int64; basic unit of data length for fft (sec)
cc_step = 60 * 1                                        # Int64; overlapping between each cc_len (sec)

# criteria for data selection in time_norm
factor_clip_std = 10.0                                  # Float64; 'clip' (default value), clipping the data whose std exceeds the 'factor_clip_std' to max/min mum value. || 中文: NoisePy与之对应的参数是 'max_over_std', 其默认值设置为10, 时间序列超过10倍均方差std的时间序列，把振幅设置成上下线（最大值），限制在一个范围内.
factor_mute = 3.0                                       # Float64; 'mute' (default value), mute the data whose amp exceed 'factor_mute' times the median to 0, for suppressing high-amplitude signals. || 中文: 时间序列超过3倍于波形的平均中值(使用 hibert 变换), 把振幅设置成0. median of envelope of `A` to find outliers.
time_half_win = 3                                       # Int64; 'smooth' (default value), the length(sample points) of the window for running mean normalization, using smooth(A, half_win) -> Smooth array A with half-window half_win (defaults to 3). || 中文: 将一段时间序列分成许多小的时窗，对每一个小的时窗内的波形进行滑动平均，NoisePy对应的参数是 "rma", 其窗口长度的默认值 smooth_N=10; 这也是bensen(2007)推荐的 time_norm 的方法.

# criteria for data selection in freq_norm
freq_half_win = 20                                      # Int64; 'coherence' and 'deconvolution' (default value), the length(sample points) of the window for running mean normalization, using coherence!(F,20) or deconvolution!(F,20) -> Smooth array F with half-window half_win (defaults to 20). || 中文: 频率域内的滑动平均，NoisePy对应的参数是 "rma", 其窗口长度的默认值 smooth_N=10; NoisePy 中 time_norm and freq_norm 使用的参数是同一个，都是smooth_N


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
addprocs(cores, topology = :master_worker)
@everywhere include("/Users/yf/3.Project/2-1.Velocity_Change/SeisNoise_Code_YF/code_pmap_thread/NoiseCC.jl")


NoiseCC.FFT(
    cores,
    flag,
    precompile,
    FFTDIR,
    select_file,
    FILETXT,
    DATADIR,
    select_sta,
    STATXT,
    input_fmt,
    time_norm,
    freq_norm,
    cc_method,
    channel_regular,
    rm_resp,
    respdir,
    freqmin,
    freqmax,
    cc_len,
    cc_step,
    factor_clip_std,
    factor_mute,
    time_half_win,
    freq_half_win,
)
