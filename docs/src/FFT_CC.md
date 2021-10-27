# FFT_CC
Use Multi-process and Multi-thread to accelerate FFT and CC from raw-data (sac/mseed) to corr-data (CORRData).


```julia
############################################
## 1.Set parameters 
############################################
cores = 1
threads = 1
flag = true
precompile = false
...
# Set your own parameters here...
# Then run following code

############################################
# 2.run NoiseCC.FFT
############################################
using Distributed
addprocs(cores, topology = :master_worker)
@everywhere using NoiseCC

NoiseCC.FFT_CC(
    # processes and threads parameters
    cores,threads,flag,precompile,
    # absolute path parameters
    DATADIR,CCFDIR,
    # select time-chunck folder for cross-correlate
    select_time,TIMETXT,
    # select station for cross
    select_sta,STATXT,SRTXT,channel_regular,
    # instrument parameters
    rm_resp,respdir,
    # some control parameters
    input_fmt,time_norm,freq_norm,cc_method,auto_corr,cross_corr,
    # pre-processing parameters
    freqmin,freqmax,cc_len,cc_step,
    # cross-correlation parameters
    maxlag,substack,substack_method,
    # criteria for data selection in time_norm
    factor_clip_std,factor_mute,time_half_win,
    # criteria for data selection in freq_norm
    freq_half_win,
)
```

- `cores`: Int64; total number of processes. cores number should be less than time-chunck number corresponded with "select_time".
- `threads`: Int64; total number of threads. threads number should be less than N, N is the number of sta-pairs in each time-chunck.
- `flag`: true of false; print computing time for debugging purpose.
- `precompile`: true or false; precompile or not.
- `DATADIR`: absolute dir where sac/mseed files are stored.
- `CCFDIR`: dir to store CC data.
- `select_time`: true or false; whether specify time-chunck folder path.
- `TIMETXT`: time-chunck folder path.
- `select_sta`: "no", "sta" or "sr"; "sta" cross-correlate only works for selected net-sta pairs base on STATXT; "sr" cross-correlate only works for specified source and receiver based on SRTXT.
- `STATXT`: station infomation including: network-name and station-name.
- `SRTXT`: source and receiver information including: source-1 receiver-1 r-2 r-3...\nsource-2 r-1 r-2 r-3...
- `channel_regular`: channel(regular expression) of FFT-data; we will use 'glob(comp,FFTDIR)' to select data in 'FFTDIR' path.
- `rm_resp`: "no", "xml", "RESP", and "PZs"; select 'no' to not remove response and use "xml", "RESP", or "PZs" to remove response
- `respdir`: path where resp files are located (required if rm_resp != "no")
- `input_fmt`: "sac" or "mseed".
- `time_norm`: "no", "one_bit", "clip", "mute", and "smooth"; "no" for no normalization in time domain.
- `freq_norm`: "no", "whiten", "coherence", and "deconvolution"; "no" for no whitening, "whiten" for preserving only the phase component of the signal, "coherence" for smoothing absolute value of itself, and "deconvolution" for smoothing absolute value squared of itself.
- `cc_method`: "CC" and "PCC"; FFT transformation depends on different cross-correlation methods, "CC" for pure cross correlation, and "PCC" for Phase cross-correlation, see [Ventosa et al., 2019] for "PCC".
- `auto_corr`: true or false; perform auto-correlation or not.
- `cross_corr`: true or false; perform cross-correlation or not.
- `freqmin`: Float64; for 'bandpass!' of raw-data and "whiten" parameter 
- `freqmax`: Float64; for 'bandpass!' of raw-data and "whiten" parameter
- `cc_len`: Int64; basic unit of data length for fft (sec)
- `cc_step`: Int64; overlapping between each cc_len (sec)
- `maxlag`: Float64; lags of cross-correlation to save (sec).
- `substack`: true or false; sub-stack one time-chunck result of cross-correlation or not.
- `substack_method`: "mean", "pws", and "robust"; Vector{String}, you can choose one or more like ["mean","pws","robust"].
- `factor_clip_std`: Float64; 'clip' (default value), clipping the data whose std exceeds the 'factor_clip_std' to max/min mum value. || 中文: NoisePy与之对应的参数是 'max_over_std', 其默认值设置为10, 时间序列超过10倍均方差std的时间序列，把振幅设置成上下线（最大值），限制在一个范围内.
- `factor_mute`: Float64; 'mute' (default value), mute the data whose amp exceed 'factor_mute' times the median to 0, for suppressing high-amplitude signals. || 中文: 时间序列超过3倍于波形的平均中值(使用 hibert 变换), 把振幅设置成0. median of envelope of `A` to find outliers.
- `time_half_win`: Int64; 'smooth' (default value), the length(sample points) of the window for running mean normalization, using smooth(A, half_win) -> Smooth array A with half-window half_win (defaults to 3). || 中文: 将一段时间序列分成许多小的时窗，对每一个小的时窗内的波形进行滑动平均，NoisePy对应的参数是 "rma", 其窗口长度的默认值 smooth_N=10; 这也是bensen(2007)推荐的 time_norm 的方法.
- `freq_half_win`: Int64; 'coherence' and 'deconvolution' (default value), the length(sample points) of the window for running mean normalization, using coherence!(F,20) or deconvolution!(F,20) -> Smooth array F with half-window half_win (defaults to 20). || 中文: 频率域内的滑动平均，NoisePy对应的参数是 "rma", 其窗口长度的默认值 smooth_N=10; NoisePy 中 time_norm and freq_norm 使用的参数是同一个，都是smooth_N
