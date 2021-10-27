# CC
Use Multi-process and Multi-thread to accelerate cross-correlation from fft-data (FFTData) to corr-data (CORRData).


```julia
############################################
## 1.Set parameters 
############################################
cores = 2 
threads = 1
flag = true
precompile = false
... 
# Set your own parameters here...
# Then run following code

############################################
# 2.run NoiseCC.CC
############################################
using Distributed
addprocs(cores, topology = :master_worker)
@everywhere using NoiseCC

NoiseCC.CC(
    # processes and threads parameters
    cores,threads,flag,precompile,
    # absolute path parameters
    FFTDIR,CCFDIR,
    # select time-chunck folder for cross-correlate
    select_time,TIMETXT,
    # select station for cross
    select_sta,STATXT,SRTXT,channel_regular,comp,
    # some control parameters
    cc_method,auto_corr,cross_corr,
    # cross-correlation parameters
    maxlag,substack,substack_method,
)
```

- `cores`: Int64; total number of processes. cores number should be less than time-chunck number corresponded with "select_time".
- `threads`: Int64; total number of threads. threads number should be less than N, N is the number of sta-pairs in each time-chunck.
- `flag`: true of false; print computing time for debugging purpose.
- `precompile`: true or false; precompile or not.
- `FFTDIR`: dir where FFT data are stored.
- `CCFDIR`: dir to store CC data.
- `select_time`: true or false; whether specify time-chunck folder path.
- `TIMETXT`: time-chunck folder path.
- `select_sta`: "no", "sta" or "sr"; "sta" cross-correlate only works for selected net-sta pairs base on STATXT; "sr" cross-correlate only works for specified source and receiver based on SRTXT.
- `STATXT`: station infomation including: network-name and station-name.
- `SRTXT`: source and receiver information including: source-1 receiver-1 r-2 r-3...\nsource-2 r-1 r-2 r-3...
- `channel_regular`: channel(regular expression) of FFT-data; we will use 'glob(comp,FFTDIR)' to select data in 'FFTDIR' path.
- `comp`: component name corresponds to "channel_regular"; we will use 'load_fft(file,comp)' to read FFTData in 'FFTDIR' path, and "comp" is the group name of jld2 file(HDF5 in julia).
- `cc_method`: "CC" and "PCC"; "CC" for pure cross correlation, or "PCC" for Phase cross-correlation, see [Ventosa et al., 2019].
- `auto_corr`: true or false; perform auto-correlation or not.
- `cross_corr`: true or false; perform cross-correlation or not.

- `maxlag`: Float64; lags of cross-correlation to save (sec).
- `substack`: true or false; sub-stack one time-chunck result of cross-correlation or not.
- `substack_method`: "mean", "pws", and "robust"; Vector{String}, you can choose one or more like ["mean","pws","robust"].
