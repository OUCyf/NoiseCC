# STACK
Use Multi-process and Multi-thread to accelerate STACK with corr-data (CORRData).


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
# 2.run NoiseCC.STACK
############################################
using Distributed
addprocs(cores, topology = :master_worker)
@everywhere using NoiseCC

NoiseCC.STACK(
    # process and threads parameters
    cores,threads,flag,precompile,
    # absolute path parameters
    CCFDIR,STACKDIR,
    # select time-chunck folder for stack
    select_time,TIMETXT,
    # select station for stack
    select_sta,STATXT,SRTXT,channel_regular,comp,
    # stack time interval
    stack_all,time_interval,min_stack_chunck,
    # stack method
    substack_method,stack_method,auto_stack,cross_stack,
    # filter 
    select_filter,freqmin,freqmax,corners,zerophase,
    # criteria for corr-data selection
    select_corr,median_high,median_low,
)
```

- `cores`: Int64; total number of processes. cores number should be less than "N_SR_pairs" number.
- `threads`: Int64; total number of threads. threads number should be less than N piece of data when "stack_all=false". please set to be "1" when "stack_all=true"
- `flag`: true of false; print computing time for debugging purpose.
- `precompile`: true or false; precompile or not.
- `CCFDIR`: dir where CC data are stored.
- `STACKDIR`: dir to store STACK data.
- `select_time`: true or false; whether specify time-chunck folder path.
- `TIMETXT`: time-chunck folder path.
- `select_sta`: "sta" or "sr"; "sta" cross-correlate only works for selected net-sta pairs base on STATXT; "sr" cross-correlate only works for specified source and receiver based on SRTXT.
- `STATXT`: station infomation including: network-name and station-name.
- `SRTXT`: source and receiver information including: source-1 receiver-1 r-2 r-3...\nsource-2 r-1 r-2 r-3...
- `channel_regular`: channel(regular expression) of CCF-data; we will use 'glob(comp,CCFDIR)' to select data in 'CCFDIR' path.
- `comp`: component name corresponds to "channel_regular"; different from NoiseCC.CC, we will use 'load_corr(file,comp)' to read CCFData in 'CCFDIR' path, and "comp" is the group name of jld2 file(HDF5 in julia).
- `stack_all`: true or false; "true" for stack corr-data over all the selected time chunck.
- `time_interval`: Int64 (sec); stack every "time_interval" time, and it works when stack_all = false || 中文: 每间隔多少秒叠加一次，要根据corr-data文件之间的时间间隔决定，例如corr-data文件是6个小时一个，那么time_interval=12*60*60，表示每12个小时叠加一次，即需要2个corr文件叠加在一起。该参数工作时，先搜索指定台站对的全部时间的文件，找到其中最小的开始时间与最大的结束时间，利用range函数生成间隔为""的叠加时窗，叠加每个时窗内的corr-data文件，叠加后的文件名称为该时窗的起始与终止时刻。
- `min_stack_chunck`: Int64; each stack window needs at least "min_stack_chunck" to be able to stack, and it works for both stack_all = true and stack_all = false. || 中文:每个时间窗内叠加，至少需要 "min_stack_chunck" 个corr-data文件才会叠加,我们在读取corr-data文件后，会进行质量评估，剔除异常文件，导致corr-data文件数量减少。
- `substack_method`: "mean", "pws", and "robust"; substack method corresponding with NoiseCC.CC process, but can only choose one type.
- `stack_method`: "mean", "pws", and "robust"; stack method, but can only choose one type.
- `auto_stack`: true or false; perform auto-correlation or not.
- `cross_stack`: true or false; perform cross-correlation or not.
- `select_filter`: true or false; whether filter before stacking.
- `freqmin`: Float64;
- `freqmax`: Float64;
- `corners`: Int64; number of corners in Butterworth filter.
- `zerophase`: true or false; if true, apply filter twice for no phase shifting.
- `select_corr`: true or false; whether to perform quality selecting based on median value of amp of corr-data.
- `median_high`: Float64; max median value (default value).
- `median_low`: Float64; min median value (default value).

