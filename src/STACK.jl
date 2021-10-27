
using SeisIO, SeisNoise, Dates, Glob, Statistics, Parsers, JLD2, Base.Threads
export STACK


function timestr2unix(timestr::String)
    year = timestr[1:4]
    year_1 = parse(Int64, year)
    yearday = timestr[5:7]
    yearday_1 = parse(Int64, yearday)
    hour = timestr[9:10]
    hour_1 = parse(Int64, hour)
    minute = timestr[12:13]
    minute_1 = parse(Int64, minute)
    second = timestr[15:16]
    second_1 = parse(Int64, second)
    t = datetime2unix(
        DateTime(year_1, 1, 1) +
        Second(yearday_1 * 86400 + hour_1 * 3600 + minute_1 * 60 + second_1),
    )
    return t
end


function datetime2timestr(time_1::DateTime)
    year_1 = Dates.value(Year(time_1))
    month_1 = Dates.value(Month(time_1))
    day_1 = Dates.value(Day(time_1))
    hour_1 = Dates.value(Hour(time_1))
    minute_1 = Dates.value(Minute(time_1))
    second_1 = Dates.value(Second(time_1))
    yearstr = string(year_1, base = 10, pad = 4)
    monthstr = string(month_1, base = 10, pad = 2)
    daystr = string(day_1, base = 10, pad = 2)
    hourstr = string(hour_1, base = 10, pad = 2)
    minutestr = string(minute_1, base = 10, pad = 2)
    secondstr = string(second_1, base = 10, pad = 2)
    timestr = string(
        yearstr,
        '_',
        monthstr,
        '_',
        daystr,
        '_',
        hourstr,
        '_',
        minutestr,
        '_',
        secondstr,
    )
    return timestr
end


function get_SR_pairs_stack(select_sta::String, STATXT::String, SRTXT::String)
    All_NET_STA = Vector{Vector{String}}(undef, 0)
    SR_pairs = Vector{Vector{String}}(undef, 0)
    if select_sta == "no"
        nothing
    elseif select_sta == "sta"
        # a.read STATXT --> station_txt
        station_txt = readlines(STATXT)  # read sta.txt
        # b.make all stations --> rAll_NET_STA (for reading data[sac/mseed] from database and auto-corr)
        N_sta = size(station_txt, 1)
        rAll_NET_STA = Vector{Vector{String}}(undef, N_sta)
        for ll = 1:N_sta
            D_net, D_sta = split(station_txt[ll])
            rAll_NET_STA[ll] = [string(D_net), string(D_sta)]
        end
        # c.make s-r pairs --> rSR_pairs (for cross-corr)
        N_SR = Int(N_sta * (N_sta - 1) / 2)
        rSR_pairs = Vector{Vector{String}}(undef, N_SR)
        num = 0
        for ss = 1:N_sta-1
            for rr = ss+1:N_sta
                num += 1
                S_net, S_sta = split(station_txt[ss])
                R_net, R_sta = split(station_txt[rr])
                rSR_pairs[num] =
                    [string(S_net), string(S_sta), string(R_net), string(R_sta)]
            end
        end
        # d.return
        All_NET_STA = rAll_NET_STA
        SR_pairs = rSR_pairs
    elseif select_sta == "sr"
        # a.read SRTXT --> station_txt
        sr_txt = readlines(SRTXT)
        # b.generate pSR_pairs with repeated sta-pairs --> pSR_pairs
        source_mark = 1   # "1" mark the source of one read-line
        S_net = ""
        S_sta = ""
        pSR_pairs = Vector{Vector{String}}(undef, 0)
        for i = 1:size(sr_txt, 1)
            if sr_txt[i] != ""
                if source_mark == 1
                    s_net, s_sta = split(sr_txt[i])
                    S_net = s_net
                    S_sta = s_sta
                    source_mark = -1   # "-1" mark next line is not the source
                else
                    R_net, R_sta = split(sr_txt[i])
                    if !((S_net == R_net) && (S_sta == R_sta))
                        push!(
                            pSR_pairs,
                            [string(S_net), string(S_sta), string(R_net), string(R_sta)],
                        )
                    end
                end
            else
                source_mark = 1   # mark next line is the source
            end
        end
        # c.delete the repeated pSR_pairs --> SR_pairs
        rSR_pairs = Vector{Vector{String}}(undef, 0)
        repeat_mark = -1   # "-1" means not repeat
        for i = 1:size(pSR_pairs, 1)-1
            S_net_1, S_sta_1, R_net_1, R_sta_1 =
                pSR_pairs[i][1], pSR_pairs[i][2], pSR_pairs[i][3], pSR_pairs[i][4]
            for j = i+1:size(pSR_pairs, 1)
                S_net_2, S_sta_2, R_net_2, R_sta_2 =
                    pSR_pairs[j][1], pSR_pairs[j][2], pSR_pairs[j][3], pSR_pairs[j][4]
                condition_1 =
                    (S_net_1 == S_net_2) &&
                    (S_sta_1 == S_sta_2) &&
                    (R_net_1 == R_net_2) &&
                    (R_sta_1 == R_sta_2)   # Completely corresponding
                condition_2 =
                    (S_net_1 == R_net_2) &&
                    (S_sta_1 == R_sta_2) &&
                    (R_net_1 == S_net_2) &&
                    (R_sta_1 == S_sta_2)   # Reverse S and R, then corresponde
                if condition_1 || condition_2
                    repeat_mark = 1
                    break
                else
                    repeat_mark = -1
                end
            end
            if repeat_mark != 1
                push!(
                    rSR_pairs,
                    [string(S_net_1), string(S_sta_1), string(R_net_1), string(R_sta_1)],
                )
            end
        end
        # The last group will never repeat
        S_net_1, S_sta_1, R_net_1, R_sta_1 =
            pSR_pairs[end][1], pSR_pairs[end][2], pSR_pairs[end][3], pSR_pairs[end][4]
        push!(
            rSR_pairs,
            [string(S_net_1), string(S_sta_1), string(R_net_1), string(R_sta_1)],
        )

        # d.make all stations --> All_NET_STA (for read data[sac/mseed] from FFT-database and auto-corr)
        double_SR_pairs = Vector{Vector{String}}(undef, 0)
        for i = 1:size(rSR_pairs, 1)
            S_net, S_sta, R_net, R_sta =
                rSR_pairs[i][1], rSR_pairs[i][2], rSR_pairs[i][3], rSR_pairs[i][4]
            push!(double_SR_pairs, [S_net, S_sta])
            push!(double_SR_pairs, [R_net, R_sta])
        end
        rAll_NET_STA = Vector{Vector{String}}(undef, 0)
        NET_1 = ""
        STA_1 = ""
        repeat_mark = -1
        for i = 1:size(double_SR_pairs, 1)-1
            net_1, sta_1 = double_SR_pairs[i][1], double_SR_pairs[i][2]
            NET_1, STA_1 = net_1, sta_1
            for j = i+1:size(double_SR_pairs, 1)
                net_2, sta_2 = double_SR_pairs[j][1], double_SR_pairs[j][2]
                if NET_1 == net_2 && STA_1 == sta_2
                    repeat_mark = 1
                    break
                else
                    repeat_mark = -1
                end
            end
            if repeat_mark != 1
                push!(rAll_NET_STA, [string(NET_1), string(STA_1)])
            end
        end
        NET_1, STA_1 = double_SR_pairs[end][1], double_SR_pairs[end][2]
        push!(rAll_NET_STA, [string(NET_1), string(STA_1)])
        # e.conclusion --> All_NET_STA,SR_pairs
        All_NET_STA = rAll_NET_STA
        SR_pairs = rSR_pairs
    end

    ### modify sequence of SR_pairs
    net_1,sta_1,net_2,sta_2 = "","","",""
    for ii = 1:size(SR_pairs,1)
        net_1,sta_1,net_2,sta_2 = SR_pairs[ii][1],SR_pairs[ii][2],SR_pairs[ii][3],SR_pairs[ii][4]
        if net_1*sta_1 > net_2*sta_2
            SR_pairs[ii][1] = net_2
            SR_pairs[ii][2] = sta_2
            SR_pairs[ii][3] = net_1
            SR_pairs[ii][4] = sta_1
        end
    end
    return sort(All_NET_STA), sort(SR_pairs)
end


# Output: All_NET_STA, SR_pairs --> STA_PAIRs
function get_sta_pairs_stack(select_sta::String, STATXT::String, SRTXT::String,auto_stack::Bool, cross_stack::Bool)
    All_NET_STA, SR_pairs = get_SR_pairs_stack(select_sta, STATXT, SRTXT)
    N = size(All_NET_STA, 1)
    SS_pairs = Vector{Vector{String}}(undef, N)
    for i = 1:N
        SS_pairs[i] =
            [All_NET_STA[i][1], All_NET_STA[i][2], All_NET_STA[i][1], All_NET_STA[i][2]]
    end
    if auto_stack == true && cross_stack == true
        STA_PAIRs = vcat(SS_pairs, SR_pairs)
    elseif auto_stack == true && cross_stack != true
        STA_PAIRs = SS_pairs
    elseif auto_stack != true && cross_stack == true
        STA_PAIRs = SR_pairs
    elseif auto_stack != true && cross_stack != true
        STA_PAIRs = []
    end
    return sort(STA_PAIRs)
end


function get_alltime_path_stack(CCFDIR::String, select_time::Bool, TIMETXT::String)
    if select_time == true
        alltime_path = readlines(TIMETXT)
    else
        alltime_path = sort(glob("*/*/*/*", CCFDIR))
    end
    return alltime_path
end


function get_allfile_list_stack(
    SR_one::Vector{String},
    alltime_path::Vector{String},
    channel_regular::String,
)
    SR_regular = string(".*", SR_one[1], ".*", SR_one[2], ".*", SR_one[3], ".*", SR_one[4])
    # a.find 'channel_regular-corr' path from 'alltime_path'
    files = Vector{String}(undef, 0)
    for i in alltime_path
        files_cha = glob(channel_regular, joinpath(i, SR_one[1] * '.' * SR_one[2]))   # regular channel
        for j in files_cha
            if occursin(Regex("$SR_regular"), basename(j))   # specify source-receiver name
                push!(files, j)
            end
        end
    end
    # b.get time_beg time_end
    N = size(files, 1)
    if N == 0   # if find 0 files
        N_glob = 0
        time_min = unix2datetime(timestr2unix("1970001_00_00_00"))
        time_max = unix2datetime(timestr2unix("1970001_00_00_00"))
        return files, N_glob, time_min, time_max
    else
        N_glob = N
        timestr_min = string(split(splitpath(files[1])[end-2], 'T')[1])
        timestr_max = string(split(splitpath(files[end])[end-2], 'T')[2])
        time_min = unix2datetime(timestr2unix(timestr_min))
        time_max = unix2datetime(timestr2unix(timestr_max))
        return files, N_glob, time_min, time_max
    end
end


# files will be divided into N pieces for N time-chunck stack
function get_Nfile_list_stack(
    SR_one::Vector{String},
    alltime_path::String,
    time_interval::Int64,
    min_stack_chunck::Int64,
    channel_regular::String,
)
    ##################
    # 1. find begin-time and end-time for all files -> files, time_min, time_max
    files, _, time_min, time_max = get_allfile_list_stack(SR_one, alltime_path, channel_regular)
    ##################
    # 2. divide all files into N pieces with 'time_interval' -> time_range, NTime_min, NTime_max, N_chunk
    t_1 = datetime2unix(time_min)
    t_2 = datetime2unix(time_max)
    time_range = collect(range(t_1, step = Float64(time_interval), stop = t_2))
    Ntime_min = unix2datetime.(time_range)
    Ntime_max = unix2datetime.(time_interval .+ time_range)
    N_chunk = size(time_range, 1)
    ##################
    # 3. find N pieces corr-data into -> Nfiles, count_chunk
    Nfiles = Vector{Vector{String}}(undef, N_chunk) # every element store one time-chunck corr-data
    for ii = 1:size(Nfiles, 1)
        Nfiles[ii] = []
    end
    count_chunk = zeros(Int64, N_chunk) # mark count for enough corr-data in one time-chunck
    for i = 1:size(files, 1)  # loop all files
        filetime = split(splitpath(files[i])[end-2], 'T')
        t_min = timestr2unix(string(filetime[1]))
        t_max = timestr2unix(string(filetime[2]))
        kk = -1   # mark the location to store corr-data in Nfiles
        for jj = 1:(N_chunk-1)   # loop day_range, find the location in time_range for every corr-data, and count the number for every time_range
            if (t_min >= time_range[jj]) && (t_max <= time_range[jj+1])
                count_chunk[jj] += 1
                kk = jj
                break
            end
        end
        if kk != -1
            push!(Nfiles[kk], files[i])
        end
    end
    ##################
    # 4. select corr-data for every time_range with 'min_stack_chunck' -> Nfiles, count_chunk, Ntime_min, Ntime_max 
    index = findall(>=(min_stack_chunck), count_chunk)
    if size(index, 1) == 0
        count_chunk = [0]
        Nfiles = [0]
        Ntime_min = [unix2datetime(timestr2unix("1970001_00_00_00"))]
        Ntime_max = [unix2datetime(timestr2unix("1970001_00_00_00"))]
    else
        count_chunk = count_chunk[index]
        Nfiles = Nfiles[index]
        Ntime_min = Ntime_min[index]
        Ntime_max = Ntime_max[index]
    end
    ##################
    return Nfiles, count_chunk, Ntime_min, Ntime_max
end


function get_CorrData_stack(
    files::Vector{String},
    comp::String,
    stack_method::String,
    substack_method::String,
    select_filter::Bool,
    freqmin::Float64,
    freqmax::Float64,
    corners::Int64,
    zerophase::Bool,
)
    ##################
    # 1. read corr -> corr Array
    # try load one file and allocate; if find 0 file, N_read=0 will be returned, then main-function will kill the stacking for this station-pair.
    C_mark = 0
    C_one = CorrData()
    for f in files
        try
            C = load_corr(f, comp)
            C_one = deepcopy(C)
            C_mark = 1
        catch e
            println("Errorpass(init) in function get_CorrData:", e)
        end
        if C_mark == 1  # if read one C, break
            break
        end
    end
    # judge whether read one C-files successfully
    if C_mark == 0
        N_read = 0
        return 0, N_read, 0
    end
    ##################
    # 2.loop readcorr and stack in raw one time-chunck (if substack_method=="no")
    index = Vector{Int64}(undef, 0)
    corr = Array{eltype(C_one.corr)}(undef, size(C_one.corr, 1), length(files))
    for ii = 1:size(files, 1)
        # a.try read and remove nan
        C_one = CorrData()
        try
            C = load_corr(files[ii], comp)
            C_one = deepcopy(C)
            remove_nan!(C_one)
        catch e
            println("Errorpass in function get_CorrData:", e)
            push!(index, ii)   # count NAN in corr
            continue
        end
        # b.filter
        if select_filter == true
            bandpass!(C_one, freqmin, freqmax, corners = corners, zerophase = zerophase)
        end
        # c.substack_method
        if substack_method == "no"
            if stack_method == "mean"
                stack!(C_one, allstack = true)
            elseif stack_method == "pws"
                pws!(C_one)
            elseif stack_method == "robust"
                C_one.corr = robustpws(C_one.corr)
            end
            corr[:, ii] = C_one.corr[:, 1]
        else
            corr[:, ii] = C_one.corr[:, 1]
        end
    end
    corr_no_NAN = corr[:, setdiff(1:end, index)]   # delete column with NAN
    N_read = size(corr_no_NAN, 2)
    ##################
    # 3. get station-pair infomation
    stack_paras = Dict(
        "azi" => C_one.azi,
        "baz" => C_one.baz,
        "dist" => C_one.dist,
        "lon" => C_one.loc.lon,
        "lat" => C_one.loc.lat,
        "maxlag" => C_one.maxlag,
        "cc_len" => C_one.cc_len,
        "cc_step" => C_one.cc_step,
        "corr_type" => C_one.corr_type,
        "comp" => C_one.comp,
        "name" => C_one.name,
    )
    return corr_no_NAN, N_read, stack_paras
end


function mute_outlier_corr_stack(A::AbstractArray, high::Float64, low::Float64)
    @assert low < high "low must be less than high"
    maxamp = vec(maximum(abs.(A), dims = 1))
    # remove Nans
    maxamp[isnan.(maxamp)] .= Inf
    medianmax = median(maxamp)
    ind = findall(x -> low * medianmax <= x <= high * medianmax, maxamp)
    A = A[:, ind]
    N_good = size(A, 2)
    return A, N_good
end


function get_stack_stack(corr::Array{Float32}, stack_method::String)
    if stack_method == "mean"
        cc_satck = mean(corr, dims = 2)
    elseif stack_method == "pws"
        cc_satck = pws(corr)
    elseif stack_method == "robust"
        cc_satck = robustpws(corr)
    end
    return cc_satck
end


function get_savestack_stack(
    SR_one::Vector{String},
    cc_stack::AbstractArray,
    N_info::Vector{Int64},
    stack_paras::Dict{String,Any},
    time_min::DateTime,
    time_max::DateTime,
    STACKDIR::String,
    stack_method::String,
    stack_all::Bool,
    comp::String,
)
    time_min_str = datetime2timestr(time_min)
    time_max_str = datetime2timestr(time_max)
    setindex!(stack_paras, N_info[1], "N_glob")
    setindex!(stack_paras, N_info[2], "N_read")
    setindex!(stack_paras, N_info[3], "N_good")
    # 1.outpath
    S_name = string(SR_one[1], '.', SR_one[2])
    S_R_name = string(SR_one[1], '.', SR_one[2], "_", SR_one[3], '.', SR_one[4])
    if stack_method == "mean"
        outpath = joinpath(STACKDIR, "stack_mean", S_name, S_R_name)
    elseif stack_method == "pws"
        outpath = joinpath(STACKDIR, "stack_pws", S_name, S_R_name)
    else
        stack_method == "robust"
        outpath = joinpath(STACKDIR, "stack_robust", S_name, S_R_name)
    end
    if !isdir(outpath)
        mkpath(outpath)
    end
    # 2.output
    if stack_all == true
        name = string(S_R_name, "__", time_min_str, "T", time_max_str, "__all", ".jld2")
    else
        name = string(S_R_name, "__", time_min_str, "T", time_max_str, ".jld2")
    end
    filename = joinpath(outpath, name)
    file = jldopen(filename, "a+")
    rtime = string(time_min_str, "T", time_max_str)  # start time and end time
    if !(comp in keys(file))
        group = JLD2.Group(file, comp)
        group["data"] = cc_stack
        group["time"] = rtime
        for keyword in keys(stack_paras)
            group[keyword] = stack_paras[keyword]
        end
    else
        file[comp]["data"] = cc_stack
        file[comp]["time"] = rtime
        for keyword in keys(stack_paras)
            file[comp][keyword] = stack_paras[keyword]
        end
    end
    close(file)
end


function main_stack(
    SR_one::Vector{String},
    alltime_path::Vector{String},
    Parameters::Dict{String,Any},
)
    flag = Parameters["flag"]
    stack_all = Parameters["stack_all"]
    select_corr = Parameters["select_corr"]
    time_interval = Parameters["time_interval"]
    min_stack_chunck = Parameters["min_stack_chunck"]
    channel_regular = Parameters["channel_regular"]
    comp = Parameters["comp"]
    stack_method = Parameters["stack_method"]
    substack_method = Parameters["substack_method"]
    select_filter = Parameters["select_filter"]
    freqmin = Parameters["freqmin"]
    freqmax = Parameters["freqmax"]
    corners = Parameters["corners"]
    zerophase = Parameters["zerophase"]
    median_high = Parameters["median_high"]
    median_low = Parameters["median_low"]
    STACKDIR = Parameters["STACKDIR"]

    SR_one_name = string(SR_one[1], '.', SR_one[2], '_', SR_one[3], '.', SR_one[4])
    println("Stacking... $(SR_one_name) || now:", now())

    if flag == true
        if stack_all == true
            ###### 1.read path for one sta-pairs ######
            tALL = @elapsed begin
                files, N_glob, time_min, time_max =
                    get_allfile_list_stack(SR_one, alltime_path, channel_regular)
                if N_glob < min_stack_chunck
                    println(
                        "   ",
                        "N_glob: $(N_glob) < min_stack_chunck: $(min_stack_chunck)",
                        "  -->  $(SR_one_name)__$(time_min)__$(time_max)__all.jld2",
                    )
                    return nothing
                end
            end
            println("'Func:get_allfile_list' took $tALL seconds")

            ###### 2.read corr-data for one sta-pairs ######
            tALL = @elapsed begin
                corr, N_read, stack_paras = get_CorrData_stack(
                    files,
                    comp,
                    stack_method,
                    substack_method,
                    select_filter,
                    freqmin,
                    freqmax,
                    corners,
                    zerophase,
                )
                if N_read < min_stack_chunck
                    println(
                        "   ",
                        "N_read: $(N_read) < min_stack_chunck: $(min_stack_chunck)",
                        "  -->  $(SR_one_name)__$(time_min)__$(time_max)__all.jld2",
                    )
                    return nothing
                end
            end
            println("'Func:get_CorrData' took $tALL seconds")
            ###### 3.select corr-data ######
            tALL = @elapsed begin
                N_good = -1
                if select_corr == true
                    corr, N_good = mute_outlier_corr_stack(corr, median_high, median_low)
                    if N_good < min_stack_chunck
                        println(
                            "   ",
                            "N_good: $(N_good) < min_stack_chunck: $(min_stack_chunck)",
                            "  -->  $(SR_one_name)__$(time_min)__$(time_max)__all.jld2",
                        )
                        return nothing
                    end
                end
            end
            println("'Func:mute_outlier_corr' took $tALL seconds")
            ###### 4.stack ######
            tALL = @elapsed begin
                cc_stack = get_stack_stack(corr, stack_method)
            end
            println("'Func:get_stack' took $tALL seconds")
            ###### 5.save stack-data ######
            tALL = @elapsed begin
                N_info = [N_glob, N_read, N_good]
                get_savestack_stack(
                    SR_one,
                    cc_stack,
                    N_info,
                    stack_paras,
                    time_min,
                    time_max,
                    STACKDIR,
                    stack_method,
                    stack_all,
                    comp,
                )
            end
            println("'Func:get_savestack' took $tALL seconds")
        else
            ###### 1.read Npath for one sta-pairs ######
            tALL = @elapsed begin
                Nfiles, NN_glob, Ntime_min, Ntime_max = get_Nfile_list_stack(
                    SR_one,
                    alltime_path,
                    time_interval,
                    min_stack_chunck,
                    channel_regular,
                )
                println("'Func:get_Nfile_list' took $tALL seconds")
            end
            ###### 2.loop time-chunck ######
            tALL = @elapsed begin
                Threads.@threads for n = 1:size(NN_glob, 1)
                    ###### 2.1.judge enough corr-data in one time-chunck glob ######
                    N_glob = NN_glob[n]
                    if N_glob < min_stack_chunck
                        println(
                            "   ",
                            "N_glob: $(N_glob) < min_stack_chunck: $(min_stack_chunck)",
                            "  -->  $(SR_one_name)__$(Ntime_min[n])__$(Ntime_max[n])",
                        )
                        continue
                    end

                    ###### 2.2.read corr-data for one sta-pairs ######
                    corr, N_read, stack_paras = get_CorrData_stack(
                        Nfiles[n],
                        comp,
                        stack_method,
                        substack_method,
                        select_filter,
                        freqmin,
                        freqmax,
                        corners,
                        zerophase,
                    )
                    if N_read < min_stack_chunck
                        println(
                            "   ",
                            "N_read: $(N_read) < min_stack_chunck: $(min_stack_chunck)",
                            "  -->  $(SR_one_name)__$(Ntime_min[n])__$(Ntime_max[n])",
                        )
                        continue
                    end

                    ###### 2.3.select corr-data ######
                    N_good = -1
                    if select_corr == true
                        corr, N_good = mute_outlier_corr_stack(corr, median_high, median_low)
                        if N_good < min_stack_chunck
                            println(
                                "   ",
                                "N_good: $(N_good) < min_stack_chunck: $(min_stack_chunck)",
                                "  -->  $(SR_one_name)__$(Ntime_min[n])__$(Ntime_max[n])",
                            )
                            continue
                        end
                    end

                    ###### 2.4.stack ######
                    cc_stack = get_stack_stack(corr, stack_method)

                    ###### 2.5.save stack-data ######
                    N_info = [N_glob, N_read, N_good]
                    get_savecorr_stack(
                        SR_one,
                        cc_stack,
                        N_info,
                        stack_paras,
                        Ntime_min[n],
                        Ntime_max[n],
                        STACKDIR,
                        stack_method,
                        stack_all,
                        comp,
                    )
                end
            end
            println("'Func:process...1-time-chunck' took $tALL seconds")
        end
    else
        if stack_all == true
            ###### 1.read path for one sta-pairs ######
            files, N_glob, time_min, time_max =
                get_allfile_list_stack(SR_one, alltime_path, channel_regular)
            if N_glob < min_stack_chunck
                println(
                    "   ",
                    "N_glob: $(N_glob) < min_stack_chunck: $(min_stack_chunck)",
                    "  -->  $(SR_one_name)__$(time_min)__$(time_max)__all.jld2",
                )
                return nothing
            end


            ###### 2.read corr-data for one sta-pairs ######
            corr, N_read, stack_paras = get_CorrData_stack(
                files,
                comp,
                stack_method,
                substack_method,
                select_filter,
                freqmin,
                freqmax,
                corners,
                zerophase,
            )
            if N_read < min_stack_chunck
                println(
                    "   ",
                    "N_read: $(N_read) < min_stack_chunck: $(min_stack_chunck)",
                    "  -->  $(SR_one_name)__$(time_min)__$(time_max)__all.jld2",
                )
                return nothing
            end

            ###### 3.select corr-data ######
            N_good = -1
            if select_corr == true
                corr, N_good = mute_outlier_corr_stack(corr, median_high, median_low)
                if N_good < min_stack_chunck
                    println(
                        "   ",
                        "N_good: $(N_good) < min_stack_chunck: $(min_stack_chunck)",
                        "  -->  $(SR_one_name)__$(time_min)__$(time_max)__all.jld2",
                    )
                    return nothing
                end
            end

            ###### 4.stack ######
            cc_stack = get_stack_stack(corr, stack_method)

            ###### 5.save stack-data ######
            N_info = [N_glob, N_read, N_good]
            get_savestack_stack(
                SR_one,
                cc_stack,
                N_info,
                stack_paras,
                time_min,
                time_max,
                STACKDIR,
                stack_method,
                stack_all,
                comp,
            )

        else
            ###### 1.read Npath for one sta-pairs ######
            Nfiles, NN_glob, Ntime_min, Ntime_max = get_Nfile_list_stack(
                SR_one,
                alltime_path,
                time_interval,
                min_stack_chunck,
                channel_regular,
            )

            ###### 2.loop time-chunck ######
            Threads.@threads for n = 1:size(NN_glob, 1)
                ###### 2.1.judge enough corr-data in one time-chunck glob ######
                N_glob = NN_glob[n]
                if N_glob < min_stack_chunck
                    println(
                        "   ",
                        "N_glob: $(N_glob) < min_stack_chunck: $(min_stack_chunck)",
                        "  -->  $(SR_one_name)__$(Ntime_min[n])__$(Ntime_max[n])",
                    )
                    continue
                end

                ###### 2.2.read corr-data for one sta-pairs ######
                corr, N_read, stack_paras = get_CorrData_stack(
                    Nfiles[n],
                    comp,
                    stack_method,
                    substack_method,
                    select_filter,
                    freqmin,
                    freqmax,
                    corners,
                    zerophase,
                )
                if N_read < min_stack_chunck
                    println(
                        "   ",
                        "N_read: $(N_read) < min_stack_chunck: $(min_stack_chunck)",
                        "  -->  $(SR_one_name)__$(Ntime_min[n])__$(Ntime_max[n])",
                    )
                    continue
                end

                ###### 2.3.select corr-data ######
                N_good = -1
                if select_corr == true
                    corr, N_good = mute_outlier_corr_stack(corr, median_high, median_low)
                    if N_good < min_stack_chunck
                        println(
                            "   ",
                            "N_good: $(N_good) < min_stack_chunck: $(min_stack_chunck)",
                            "  -->  $(SR_one_name)__$(Ntime_min[n])__$(Ntime_max[n])",
                        )
                        continue
                    end
                end

                ###### 2.4.stack ######
                cc_stack = get_stack_stack(corr, stack_method)

                ###### 2.5.save stack-data ######
                N_info = [N_glob, N_read, N_good]
                get_savecorr_stack(
                    SR_one,
                    cc_stack,
                    N_info,
                    stack_paras,
                    Ntime_min[n],
                    Ntime_max[n],
                    STACKDIR,
                    stack_method,
                    stack_all,
                    comp,
                )
            end
        end

    end
end


function STACK(
    cores::Int64,
    threads::Int64,
    flag::Bool,
    precompile::Bool,
    CCFDIR::String,
    STACKDIR::String,
    select_time::Bool,
    TIMETXT::String,
    select_sta::String,
    STATXT::String,
    SRTXT::String,
    channel_regular::String,
    comp::String,
    stack_all::Bool,
    time_interval::Int64,
    min_stack_chunck::Int64,
    substack_method::String,
    stack_method::String,
    auto_stack::Bool,
    cross_stack::Bool,
    select_filter::Bool,
    freqmin::Float64,
    freqmax::Float64,
    corners::Int64,
    zerophase::Bool,
    select_corr::Bool,
    median_high::Float64,
    median_low::Float64,
)
    ### a. parameters --> Dict-Parameters
    Parameters = Dict(
        "cores" => cores,
        "threads" => threads,
        "flag" => flag,
        "CCFDIR" => CCFDIR,
        "STACKDIR" => STACKDIR,
        "select_time" => select_time,
        "TIMETXT" => TIMETXT,
        "select_sta" => select_sta,
        "STATXT" => STATXT,
        "SRTXT" => SRTXT,
        "channel_regular" => channel_regular,
        "comp" => comp,
        "stack_all" => stack_all,
        "time_interval" => time_interval,
        "min_stack_chunck" => min_stack_chunck,
        "substack_method" => substack_method,
        "stack_method" => stack_method,
        "auto_stack" => auto_stack,
        "cross_stack" => cross_stack,
        "select_filter" => select_filter,
        "freqmin" => freqmin,
        "freqmax" => freqmax,
        "corners" => corners,
        "zerophase" => zerophase,
        "select_corr" => select_corr,
        "median_high" => median_high,
        "median_low" => median_low,
    )


    ### b. generate STA_PAIRs,alltime_path
    alltime_path = get_alltime_path_stack(CCFDIR, select_time, TIMETXT)
    STA_PAIRs = get_sta_pairs_stack(select_sta, STATXT, SRTXT, auto_stack, cross_stack)
    N_sta_pairs = size(STA_PAIRs, 1)


    ### c. precompile
    if precompile == true
        println("precompile begin...")
        tALL_pre = @elapsed begin
            pmap(main_stack, [STA_PAIRs[1]], [alltime_path], [Parameters])
            stack_name = string("stack_", Parameters["stack_method"])
            toremove = glob("$(stack_name)/*/*/*", Parameters["STACKDIR"])
            rm.(toremove)
        end
        println("      precompile took $tALL_pre seconds")
        println("precompile end...")
    end

    ### d. run-loop all SR_pairs (pmap)
    println("\n\n******************* START ******************")
    println(
        "Beginning STACK with $(cores) processors and each core with $(threads) threads... || N_SR_pairs: $(N_sta_pairs)",
    )
    tALL = @elapsed begin
        pmap(
            main_stack,
            STA_PAIRs,
            fill(alltime_path, N_sta_pairs),
            fill(Parameters, N_sta_pairs),
        )
    end
    println("STACK took $tALL seconds")
    println("******************* END ******************")
end


