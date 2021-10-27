# Greetings! from path ~/.julia/config/startup.jl
using SeisIO, SeisNoise, Dates, Glob, Base.Threads, Distributed
export FFT_CC


function get_alltime_path_fftcc(FFTDIR::String, select_time::Bool, TIMETXT::String)
    if select_time == true
        alltime_path = readlines(TIMETXT)
    else
        alltime_path = sort(glob("*/*/*", FFTDIR))
    end
    return alltime_path
end


function get_SR_pairs_fftcc(select_sta::String, STATXT::String, SRTXT::String)
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


function get_allfile_list_fftcc(Onetime_path::String, channel_regular::String)
    files = Vector{String}(undef, 0)
    files_cha = glob(channel_regular, Onetime_path)
    append!(files, files_cha)
    sort!(files)
    return files
end


function get_S_dict_fftcc(S::SeisData)
    S_dict = Dict{String,SeisChannel}()
    for i = 1:S.n
        net, sta, kel, cha = split(S[i].id, '.')
        key = string(net, '.', sta)
        setindex!(S_dict, S[i], key)
    end
    return S_dict
end


function get_SeisData_fftcc(
    Onetime_path::String,
    channel_regular::String,
    select_sta::String,
    STATXT::String,
    SRTXT::String,
    input_fmt::String,
)
    # 1.read one time-chunck paths
    files = get_allfile_list(Onetime_path, channel_regular)
    # 2.read data
    N_files = size(files, 1)
    All_NET_STA = Vector{Vector{String}}(undef, 0)
    SR_pairs = Vector{Vector{String}}(undef, 0)
    S = SeisData()
    if select_sta == "no"
        # a.read all data --> S::SeisData
        for j = 1:N_files
            file = files[j]
            S_one = ""
            try
                s_one = read_data(input_fmt, file)
                S_one = deepcopy(s_one)
            catch e
                println("Errorpass in function 'get_SeisData':", e)
                continue
            end
            append!(S, S_one)
        end
        if S.n == 0   # one time-chunck include 0 files
            return S, 0, 0
        end
        # b.make all stations --> All_NET_STA (for read data[sac/mseed] from database and auto-corr)
        N_sta = S.n
        rAll_NET_STA = Vector{Vector{String}}(undef, 0)
        for j = 1:N_sta
            net, sta, kno, cha = split(S[j].id, '.')
            push!(rAll_NET_STA, [string(net), string(sta)])
        end
        All_NET_STA = deepcopy(rAll_NET_STA)
        # c.make s-r pairs --> SR_pairs (for cross-corr)
        N_SR = Int(N_sta * (N_sta - 1) / 2)
        rSR_pairs = Vector{Vector{String}}(undef, N_SR)
        num = 0
        for ss = 1:N_sta-1
            for rr = ss+1:N_sta
                num += 1
                S_net, S_sta = rAll_NET_STA[ss]
                R_net, R_sta = rAll_NET_STA[rr]
                rSR_pairs[num] = [S_net, S_sta, R_net, R_sta]
            end
        end
        SR_pairs = deepcopy(rSR_pairs)
    else
        # a.read sta-pairs --> rAll_NET_STA, rSR_pairs
        rAll_NET_STA, rSR_pairs = get_SR_pairs(select_sta, STATXT, SRTXT)
        # b.read data from path with selected 'All_NET_STA' --> error_net_sta, S
        error_net_sta = Vector{Vector{String}}(undef, 0)
        read_mark = 1
        for ii = 1:size(rAll_NET_STA, 1)
            for j = 1:N_files
                file = files[j]
                if occursin(rAll_NET_STA[ii][1], basename(file)) &&
                   occursin(rAll_NET_STA[ii][2], basename(file))
                    read_mark = 1
                    try
                        S_one = read_data(input_fmt, file)
                        append!(S, S_one)
                        break
                    catch e
                        push!(error_net_sta, [rAll_NET_STA[ii][1], rAll_NET_STA[ii][2]])   # read error
                        println("Errorpass in function 'get_SeisData':", e)
                        continue
                    end
                else
                    read_mark = -1
                end
            end
            if read_mark == -1   # could not find file with net-sta provided by sta.txt
                push!(error_net_sta, [rAll_NET_STA[ii][1], rAll_NET_STA[ii][2]])
            end
        end
        if S.n == 0   # one time-chunck include 0 files
            return S, 0, 0
        end
        # c.modify rAll_NET_STA, rSR_pairs --> All_NET_STA, SR_pairs
        error_num = length(error_net_sta)
        if error_num == 0
            All_NET_STA = rAll_NET_STA
            SR_pairs = rSR_pairs
        else
            # modify All_NET_STA
            for i = 1:size(rAll_NET_STA, 1)
                net_1, sta_1 = rAll_NET_STA[i][1], rAll_NET_STA[i][2]
                mark_sta = -1
                for j = 1:error_num
                    net_2, sta_2 = error_net_sta[j][1], error_net_sta[j][2]
                    if net_1 == net_2 && sta_1 == sta_2
                        mark_sta = 1 
                        break
                    end
                end
                if mark_sta == -1
                    push!(All_NET_STA, [net_1, sta_1])
                end
            end
            # modify SR_pairs
            for i = 1:size(rSR_pairs, 1)
                s_net_1, s_sta_1, r_net_1, r_sta_1 =
                    rSR_pairs[i][1], rSR_pairs[i][2], rSR_pairs[i][3], rSR_pairs[i][4]
                mark_sta = -1
                for j = 1:error_num
                    net_2, sta_2 = error_net_sta[j][1], error_net_sta[j][2]
                    if (s_net_1 == net_2 && s_sta_1 == sta_2) ||
                        (r_net_1 == net_2 && r_sta_1 == sta_2)
                        mark_sta = 1 
                        break
                    end
                end
                if mark_sta == -1
                    push!(SR_pairs, [s_net_1, s_sta_1, r_net_1, r_sta_1])
                end
            end
        end
    end
    # 3.convert S --> S_dict
    S_dict = get_S_dict(S)

    return S_dict, All_NET_STA, SR_pairs
end


function get_SeisData2RawData_fftcc(S_dict::Dict{String,SeisChannel}, cc_len, cc_step)
    R_dict = Dict{String,RawData}()
    for (k, v) in S_dict
        R = RawData(v, cc_len, cc_step)
        demean!(R)
        taper!(R)
        setindex!(R_dict, R, k)
    end
    return R_dict
end


function get_time_norm_fftcc(
    R_dict::Dict{String,RawData},
    time_norm::String,
    factor_clip_std::Float64,
    factor_mute::Float64,
    time_half_win::Int64,
)
    if time_norm == "no"
        nothing
    elseif time_norm == "one_bit"
        for (k, _) in R_dict
            onebit!(R_dict[k])
        end
    elseif time_norm == "clip"
        for (k, _) in R_dict
            clip!(R_dict[k], factor_clip_std)
        end
    elseif time_norm == "mute"
        for (k, _) in R_dict
            mute!(R_dict[k], factor_mute)
        end
    elseif time_norm == "smooth"
        for (k, _) in R_dict
            smooth!(R_dict[k].x, time_half_win)
        end
    end
    return R_dict
end


function get_rfft_fftcc(R_dict::Dict{String,RawData}, cc_method::String)
    # N_station = R_dict.count
    Fs_dict = Dict{String,FFTData}()
    for (k, _) in R_dict
        if cc_method == "CC"
            Fs = rfft(R_dict[k])
        else
            Fs = phase(R_dict[k])
        end
        setindex!(Fs_dict, Fs, k)
    end
    return Fs_dict
end


function get_freq_norm_fftcc(
    Fs_dict::Dict{String,FFTData},
    freq_norm::String,
    freqmin::Float64,
    freqmax::Float64,
    freq_half_win::Int64,
)
    if freq_norm == "no"
        nothing
    elseif freq_norm == "whiten"
        for (k, _) in Fs_dict
            whiten!(Fs_dict[k], freqmin, freqmax)
        end
    elseif freq_norm == "coherence"
        for (k, _) in Fs_dict
            coherence!(Fs_dict[k], freq_half_win)
        end
    elseif freq_norm == "deconvolution"
        for (k, _) in Fs_dict
            deconvolution!(Fs_dict[k], freq_half_win)
        end
    end
    return Fs_dict
end


function get_cc_substack_fftcc(
    Onetime_path::String,
    Fs_dict::Dict{String,FFTData},
    All_NET_STA::Vector{Vector{String}},
    SR_pairs::Vector{Vector{String}},
    auto_corr::Bool,
    cross_corr::Bool,
    cc_method::String,
    CCFDIR::String,
    maxlag::Float64,
    substack::Bool,
    substack_method::Vector{String},
)
    year = splitpath(Onetime_path)[end-2]
    yearday = splitpath(Onetime_path)[end-1]
    daytime = splitpath(Onetime_path)[end]

    if auto_corr == true
        Threads.@threads for i = 1:size(All_NET_STA, 1)
            # debug:
            println("thread num = ", Threads.threadid())
            #
            net, sta = All_NET_STA[i][1], All_NET_STA[i][2]
            key = string(net, '.', sta)
            C = correlate(Fs_dict[key], Fs_dict[key], maxlag, corr_type = cc_method)
            if substack == true
                for jj in substack_method
                    if jj == "mean"
                        C_mean = stack(C, allstack = true)
                        SPLITOUT_mean =
                            joinpath(CCFDIR, "substack_mean", year, yearday, daytime, key)
                        if !isdir(SPLITOUT_mean)
                            mkpath(SPLITOUT_mean)
                        end
                        save_corr(C_mean, SPLITOUT_mean)
                    end
                    if jj == "pws"
                        C_pws = pws(C)
                        SPLITOUT_pws =
                            joinpath(CCFDIR, "substack_pws", year, yearday, daytime, key)
                        if !isdir(SPLITOUT_pws)
                            mkpath(SPLITOUT_pws)
                        end
                        save_corr(C_pws, SPLITOUT_pws)
                    end
                    if jj == "robust"
                        C_robust = deepcopy(C)
                        C_robust.corr = robustpws(C.corr)
                        SPLITOUT_robust =
                            joinpath(CCFDIR, "substack_robust", year, yearday, daytime, key)
                        if !isdir(SPLITOUT_robust)
                            mkpath(SPLITOUT_robust)
                        end
                        save_corr(C_robust, SPLITOUT_robust)
                    end
                end
            else
                SPLITOUT_no = joinpath(CCFDIR, "substack_no", year, yearday, daytime, key)
                if !isdir(SPLITOUT_no)
                    mkpath(SPLITOUT_no)
                end
                save_corr(C, SPLITOUT_no)
            end
        end
    end

    if cross_corr == true
        Threads.@threads for i = 1:size(SR_pairs, 1)
            s_net, s_sta, r_net, r_sta =
                SR_pairs[i][1], SR_pairs[i][2], SR_pairs[i][3], SR_pairs[i][4]
            key_1 = string(s_net, '.', s_sta)
            key_2 = string(r_net, '.', r_sta)
            C = correlate(Fs_dict[key_1], Fs_dict[key_2], maxlag, corr_type = cc_method)
            if substack == true
                for jj in substack_method
                    if jj == "mean"
                        C_mean = stack(C, allstack = true)
                        SPLITOUT_mean =
                            joinpath(CCFDIR, "substack_mean", year, yearday, daytime, key_1)
                        if !isdir(SPLITOUT_mean)
                            mkpath(SPLITOUT_mean)
                        end
                        save_corr(C_mean, SPLITOUT_mean)
                    end
                    if jj == "pws"
                        C_pws = pws(C)
                        SPLITOUT_pws =
                            joinpath(CCFDIR, "substack_pws", year, yearday, daytime, key_1)
                        if !isdir(SPLITOUT_pws)
                            mkpath(SPLITOUT_pws)
                        end
                        save_corr(C_pws, SPLITOUT_pws)
                    end
                    if jj == "robust"
                        C_robust = deepcopy(C)
                        C_robust.corr = robustpws(C.corr)
                        SPLITOUT_robust = joinpath(
                            CCFDIR,
                            "substack_robust",
                            year,
                            yearday,
                            daytime,
                            key_1,
                        )
                        if !isdir(SPLITOUT_robust)
                            mkpath(SPLITOUT_robust)
                        end
                        save_corr(C_robust, SPLITOUT_robust)
                    end
                end
            else
                SPLITOUT_no = joinpath(CCFDIR, "substack_no", year, yearday, daytime, key_1)
                if !isdir(SPLITOUT_no)
                    mkpath(SPLITOUT_no)
                end
                save_corr(C, SPLITOUT_no)
            end
        end
    end

    return nothing
end


function main_fft_cc(Onetime_path::String, Parameters::Dict{String,Any})
    flag = Parameters["flag"]
    channel_regular = Parameters["channel_regular"]
    select_sta = Parameters["select_sta"]
    STATXT = Parameters["STATXT"]
    SRTXT = Parameters["SRTXT"]
    input_fmt = Parameters["input_fmt"]
    rm_resp = Parameters["rm_resp"]
    respdir = Parameters["respdir"]
    cc_len = Parameters["cc_len"]
    cc_step = Parameters["cc_step"]
    time_norm = Parameters["time_norm"]
    factor_clip_std = Parameters["factor_clip_std"]
    factor_mute = Parameters["factor_mute"]
    time_half_win = Parameters["time_half_win"]
    cc_method = Parameters["cc_method"]
    freq_norm = Parameters["freq_norm"]
    freqmin = Parameters["freqmin"]
    freqmax = Parameters["freqmax"]
    freq_half_win = Parameters["freq_half_win"]
    auto_corr = Parameters["auto_corr"]
    cross_corr = Parameters["cross_corr"]
    CCFDIR = Parameters["CCFDIR"]
    maxlag = Parameters["maxlag"]
    substack = Parameters["substack"]
    substack_method = Parameters["substack_method"]

    println("Cross-correlating Onetime_day: $(Onetime_path) || now:", now())

    if flag == true
        ###### 1.read sac/mseed -> SeisData() ######
        tALL = @elapsed begin
            S_dict, All_NET_STA, SR_pairs = get_SeisData_fftcc(
                Onetime_path,
                channel_regular,
                select_sta,
                STATXT,
                SRTXT,
                input_fmt,
            )
        end
        println("'Func:get_SeisData' took $tALL seconds")

        ###### 2.remove response ######
        #...

        ###### 3.SeisData() -> RawData() ######
        tALL = @elapsed begin
            R_dict = get_SeisData2RawData_fftcc(S_dict, cc_len, cc_step)
        end
        println("'Func:get_SeisData2RawData' took $tALL seconds")

        ###### 4.time_norm ######
        tALL = @elapsed begin
            R_dict = get_time_norm_fftcc(
                R_dict,
                time_norm,
                factor_clip_std,
                factor_mute,
                time_half_win,
            )
        end
        println("'Func:get_time_norm' took $tALL seconds")

        ###### 5.fft ######
        tALL = @elapsed begin
            Fs_dict = get_rfft_fftcc(R_dict, cc_method)
        end
        println("'Func:get_rfft' took $tALL seconds")

        ###### 6.freq_norm ######
        tALL = @elapsed begin
            Fs_dict = get_freq_norm_fftcc(Fs_dict, freq_norm, freqmin, freqmax, freq_half_win)
        end
        println("'Func:get_freq_norm' took $tALL seconds")

        ###### 7.cc_substack ######
        tALL = @elapsed begin
            get_cc_substack_fftcc(
                Onetime_path,
                Fs_dict,
                All_NET_STA,
                SR_pairs,
                auto_corr,
                cross_corr,
                cc_method,
                CCFDIR,
                maxlag,
                substack,
                substack_method,
            )
        end
        println("'Func:get_cc_substack' took $tALL seconds")
    else
        ###### 1.read sac/mseed -> SeisData() ######
        S_dict, All_NET_STA, SR_pairs = get_SeisData_fftcc(
            Onetime_path,
            channel_regular,
            select_sta,
            STATXT,
            SRTXT,
            input_fmt,
        )

        ###### 2.remove response ######
        #...

        ###### 3.SeisData() -> RawData() ######
        R_dict = get_SeisData2RawData_fftcc(S_dict, cc_len, cc_step)

        ###### 4.time_norm ######
        R_dict =
            get_time_norm_fftcc(R_dict, time_norm, factor_clip_std, factor_mute, time_half_win)

        ###### 5.fft ######
        Fs_dict = get_rfft_fftcc(R_dict, cc_method)

        ###### 6.freq_norm ######
        Fs_dict = get_freq_norm_fftcc(Fs_dict, freq_norm, freqmin, freqmax, freq_half_win)

        ###### 7.cc_substack ######
        get_cc_substack_fftcc(
            Onetime_path,
            Fs_dict,
            All_NET_STA,
            SR_pairs,
            auto_corr,
            cross_corr,
            cc_method,
            CCFDIR,
            maxlag,
            substack,
            substack_method,
        )

    end
    return nothing
end


function FFT_CC(
    cores::Int64,
    threads::Int64,
    flag::Bool,
    precompile::Bool,
    DATADIR::String,
    CCFDIR::String,
    select_time::Bool,
    TIMETXT::String,
    select_sta::String,
    STATXT::String,
    SRTXT::String,
    channel_regular::String,
    rm_resp::String,
    respdir::String,
    input_fmt::String,
    time_norm::String,
    freq_norm::String,
    cc_method::String,
    auto_corr::Bool,
    cross_corr::Bool,
    freqmin::Float64,
    freqmax::Float64,
    cc_len::Int64,
    cc_step::Int64,
    maxlag::Float64,
    substack::Bool,
    substack_method::Vector{String},
    factor_clip_std::Float64,
    factor_mute::Float64,
    time_half_win::Int64,
    freq_half_win::Int64,
)

    ### a. parameters --> Dict-Parameters
    Parameters = Dict(
        "cores" => cores,
        "threads" => threads,
        "flag" => flag,
        "DATADIR" => DATADIR,
        "CCFDIR" => CCFDIR,
        "select_time" => select_time,
        "TIMETXT" => TIMETXT,
        "select_sta" => select_sta,
        "STATXT" => STATXT,
        "SRTXT" => SRTXT,
        "rm_resp" => rm_resp,
        "respdir" => respdir,
        "input_fmt" => input_fmt,
        "time_norm" => time_norm,
        "freq_norm" => freq_norm,
        "cc_method" => cc_method,
        "auto_corr" => auto_corr,
        "cross_corr" => cross_corr,
        "channel_regular" => channel_regular,
        "cc_len" => cc_len,
        "cc_step" => cc_step,
        "maxlag" => maxlag,
        "freqmin" => freqmin,
        "freqmax" => freqmax,
        "substack" => substack,
        "substack_method" => substack_method,
        "factor_clip_std" => factor_clip_std,
        "factor_mute" => factor_mute,
        "time_half_win" => time_half_win,
        "freq_half_win" => freq_half_win,
    )

    ### b. get_alltime_path
    alltime_path = get_alltime_path_fftcc(DATADIR, select_time, TIMETXT)
    N_time_chunks = size(alltime_path, 1)

    ### c. precompile
    if precompile == true
        println("precompile begin...")
        tALL_pre = @elapsed begin
            pmap(main_fft_cc, [alltime_path[1]], [Parameters])
            toremove = glob("*/*/*/*/*/*", Parameters["CCFDIR"])
            rm.(toremove)
        end
        println("      precompile took $tALL_pre seconds")
        println("precompile end...")
    end

    ### d. run
    println("\n\n******************* START ******************")
    println(
        "Beginning XCORR with $(cores) processors and each core with $(threads) threads... || N_time_chunks: $(N_time_chunks)",
    )
    tALL = @elapsed begin
        pmap(main_fft_cc, alltime_path, fill(Parameters, N_time_chunks))
    end
    println("XCORR took $tALL seconds")
    println("******************* END ******************\n\n")
end