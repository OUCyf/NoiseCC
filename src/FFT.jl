# Greetings! from path ~/.julia/config/startup.jl
# Greetings! from path ~/.julia/config/startup.jl
using SeisIO, SeisNoise, Dates, Glob
export FFT




"""
    get_allfile_path(Parameters)

get all files path --> files

# Arguments
- `Parameters::Dict{String,Any}`: ...
"""
function get_allfile_path_fft(
    select_file::Bool,
    FILETXT::String,
    DATADIR::String,
    select_sta::String,
    STATXT::String,
    channel_regular::String,
)
    # get all files path --> files
    if select_file == true
        files = readlines(FILETXT)
    else
        if select_sta == "no"
            files = sort(glob("*/*/*/$(channel_regular)", DATADIR))
        elseif select_sta == "sta"
            # a.glob DATADIR --> rfiles, N_rfiles
            rfiles = sort(glob("*/*/*/$(channel_regular)", DATADIR))
            N_rfiles = size(rfiles, 1)
            # b.read STATXT --> station_txt
            station_txt = readlines(STATXT)
            N_sta = size(station_txt, 1)
            # c.get all stations --> All_NET_STA 
            All_NET_STA = Vector{Vector{String}}(undef, N_sta)
            for ll = 1:N_sta
                D_net, D_sta = split(station_txt[ll])
                All_NET_STA[ll] = [string(D_net), string(D_sta)]
            end
            # d.select data path with selected 'All_NET_STA'  --> files
            files = Vector{String}(undef, 0)
            read_mark = -1
            for jj = 1:N_rfiles
                for ii = 1:N_sta
                    file = rfiles[jj]
                    if occursin(All_NET_STA[ii][1], basename(file)) &&
                       occursin(All_NET_STA[ii][2], basename(file))
                        read_mark = 1
                        push!(files, file)
                    end
                end
            end
            # e.no stations in DATADIR --> files=[]
            if read_mark == -1
                files = []
            end
        end
    end
    return files
end


function get_SeisData_fft(Onefile_path::String, input_fmt::String,freqmin::Float64, freqmax::Float64)
    S = SeisData()
    try
        S_one = read_data(input_fmt, Onefile_path)
        detrend!(S_one)
        taper!(S_one)
        bandpass!(S_one, freqmin, freqmax, corners = 4, zerophase = false)
        append!(S, S_one)
    catch e
        println("Errorpass in function 'get_SeisData_fft':", e)
    end
    return S
end


function get_SeisData2RawData_fft(S::SeisData, cc_len::Int64, cc_step::Int64)
    R = RawData(S, cc_len, cc_step)
    return R
end


function get_time_norm_fft(
    R::RawData,
    time_norm::String,
    factor_clip_std::Float64,
    factor_mute::Float64,
    time_half_win::Int64,
)
    if time_norm == "no"
        nothing
    elseif time_norm == "one_bit"
        onebit!(R)
    elseif time_norm == "clip"
        clip!(R, factor_clip_std)
    elseif time_norm == "mute"
        mute!(R, factor_mute)
    elseif time_norm == "smooth"
        smooth!(R.x, time_half_win)
    end
    return R
end


function get_rfft_fft(R::RawData, cc_method::String)
    Fs = FFTData()
    if cc_method == "CC"
        Fs = rfft(R)
    elseif cc_method == "PCC"
        Fs = phase(R)
    end
    return Fs
end


function get_freq_norm_fft(
    Fs::FFTData,
    freq_norm::String,
    freqmin::Float64,
    freqmax::Float64,
    freq_half_win::Int64,
)
    if freq_norm == "no"
        nothing
    elseif freq_norm == "whiten"
        whiten!(Fs, freqmin, freqmax)
    elseif freq_norm == "coherence"
        coherence!(Fs, freq_half_win)
    elseif freq_norm == "deconvolution"
        deconvolution!(Fs, freq_half_win)
    end
    return Fs
end


function get_save_fft_fft(Fs::FFTData, Onefile_path::String, FFTDIR::String)
    raw_filepath = splitpath(Onefile_path)
    filepath =
        joinpath(FFTDIR, raw_filepath[end-3], raw_filepath[end-2], raw_filepath[end-1])
    if !isdir(filepath)
        mkpath(filepath)
    end
    save_fft(Fs, filepath)
    return nothing
end


function main_fft(Onefile_path::String, Parameters::Dict{String,Any})
    flag = Parameters["flag"]
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
    FFTDIR = Parameters["FFTDIR"]

    println("FFT file Onefile: $(Onefile_path) || now:", now())
    if flag == true
        ###### 1.read sac/mseed -> SeisData() ######
        tALL = @elapsed begin
            S = get_SeisData_fft(Onefile_path, input_fmt, freqmin, freqmax)
        end
        println("'Func:get_SeisData' took $tALL seconds")
        ###### 2.remove response ######
        #...
        ###### 3.SeisData() -> RawData() ######
        tALL = @elapsed begin
            R = get_SeisData2RawData_fft(S, cc_len, cc_step)
        end
        println("'Func:get_SeisData2RawData' took $tALL seconds")
        ###### 4.time_norm ######
        tALL = @elapsed begin
            R = get_time_norm_fft(R, time_norm, factor_clip_std, factor_mute, time_half_win)
        end
        println("'Func:get_time_norm' took $tALL seconds")
        ###### 5.fft ######
        tALL = @elapsed begin
            Fs = get_rfft_fft(R, cc_method)
        end
        println("'Func:get_rfft' took $tALL seconds")
        ###### 6.freq_norm ######
        tALL = @elapsed begin
            Fs = get_freq_norm_fft(Fs, freq_norm, freqmin, freqmax, freq_half_win)
        end
        println("'Func:get_freq_norm' took $tALL seconds")
        ###### 7.save_fft ######
        tALL = @elapsed begin
            get_save_fft_fft(Fs, Onefile_path, FFTDIR)
        end
        println("'Func:get_save_fft' took $tALL seconds")
    else
        ###### 1.read sac/mseed -> SeisData() ######
        S = get_SeisData_fft(Onefile_path, input_fmt, freqmin, freqmax)
        ###### 2.remove response ######
        #...
        ###### 3.SeisData() -> RawData() ######
        R = get_SeisData2RawData_fft(S, cc_len, cc_step)
        ###### 4.time_norm ######
        R = get_time_norm_fft(R, time_norm, factor_clip_std, factor_mute, time_half_win)
        ###### 5.fft ######
        Fs = get_rfft_fft(R, cc_method)
        ###### 6.freq_norm ######
        Fs = get_freq_norm_fft(Fs, freq_norm, freqmin, freqmax, freq_half_win)
        ###### 7.save_fft ######
        get_save_fft_fft(Fs, Onefile_path, FFTDIR)
    end
    return nothing
end

function FFT(
    cores::Int64,
    flag::Bool,
    precompile::Bool,
    FFTDIR::String,
    select_file::Bool,
    FILETXT::String,
    DATADIR::String,
    select_sta::String,
    STATXT::String,
    input_fmt::String,
    time_norm::String,
    freq_norm::String,
    cc_method::String,
    channel_regular::String,
    rm_resp::String,
    respdir::String,
    freqmin::Float64,
    freqmax::Float64,
    cc_len::Int64,
    cc_step::Int64,
    factor_clip_std::Float64,
    factor_mute::Float64,
    time_half_win::Int64,
    freq_half_win::Int64,
)

    ### a. parameters --> Dict-Parameters
    Parameters = Dict(
        "cores" => cores,
        "flag" => flag,
        "FFTDIR" => FFTDIR,
        "select_file" => select_file,
        "FILETXT" => FILETXT,
        "DATADIR" => DATADIR,
        "select_sta" => select_sta,
        "STATXT" => STATXT,
        "input_fmt" => input_fmt,
        "time_norm" => time_norm,
        "freq_norm" => freq_norm,
        "cc_method" => cc_method,
        "channel_regular" => channel_regular,
        "rm_resp" => rm_resp,
        "respdir" => respdir,
        "freqmin" => freqmin,
        "freqmax" => freqmax,
        "cc_len" => cc_len,
        "cc_step" => cc_step,
        "factor_clip_std" => factor_clip_std,
        "factor_mute" => factor_mute,
        "time_half_win" => time_half_win,
        "freq_half_win" => freq_half_win,
    )

    ### b. get_allfile_path
    files =
        get_allfile_path_fft(select_file, FILETXT, DATADIR, select_sta, STATXT, channel_regular)
    N_files = size(files, 1)

    ### c. precompile
    if precompile == true
        println("precompile begin...")
        tALL_pre = @elapsed begin
            pmap(main_fft, [files[1]], [Parameters])
            toremove = glob("*/*/*/*", FFTDIR)
            rm.(toremove)
        end
        println("      precompile took $tALL_pre seconds")
        println("precompile end...")
    end

    ### d. run
    println("\n\n******************* START ******************")
    println("Beginning FFT with $(cores) processors... || N_files: $(N_files)")
    tALL = @elapsed begin
        pmap(main_fft, files, fill(Parameters, N_files))
    end
    println("FFT took $tALL seconds")
    println("******************* END ******************\n\n")
end
