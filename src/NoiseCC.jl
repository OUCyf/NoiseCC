

module NoiseCC

# mode-1: fft
include("FFT.jl")

# mode-1: cc
include("CC.jl")

# mode-2: fft and cc
include("FFT_CC.jl")

# stack
include("STACK.jl")

end
