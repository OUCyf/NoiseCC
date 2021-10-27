# NoiseCC.jl: A framwork for large ambient noise cross-correlation (CC) based on SeisNoise.jl

NoiseCC.jl use Multi-process and Multi-thread to accelerate noise cross-correlation (CC) for large dataset. We designed a generic CC framework and wrap SeisNoise module into a Julia interface for the convenience of users, which is suitable for large-scale Dense-Stations, Large-N and DAS data sets. It is convenient for the subsequent study of noise-imaging and noise-monitor.

!!! note

    All the core code for cross-correlation comes from [SeisNoise.jl](https://github.com/tclements/SeisNoise.jl) and [SeisIO Documentation](https://seisio.readthedocs.io/en/latest/), and this package is an example of SeisNoise.jl application. Please read through the [SeisNoise.jl](https://github.com/tclements/SeisNoise.jl) and [SeisIO Documentation](https://seisio.readthedocs.io/en/latest/) to get familiar with seismic data processing and parameters setted for cross-correlation in Julia.


## Installation
Install [SeisNoise](https://github.com/tclements/SeisNoise.jl) firstly. Use the Julia package manager (Press `]` to enter `pkg`):

```julia
julia>]
(@v1.6) pkg> add SeisNoise
```

Then install NoiseCC:

```julia
(@v1.6) pkg> add NoiseCC
```

## Package Features
- Wrap SeisNoise module into a Julia interface for easy using.
- Multi-process and Multi-thread to accelerate noise cross-correlation.
- A standardized seismic data processing framework.

## Author Saying
**Why I wrap SeisNoise into a API in julia?**
- Seismic cross-correlation calculation is the most common task in seismology. It often requires designing different program solutions based on data in different scenarios, SeisNoise was originally designed to facilitate users to write their own programs for different problems. But writing parallel programs for large data is often repetitive and unfriendly to beginners, which means learning a new language and writing your own programs.
- My idea was to provide a common programming interface that would work with most scenarios, which would save the time to rewrite code. So, NoiseCC can be seen as an application case of SeisNoise, which meets my requirement of cross-correlation calculation in most scenarios. I want to make this part of the code available to anyone who needs it.

