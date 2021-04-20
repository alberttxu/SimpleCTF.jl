# SimpleCTF.jl

SimpleCTF.jl is a power spectrum estimation program for transmission electron microscope images (CTF estimation).

## Prerequisites
- Install [IMOD](https://bio3d.colorado.edu/imod/)
- Install [Julia](https://julialang.org/)

## Install

`julia> ] add https://github.com/alberttxu/SimpleCTF.jl`

Create a command line alias

`$ alias simplectf="julia <(echo 'import SimpleCTF; SimpleCTF.main()')"`


## Usage
```
$ simplectf --pixelsize 0.87 --voltage 200 --Cs 2.7 input.mrc
input.mrc
clip: Taking power spectrum of 1 slices...
  0.471685 seconds (36 allocations: 1.406 KiB)
preprocessing
  4.885600 seconds (12.44 M allocations: 699.103 MiB, 5.73% gc time, 99.04% compilation time)
searching
  1.367068 seconds (3.55 M allocations: 236.946 MiB, 6.93% gc time, 93.05% compilation time)
z1: 1.0678 μm
z2: 1.0678 μm
θ: -0.46 degrees
correlation: 0.764
 10.449962 seconds (24.43 M allocations: 1.400 GiB, 5.20% gc time, 58.48% compilation time)
```

For the first input file Julia has to JIT compile the code. Successive files are faster if more are given.

```
$ simplectf --pixelsize 0.87 --voltage 200 --Cs 2.7 input.mrc input.mrc
input.mrc
clip: Taking power spectrum of 1 slices...
  0.468964 seconds (36 allocations: 1.406 KiB)
preprocessing
  4.976644 seconds (12.44 M allocations: 699.103 MiB, 5.44% gc time, 99.53% compilation time)
searching
  1.384385 seconds (3.55 M allocations: 236.946 MiB, 7.24% gc time, 92.78% compilation time)
z1: 1.0678 μm
z2: 1.0678 μm
θ: -0.46 degrees
correlation: 0.764
 10.605652 seconds (24.43 M allocations: 1.400 GiB, 5.10% gc time, 58.81% compilation time)

input.mrc
clip: Taking power spectrum of 1 slices...
  0.494393 seconds (36 allocations: 1.406 KiB)
preprocessing
  0.025121 seconds (262.60 k allocations: 11.093 MiB)
searching
  0.118567 seconds (1.54 k allocations: 43.981 MiB, 36.68% gc time)
z1: 1.0678 μm
z2: 1.0678 μm
θ: -0.46 degrees
correlation: 0.764
  0.645198 seconds (264.80 k allocations: 59.013 MiB, 6.74% gc time)
```


## Options
```
$ simplectf -h
usage: 11 --pixelsize PIXELSIZE --voltage VOLTAGE --Cs CS
          [--min_resolution MIN_RESOLUTION]
          [--max_resolution MAX_RESOLUTION]
          [--amplitude_contrast AMPLITUDE_CONTRAST] [--search_phase]
          [-h] images...

positional arguments:
  images                Aligned images for ctf estimation. Can be any
                        format recognized by IMOD.

optional arguments:
  --pixelsize PIXELSIZE
                        Pixel size in Angstroms (type: Float64)
  --voltage VOLTAGE     Accelerating voltage in kV (type: Float64)
  --Cs CS               Spherical abberation in millimeters (type:
                        Float64)
  --min_resolution MIN_RESOLUTION
                        Lower resolution for search (type: Float64,
                        default: 25.0)
  --max_resolution MAX_RESOLUTION
                        Upper resolution for search (type: Float64,
                        default: 5.5)
  --amplitude_contrast AMPLITUDE_CONTRAST
                        Fraction of amplitude contrast (type: Float64,
                        default: 0.07)
  --search_phase        Additionally search for phase shift
  -h, --help            show this help message and exit
```