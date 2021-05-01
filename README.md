# SimpleCTF.jl

SimpleCTF.jl is a power spectrum estimation program for transmission electron microscope images (CTF estimation).

## Prerequisites
- Julia 1.6
- Input files must be tif format

## Install

`julia> ] add https://github.com/alberttxu/SimpleCTF.jl`

Create a command line alias

```
$ alias simplectf="julia <(echo 'import SimpleCTF; SimpleCTF.main()')"
$ export JULIA_NUM_THREADS=8
```


## Usage
```
$ simplectf --pixelsize 0.87 --voltage 200 --Cs 2.7 input.tif
processing input.tif
 14.544366 seconds (31.21 M allocations: 1.999 GiB, 7.29% gc time, 0.14% compilation time)
results written to simplectf_results.csv

$ cat simplectf_results.csv
filename,defocus1,defocus2,astigmatism angle,phase shift,correlation
input.tif,1.0671,1.0671,-2.4644,0.0,0.2242
```

For the first input file Julia has to JIT compile the code. Successive files are faster if more are given.

## Options
```
usage: 11 --pixelsize PIXELSIZE --voltage VOLTAGE --Cs CS
          [--min_resolution MIN_RESOLUTION]
          [--max_resolution MAX_RESOLUTION]
          [--amplitude_contrast AMPLITUDE_CONTRAST]
          [--do_phase_search] [--csv CSV] [-h] images...

positional arguments:
  images                Aligned images for ctf estimation. Must be tif
                        format.

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
  --do_phase_search     Additionally search for phase shift
  --csv CSV             csv to store results (default:
                        "simplectf_results.csv")
  -h, --help            show this help message and exit
```