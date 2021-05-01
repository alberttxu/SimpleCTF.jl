module SimpleCTF

export find_ctf
include("utils.jl")

using ArgParse
using LinearAlgebra: BLAS
using Base.Threads

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--pixelsize"
            help = "Pixel size in Angstroms"
            arg_type = Float64
            required = true
        "--voltage"
            help = "Accelerating voltage in kV"
            arg_type = Float64
            required = true
        "--Cs"
            help = "Spherical abberation in millimeters"
            arg_type = Float64
            required = true
        "--min_resolution"
            help = "Lower resolution for search"
            arg_type = Float64
            default = 25.
        "--max_resolution"
            help = "Upper resolution for search"
            arg_type = Float64
            default = 5.5
        "--amplitude_contrast"
            help = "Fraction of amplitude contrast"
            arg_type = Float64
            default = 0.07
        "--do_phase_search"
            help = "Additionally search for phase shift"
            action = :store_true
        "--csv"
            help = "csv to store results"
            default = "simplectf_results.csv"
        "images"
            nargs = '*'
            default = []
            help = "Aligned images for ctf estimation. Must be tif format."
            required = true
    end
    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()

    pixelsize = parsed_args["pixelsize"]
    voltage = parsed_args["voltage"]
    Cs = parsed_args["Cs"]
    min_resolution = parsed_args["min_resolution"]
    max_resolution = parsed_args["max_resolution"]
    amplitude_contrast = parsed_args["amplitude_contrast"]
    do_phase_search = parsed_args["do_phase_search"]
    out_csv = parsed_args["csv"]
    input_images = parsed_args["images"]

    BLAS.set_num_threads(1)
    FFTW.set_num_threads(1)

    results = []
    @time @threads for input_image in input_images
        println("processing $input_image")
        reference_spectrum, estimated_spectrum, search_result = find_ctf(
            input_image,
            pixelsize,
            voltage,
            Cs,
            min_resolution,
            max_resolution,
            amplitude_contrast,
            do_phase_search,
        )
        push!(results, (input_image, search_result))
        #diagnostic = [reference_spectrum[end:-1:1, :]; estimated_spectrum]
        #diagnostic ./= maximum(diagnostic)
    end

    f = open(out_csv, "w")
    write(f, "filename,defocus1,defocus2,astigmatism angle,phase shift,correlation\n")
    for (input_image, search_result) in results
        write(f, "$input_image,$search_result\n")
    end
    close(f)
    println("results written to $out_csv")

end

end
