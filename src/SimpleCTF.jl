module SimpleCTF

using ArgParse
include("utils.jl")


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--pixelsize"
            help = "Pixel size in Angstroms"
            arg_type = Float64
        "--voltage"
            help = "Accelerating voltage in kV"
            arg_type = Float64
        "--Cs"
            help = "Spherical abberation in millimeters"
            arg_type = Float64
        "--min_resolution"
            help = "Lower resolution for search"
            arg_type = Float64
            default = 25.0
        "--max_resolution"
            help = "Upper resolution for search"
            arg_type = Float64
            default = 5.5
        "--amplitude_contrast"
            help = "Fraction of amplitude contrast"
            arg_type = Float64
            default = 0.07
        "--search_phase"
            help = "Search for phase shift"
            action = :store_true
        "images"
            nargs = '*'
            default = []
            help = "Aligned images for ctf estimation. Can be any format recognized by IMOD."
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
    search_phase = parsed_args["search_phase"]
    mrc_images = parsed_args["mrc_images"]

    for mrc_image in mrc_images
        println(mrc_image)
        @time reference_spectrum, estimated_spectrum = find_ctf(
            mrc_image,
            pixelsize,
            voltage,
            Cs,
            min_resolution,
            max_resolution,
            amplitude_contrast,
            search_phase,
        )
        println()
        #diagnostic = [reference_spectrum[end:-1:1, :]; estimated_spectrum]
        #diagnostic ./= maximum(diagnostic)
    end
end


export find_ctf

end
