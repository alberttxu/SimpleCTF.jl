using LinearAlgebra
using FFTW, Images, ImageIO, TiffImages, ImageFiltering, Interpolations, NLopt


struct Range2D{T}
    x::UnitRange{T}
    y::UnitRange{T}
end


"""
Given an accelerating voltage in kV, return the electron wavelength in picometers.
Taken from ctffind4 code.
"""
function wavelength_from_voltage(voltage)
    1226.39 / sqrt(1000.0 * voltage + 0.97845e-6 * (1000.0 * voltage)^2)
end


""" Return a log averaged-spectrum using Bartlett's method. """
function get_spectrum(img, N::Integer = 1024, stride = 1024)
    m, n = size(img)
    A = abs.(rfft(@view img[1:N, 1:N]))
    for j in stride+1:stride:n-N, i in stride+1:stride:m-N
        A .+= abs.(rfft(@view img[i:i+1023, j:j+1023]))
    end
    return log.(A)
end


"""
Given the power spectrum of an image (rfft without fftshift),
return a 1-D profile in the horizontal direction.

width: number of rows to average
σ: std. deviation of gaussian smoothing
"""
function profile1D(spectrum::AbstractArray{<:Real}; σ = 2, width::Integer = 15)
    j_max = div(size(spectrum)[2], 2)
    horizontal_strip = @view spectrum[1:width, 1:j_max]
    profile = vec(sum(horizontal_strip, dims = 1)) / width
    if σ > 0
        imfilter!(profile, profile, KernelFactors.IIRGaussian((σ,)))
    end
    return profile
end


""" Linearly interpolates the min peaks of a 1-D profile.

Returns the interpolation as a function of the array indices.
"""
function get_background(search_profile::AbstractArray)
    # 2 iterations to get the min peaks to zero
    n_iters = 2
    itp = Vector{Any}(undef, n_iters)
    for i = 1:n_iters
        local_mins = [x[1] for x in findlocalminima(search_profile)]
        itp[i] = LinearInterpolation(
            local_mins,
            search_profile[local_mins],
            extrapolation_bc = Line(),
        )
        lower_envelope = itp[i].(1:length(search_profile))
        search_profile -= lower_envelope
    end
    return x -> sum(itp[i](x) for i = 1:n_iters)
end


""" Linearly interpolates the thon ring peaks of a background-noise subtracted 1-D profile.

Returns the interpolation as a function of the array indices.
"""
function get_envelope_function(subtracted_profile::AbstractArray)
    local_maxes = [x[1] for x in findlocalmaxima(subtracted_profile)]
    return LinearInterpolation(
        local_maxes,
        subtracted_profile[local_maxes],
        extrapolation_bc = Line(),
    )
end


""" Given x and y (2D meshgrid arrays), return a bandpass binary mask. """
function create_binary_mask(x::AbstractArray{Int,2}, y::AbstractArray{Int,2}, low, high)
    r = @. sqrt(x^2 + y^2)
    mask = low .< r .< high
    return mask
end


""" Given x and y (2D meshgrid arrays), return a bandpass mask that is damped by an envelope function. """
function create_damped_mask(x::AbstractArray, y::AbstractArray, low, high, envelope)
    r = @. sqrt(x^2 + y^2)
    mask = Float64.(low .< r .< high)
    @. mask *= envelope(r - low)
    return mask
end

""" Subtract noise background and calculate approximate defocus for optimization.

Returns a reference spectrum (fftshifted), rough defocus estimate, a damped mask for simulation,
and meshgrid arrays kx & ky of spatial frequencies.

`Cs`: spherical abberation in millimeters.
`λ`: electron wavelength in picometers.
`pixelsize`: real space pixel size in Angstroms.
`min_resolution`: in Angstroms.
`max_resolution`: in Angstroms.
"""
function preprocess_spectrum(
    spectrum::AbstractArray,
    Cs,
    λ,
    pixelsize,
    min_resolution,
    max_resolution,
)
    N = size(spectrum)[2]
    low_res_idx = round(Int, N * pixelsize / min_resolution)
    high_res_idx = round(Int, N * pixelsize / max_resolution)

    mid = ceil(Int, (N + 1) / 2)
    search_range2D = Range2D(1:high_res_idx, mid-(high_res_idx-1):mid+(high_res_idx-1))
    shifted_spectrum = fftshift(spectrum, 2)
    search_spectrum = shifted_spectrum[search_range2D.x, search_range2D.y]

    # calculate background and envelope functions from 1D profile
    ctfProfile = profile1D(real(spectrum))
    search_profile = ctfProfile[low_res_idx:high_res_idx]
    background_func = get_background(search_profile)
    subtracted_profile = search_profile - background_func.(1:length(search_profile))
    envelope = get_envelope_function(subtracted_profile)

    # estimate approximate defocus
    second_peak_idx = low_res_idx + envelope.itp.knots[1][2] - 1
    k_peak = second_peak_idx / (N * pixelsize)
    z_estimate = (1.5 + 0.5 * Cs * λ^3 * k_peak^4 * 10) / (λ * k_peak^2 * 100)

    kx = (
        fftshift(fftfreq(N, 1 / pixelsize))[search_range2D.y]' .*
        ones(length(search_range2D.x))
    )
    ky = ones(length(search_range2D.y))' .* rfftfreq(N, 1 / pixelsize)[search_range2D.x]
    x = round.(Int, kx * N * pixelsize)
    y = round.(Int, ky * N * pixelsize)
    r = @. sqrt(x^2 + y^2)

    background_spectrum = background_func.(r .- low_res_idx)
    subtracted_spectrum = search_spectrum - background_spectrum
    binary_mask = create_binary_mask(x, y, low_res_idx, high_res_idx)
    reference_spectrum = subtracted_spectrum .* binary_mask
    damped_mask = create_damped_mask(x, y, low_res_idx, high_res_idx, envelope)

    return (reference_spectrum, z_estimate, damped_mask, kx, ky, envelope, low_res_idx)
end


""" Calculate ctf in place

`A`: array to store result.
`z1`, `z2`: defocus values in microns.
`θ`: astigmatism angle in degrees.
`ϕ`: phase shift in degrees.
`kx`, `ky`: meshgrid spatial frequencies.
`λ`: electron wavelength in picometers.
`Cs`: spherical abberation in millimeters.
`damped_mask`: mask created by preprocessing
"""
function ctf!(
    A::AbstractArray,
    z1,
    z2,
    θ,
    ϕ,
    kx::AbstractArray,
    ky::AbstractArray,
    λ,
    Cs,
    amplitude_contrast,
    damped_mask::AbstractArray,
)
    c1 = -50 * λ
    c2 = 2.5 * Cs * λ^3
    c3 = ϕ / 360
    c4 = pi - asin(amplitude_contrast)
    semiaxes = [
        cos(θ * pi / 180) -sin(θ * pi / 180)
        sin(θ * pi / 180) cos(θ * pi / 180)
    ]
    Z = semiaxes * [z1 0; 0 z2] * semiaxes'
    for i in eachindex(A)
        quad_form = kx[i]^2 * Z[1, 1] + 2 * kx[i] * ky[i] * Z[1, 2] + ky[i]^2 * Z[2, 2]
        W = c1 * quad_form + c2 * (kx[i]^2 * ky[i]^2)^2 + c3
        A[i] = 2pi * W + c4
    end
    # calculating sin inside the above loop is slower?
    @. A = abs(sin(A)) * damped_mask
    return nothing
end


""" Returns the cross correlation between a reference power spectrum and an estimate from ctf!

Same arguments as ctf! with an additional `reference_spectrum` image.
"""
function f!(
    A::AbstractArray,
    z1,
    z2,
    θ,
    ϕ,
    kx::AbstractArray,
    ky::AbstractArray,
    λ,
    Cs,
    amplitude_contrast,
    damped_mask::AbstractArray,
    reference_spectrum::AbstractArray,
)
    ctf!(A, z1, z2, θ, ϕ, kx, ky, λ, Cs, amplitude_contrast, damped_mask)
    cross_correlation =
        dot(real(reference_spectrum), A) / (norm(A) * norm(real(reference_spectrum)))
    return cross_correlation
end


struct SearchResult
    z1::Float64
    z2::Float64
    θ::Float64
    ϕ::Float64
    correlation::Float64
end


function Base.show(io::IO, x::SearchResult)
    print(io, "$(x.z1),$(x.z2),$(x.θ),$(x.ϕ),$(x.correlation)")
end


""" Maximizes the cross correlation to determine ctf parameters.

`input_image`: path to mrc.
`pixelsize`: pixel size in Angstroms.
`voltage`: accelerating voltage in kV.
`Cs`: spherical abberation in millimeters.
`min_resolution`: min resolution cutoff.
`max_resolution`: max resolution cutoff.
`amplitude_contrast`: fraction between 0 and 1.
`do_phase_search`: boolean to also search for phase shift.
"""
function find_ctf(
    input_image,
    pixelsize,
    voltage,
    Cs,
    min_resolution,
    max_resolution,
    amplitude_contrast,
    do_phase_search::Bool,
)
    img = real.(load(input_image))
    spectrum = get_spectrum(img)

    λ = wavelength_from_voltage(voltage)
    reference_spectrum, z_estimate, damped_mask, kx, ky =
        preprocess_spectrum(spectrum, Cs, λ, pixelsize, min_resolution, max_resolution)
    A = similar(reference_spectrum, Float64)

    if do_phase_search
        opt = Opt(:LN_COBYLA, 4)
        opt.lower_bounds = [0.5, 0.5, -90, 0.0]
        opt.upper_bounds = [5.0, 5.0, 90, 180.0]
        opt.max_objective =
            (p, grad) -> f!(
                A,
                p[1],
                p[2],
                p[3],
                p[4],
                kx,
                ky,
                λ,
                Cs,
                amplitude_contrast,
                damped_mask,
                reference_spectrum,
            )
        p0 = [z_estimate, z_estimate, 0, 0]
    else
        opt = Opt(:LN_COBYLA, 3)
        opt.lower_bounds = [0.5, 0.5, -90]
        opt.upper_bounds = [5.0, 5.0, 90]
        opt.max_objective =
            (p, grad) -> f!(
                A,
                p[1],
                p[2],
                p[3],
                0,
                kx,
                ky,
                λ,
                Cs,
                amplitude_contrast,
                damped_mask,
                reference_spectrum,
            )
        p0 = [z_estimate, z_estimate, 0]
    end
    opt.xtol_rel = 1e-4
    inequality_constraint!(opt, (p, grad) -> p[2] - p[1])

    (correlation, params, ret) = optimize(opt, p0)
    numevals = opt.numevals # the number of function evaluations
    if do_phase_search
        search_result = SearchResult(round.((params[1], params[2], params[3], params[4], correlation), digits=4)...)
    else
        search_result = SearchResult(round.((params[1], params[2], params[3], 0, correlation), digits=4)...)
    end

    return reference_spectrum, A, search_result
end
