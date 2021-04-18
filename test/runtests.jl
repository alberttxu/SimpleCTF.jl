using SimpleCTF
using Test


function is_imod_clip_installed()
    try
        read(`which clip`)
    catch
        return false
    end
    return true
end


@testset "SimpleCTF.jl" begin
    # Write your tests here.
    @test is_imod_clip_installed()
    @test isapprox(SimpleCTF.wavelength_from_voltage(200), 2.5079, atol=1e-4)
    @test isapprox(SimpleCTF.wavelength_from_voltage(300), 1.9687, atol=1e-4)
end
