# using ADIOS2
using CSV
using CairoMakie
using SixelTerm
# using openPMD

# function plot_adios2()
#     filename = "/Users/eschnett/simulations/run_cake_coord_self_tests/output-0000/run_cake_coord_self_tests/run_cake_coord_self_tests.it00000000.bp5"
#     file = adios_open_serial(filename, mode_read)
#     return
# end

# function plot_openpmd()
#     filename = "/Users/eschnett/simulations/run_cake_coord_self_tests/output-0000/run_cake_coord_self_tests/run_cake_coord_self_tests.it00000000.bp5"
#     series = Series(filename, ACCESS_READ_ONLY)
#     return
# end

function plot_csv()
    filename1 = "/Users/eschnett/simulations/run_cake_coord_self_tests/output-0000/run_cake_coord_self_tests/coordinates-vertex_coords.it000000.p0000.tsv"
    filename2 = "/Users/eschnett/simulations/run_cake_coord_self_tests/output-0000/run_cake_coord_self_tests/testmultipatch-multipatch_test_gfs.it000000.p0000.tsv"
    file1 = CSV.File(filename1)
    file2 = CSV.File(filename2)

    npoints = length(file1)
    @assert length(file2) == npoints

    ps = file1[Symbol("3:patch")]
    is = file1[Symbol("6:i")]
    js = file1[Symbol("7:j")]
    ks = file1[Symbol("8:k")]
    as = file1[Symbol("9:x")]
    bs = file1[Symbol("10:y")]
    cs = file1[Symbol("11:z")]
    xs = file1[Symbol("12:vcoordx")]
    ys = file1[Symbol("13:vcoordy")]
    zs = file1[Symbol("14:vcoordz")]
    gf = file2[Symbol("12:test_gf")]

    @assert (is[begin] + is[end]) % 2 == 0
    mid_i = (is[begin] + is[end]) ÷ 2
    @assert (js[begin] + js[end]) % 2 == 0
    mid_j = (js[begin] + js[end]) ÷ 2
    @assert (ks[begin] + ks[end]) % 2 == 0
    mid_k = (ks[begin] + ks[end]) ÷ 2

    # Patches:
    #     0: center
    #     1: +x
    #     2: -x
    #     3: +y
    #     4: -y
    #     5: +z
    #     6: -z
    # Local coordinates:
    #     0: xyz
    #     1: θϕr

    # inds = Int[n for n in 1:npoints if ps[n] == 1 && isapprox(zs[n], 0; atol=sqrt(eps()))]

    indsxy = Dict{Int,AbstractVector{Int}}()
    indsxy[0] = Int[n for n in 1:npoints if ps[n] == 0 && ks[n] == mid_k]
    indsxy[1] = Int[n for n in 1:npoints if ps[n] == 1 && is[n] == mid_i]
    indsxy[2] = Int[n for n in 1:npoints if ps[n] == 2 && is[n] == mid_i]
    indsxy[3] = Int[n for n in 1:npoints if ps[n] == 3 && is[n] == mid_i]
    indsxy[4] = Int[n for n in 1:npoints if ps[n] == 4 && is[n] == mid_i]

    fig1 = Figure(; resolution=(1280, 960))
    ax = Axis(fig1[1, 1]; title="Coordinates", xlabel="x", ylabel="y")
    rmax = 5.5
    xlims!(-rmax, +rmax)
    ylims!(-rmax, +rmax)
    for patch in 0:4
        scatter!((@view xs[indsxy[patch]]), (@view ys[indsxy[patch]]); colormap=:plasma, markersize=4, label="patch $patch")
    end
    axislegend()
    colsize!(fig1.layout, 1, Aspect(1, 1.0))
    # rowsize!(fig1.layout, 1, Aspect(1, 1.0))

    fig2 = Figure(; resolution=(1280, 960))
    ax = Axis(fig2[1, 1]; title="Test grid function", xlabel="x", ylabel="y")
    rmax = 5.5
    xlims!(-rmax, +rmax)
    ylims!(-rmax, +rmax)
    for patch in 0:4
        tr = tricontourf!(fig2[1, 1], xs[indsxy[patch]], ys[indsxy[patch]], gf[indsxy[patch]]; colormap=:plasma)
        Colorbar(fig2[1, 2], tr)
    end
    # axislegend()
    colsize!(fig2.layout, 1, Aspect(1, 1.0))
    # rowsize!(fig2.layout, 1, Aspect(1, 1.0))

    return fig1, fig2
end

nothing
