import FLOWVLM as vlm
using CSV
using DataFrames
using Plots
plotly()

file = CSV.File(open("/home/aero/Dropbox/bronberg/projects/183-nobleProp/foils/propdata.csv"))

CDmax = 1.3

plt1=plot()
plt2=plot()
plt3=plot()
plt4=plot()
for i in 1:length(file["foil"])
    
    r_over_R = file["r/R"][i]
    c_over_R = file["chord/R"][i]
    ar = c_over_R/r_over_R
    #r = r_over_R*Rtip
    #cord = c_over_R*Rtip
    #V = omega*r
    Re = file["Re"][i]
    Mach = file["Mach number"][i]
    foil = "/home/aero/Dropbox/bronberg/projects/183-nobleProp/foils/$(file["foil"][i])"
    df = DataFrame(CSV.File(foil;header=false,skipto=2,delim="  "))
    x = df."Column1"
    y = df."Column2"

    polar = vlm.ap.runXFOIL(x, y, Re;
                            alphas=range(-25,25,101),
                            verbose=true, Mach=Mach,
                            iter=500, ncrit=9)

    #=
    alpha, CL = vlm.ap.get_cl(polar)
    alpha, CD = vlm.ap.get_cd(polar)
    alpha, CM = vlm.ap.get_cm(polar)
    plot!(plt1, alpha, CL, label = "out")
    plot!(plt2, alpha, CD, label = "")
    plot!(plt3, alpha, CM, label = "")
    =#

    polar = vlm.ap.correction3D(polar, r_over_R, c_over_R, nothing)
    alpha, CL = vlm.ap.get_cl(polar)
    alpha, CD = vlm.ap.get_cd(polar)
    alpha, CM = vlm.ap.get_cm(polar)
    plot!(plt1, alpha, CL, label = "RE$(Re)Mach$(Mach)")
    plot!(plt2, alpha, CD, label = "")
    plot!(plt3, alpha, CM, label = "")
    plot!(plt4, CD, CL, label = "")

    #polar = vlm.ap.extrapolate(polar, CDmax, AR=ar)
    #=
    alpha, CL = vlm.ap.get_cl(polar)
    alpha, CD = vlm.ap.get_cd(polar)
    alpha, CM = vlm.ap.get_cm(polar)
    plot!(plt1, alpha, CL, label = "360 extrapolated")
    plot!(plt2, alpha, CD, label = "")
    plot!(plt3, alpha, CM, label = "")
    =#
    #=
    this_polar = ap.extrapolate(this_polar, CDmax, AR=AR)

    polar = vlm.ap.runXFOIL(x, y, this_Re;
                            alphas=alphas,
                            verbose=verbose_xfoil, Mach=this_Ma,
                            iter=100, ncrit=ncrit)
    if verbose; println(" done."); end;

    if save_polars != nothing
        if !(isdir(save_polars)); mkdir(save_polars); end;
        vlm.ap.save_polar2(polar, save_polar_pref*"-sec$(rfli)-Re$(ceil(Int, this_Re))"; path=save_polars)
    end


    function runXFOIL(xs, ys, Re; alphas=range(-10, stop=25, step=0.5), Mach=0.0,
        iter::Int=100, verbose=false,
        npanels=160, ncrit=9, optargs...)
    =#
    vlm.ap.save_polar2(polar, "polars/$(file["foil"][i])RE$(string(Int(round(Re, digits=-3))))Mach$(string(Int(round(Mach*100, digits=0)))).csv")

end
title!(plt1, "CL")
title!(plt2, "CD")
title!(plt3, "CM")
title!(plt4, "CL vs CD")
nrows = 2
ncols = 2
plt = plot(plt1, plt2,plt3,plt4, layout=(nrows,ncols), size=(ncols*700.0, nrows*600.0), legend=true)
display(plt)
#vlm.ap.get_xsepup(polar)
#vlm.ap.get_xseplo(polar)
#vlm.ap.get_(polar)