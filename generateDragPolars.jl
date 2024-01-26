import FLOWVLM as vlm
using CSV
using DataFrames
using Plots
plotly()

file = CSV.File(open("/home/bronbergcfd/Dropbox/bronberg/projects/183-nobleProp/foils/propdata.csv"))
accountForCompressible = true
temp = 313.0
RPM = 2350
Rtip = 0.9099
roR = [0.1, 0.368421052631579, 0.6842105263157894, 0.8947368421052632, 1.0]
coR = [0.0651162790697674, 0.0651162790697674, 0.0651162790697674, 0.0651162790697674, 0.0433953488372093]
AR = coR./roR
name = "rc410TE"
airfoil = "/home/theaero/Dropbox/VPMData/database/airfoils/rc410TE0_4.csv"
CDmax = 1.3
#speed of sound
c = sqrt(1.4*278*temp)
#viscousity
mu_ = 1.85508e-5 
T_  = 15+273.15
S = 113
mu = mu_*(temp/T_)^(3/2)*(T_+S)/(temp+S)
println("mu: $mu")
#rho
rho = 1.13
#omega
omega = RPM/60*2pi
plt1=plot()
plt2=plot()
plt3=plot()
for i in 1:length(roR)
    ar = AR[i]
    r_over_R = roR[i]
    c_over_R = coR[i]
    r = r_over_R*Rtip
    cord = c_over_R*Rtip
    V = omega*r
    Re = rho*cord*V/mu
    Mach = accountForCompressible*V/c

    df = DataFrame(CSV.File(airfoil))
    x = df."x/c"
    y = df."y/c"

    polar = vlm.ap.runXFOIL(x, y, Re;
                            alphas=range(-25,25,101),
                            verbose=true, Mach=Mach,
                            iter=200, ncrit=9)

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

    polar = vlm.ap.extrapolate(polar, CDmax, AR=ar)
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
    vlm.ap.save_polar2(polar, "polars/$(name)RE$(string(Int(round(Re, digits=-3))))Mach$(string(Int(round(Mach*100, digits=0)))).csv")

end
title!(plt1, "CL")
title!(plt2, "CD")
title!(plt3, "CM")
nrows = 2
ncols = 2
plt = plot(plt1, plt2,plt3, layout=(nrows,ncols), size=(ncols*700.0, nrows*600.0), legend=true)
display(plt)
#vlm.ap.get_xsepup(polar)
#vlm.ap.get_xseplo(polar)
vlm.ap.get_(polar)