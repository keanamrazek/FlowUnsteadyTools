#=
Description:
case script for rotor only systems

Author: Keana Mrazek

=#
caseName = "VPMCaseP$var"

bulkDataSavePath = "/home/aero/Desktop/tempVPMFilesKeana/V$airSpeed/"*caseName

data_path        = "/home/aero/Dropbox/VPMData/database" #/mnt/c/Dropbox/VPMData/database/"

savePath = false

runParaview=false

plotResults=true

runXFoil = false

save_wopwopin = false

################################################################
## Case description
fidelityLevel = "Low"  # Low, Mid, High

includeViscous = "Inviscid" # "Inviscid", "CoreSpreading", "ParticleStrengthExchange"

magVinf         = 1e-8+airSpeed
magVvehicle     = 1e-8 

AOA = 0.0

sideSlip = 0.0

temperature = 15.0 #deg C ground temperature

altitude = 1368.0 #m
ISA = 25

rhoOverride = 0.984096785411496
muOverride = 1.904993200205812e-5
temperatureOverride = 10

nrevs    = 8 # Number of revolutions in simulation

domainLimits = 5.0^2 #simulation domain radius

gammaCrit = 0.005^2  #minimum vorticity

gammaCritup = 13.0^2 #maximum vorticity

wakeTreatmentStep = 20 #wake treatment interval

domainCentre = [0.0, 0.0, 0.0]

translateVehicle = [0.0, 0.0, 0.0]                                 # New position

Vinf(X, t)      = t==0 ? magVvehicle*[1,0,0] : magVinf*[1,0,0] # Freestream function

# Non-dimensional translational velocity of vehicle over time
Vvehicle(t) = [-1, 0, 0]        # <---- Vehicle is traveling in the -x direction

# Angle of the vehicle over time
anglevehicle(t) = zeros(3)

# RPM control input over time (RPM over `RPMref`)
RPMcontrol(t) = 1.0

angles = ()                                 # Angle of each tilting system (none)
RPMs = (RPMcontrol, )                       # RPM of each rotor system

Vref = magVvehicle                          # Reference velocity to scale maneuver by
magVref         = sqrt(magVinf^2 + magVvehicle^2) # (m/s) reference velocity
RPMref = 2350.0                                # Reference RPM to scale maneuver by
################################################################
#Rotor system description

ROut            = [0.90990415] #m

RInner          = [0.10] #m

numberOfBlades  = [5]

ypos            = [0.0]

xpos            = [0.0]

zpos            = [0.0] 

roll            = [0.0] # of rotor
pitch           = [180.0]
yaw             = [0.0]


CWs             = [false]

bladeFile      = ["nobleProp_blade.csv"]       # Blade geometry

collective      = [11.3519405652475+var]

coningAngle     = [0.0]
ROut = ROut.*cosd.(coningAngle)
RInner = RInner.*cosd.(coningAngle)

rotorRPM        = [2350]

numberOfRotors = length(ROut)

ncrit  = 9 # Turbulence criterion for XFOIL

################################################################
# imports

import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm        
using CSV
using DataFrames
using Plots
using FFTW

################################################################
# Pre proccessing

if fidelityLevel=="High"
    n_rotor         = 50                  # Number of blade elements per blade
    r_rotor         = 1/15                # Geometric expansion of elements
    nsteps_per_rev  = 360
    p_per_step      = 2                   # Sheds per time step
    vpm_integration = vpm.rungekutta3
elseif fidelityLevel=="Mid"
    n_rotor         = 20                  # Number of blade elements per blade
    r_rotor         = 1/13                # Geometric expansion of elements
    nsteps_per_rev  = 72
    p_per_step      = 4                   # Sheds per time step
    vpm_integration = vpm.rungekutta3
elseif fidelityLevel=="Low"
    n_rotor         = 15                  # Number of blade elements per blade
    r_rotor         = 1/11                # Geometric expansion of elements
    nsteps_per_rev  = 45
    p_per_step      = 4                   # Sheds per time step
    vpm_integration = vpm.euler
end

# time parameters nsteps/nsteps_per_rev / (RPM/60) 
nsteps = nrevs*nsteps_per_rev
ttot            = nsteps/nsteps_per_rev / (rotorRPM[1]/60)
max_particles   = Int(round(sum(numberOfBlades)*n_rotor*p_per_step*nsteps*2+1, digits=0))#nrotors*((2*n_rotor+1)*B)*nsteps*p_per_step + 1 # Maximum number of particles
# density altitude
TAlt = 15.04 - 0.00659*altitude
PAlt = 101.29*((TAlt+273.15)/288.08)^5.256
rho = PAlt/(0.2869*(TAlt+ISA+273.15))

#viscousity
mu_ = 1.85508e-5 
T_  = 15+273.15
S = 113
mu = mu_*((TAlt+273.15)/T_)^(3/2)*(T_+S)/(TAlt+273.15+ S)


if rhoOverride !== nothing
    rho = rhoOverride
end
if temperatureOverride !== nothing
    TAlt = temperatureOverride
    mu = mu_*((TAlt+273.15)/T_)^(3/2)*(T_+S)/(TAlt+273.15+ S)
end
if muOverride !== nothing
    mu = muOverride
end



#speed of sound
c = sqrt(1.4*286.9*(TAlt+273.15))

#
ReD             = ((2*pi*rotorRPM./60).*ROut* rho/mu * 2).*ROut      # Diameter-based rotor Reynolds number
Matip           = (2*pi*rotorRPM./60/c).*ROut   # Tip Mach number

println("""
    Vref:   $(round(magVref, digits=1)) m/s
    Temperature:  $(round(TAlt, digits=1))
    mu:     $(mu)
    rho:    $(round(rho, digits=3))
    RPM:    $(rotorRPM)
    Matip:  $(round.(Matip, digits=3))
    ReD:    $(round.(ReD, digits=0))
    Max particles: $(max_particles)
""")
################################################################
#solver parameters
VehicleType     = uns.UVLMVehicle
const_solution = false

# VPM particle shedding
shed_starting   = true                      # Whether to shed starting vortex
shed_unsteady   = true                      # Whether to shed vorticity from unsteady loading
unsteady_shedcrit = 0.001                   # Shed unsteady loading whenever circulation
                                            #  fluctuates by more than this ratio
# Regularization
sigma_rotor_surf= ROut[1]/100                      # Rotor-on-VPM smoothing radius
lambda_vpm      = 2.125                     # VPM core overlap
                                            # VPM smoothing radius
sigma_vpm_overwrite = lambda_vpm * 2*pi*ROut[1]/(nsteps_per_rev*p_per_step) #can set to "nothing" then, lambda_vpm will be used
sigmafactor_vpmonvlm= 1                     # Shrink particles by this factor when
                                            #  calculating VPM-on-VLM/Rotor induced velocities
#Viscousity dissapation scheme
if includeViscous == "Inviscid"
    vpm_viscous     = vpm.Inviscid()
elseif includeViscous == "CoreSpreading"
    vpm_viscous     = vpm.CoreSpreading(-1, -1, vpm.zeta_fmm; beta=100.0, itmax=20, tol=1e-1)
#elseif includeViscous == "ParticleStrengthExchange"
end

# VPM LES subfilter-scale model (turbulent decay of the wake (turbulent diffusion))
#vpm_SFS         = vpm.SFS_none              
# vpm_SFS       = vpm.SFS_Cd_twolevel_nobackscatter
# vpm_SFS       = vpm.SFS_Cd_threelevel_nobackscatter
vpm_SFS       = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive;
                                   alpha=0.999, maxC=1.0,
                                   clippings=[vpm.clipping_backscatter])
# vpm_SFS       = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive;
#                                   alpha=0.999, rlxf=0.005, minC=0, maxC=1
#                                   clippings=[vpm.clipping_backscatter],
#                                   controls=[vpm.control_sigmasensor],
#                                   )

# Rotor solver
vlm_rlx         = 0.5                       # VLM relaxation <-- this also applied to rotors
#hubtiploss_correction = ((2, 1, 0.9, 0.05), (4, 1, 1, 0.05)) # Hub and tip correction # Hub and tip correction
hubtiploss_correction       = vlm.hubtiploss_correction_prandtl # Hub and tip correction
################################################################
# fountain effect
suppress_fountain   = false

# Supress wake shedding on blade elements inboard of this r/R radial station
no_shedding_Rthreshold = suppress_fountain ? 0.35 : 0.0

# Supress wake shedding for this many time steps
no_shedding_nstepsthreshold = 3*nsteps_per_rev

omit_shedding = []          # Index of blade elements to supress wake shedding

# Function to suppress or activate wake shedding
function wake_treatment_supress(sim, args...; optargs...)

    # Case: start of simulation -> suppress shedding
    if sim.nt == 1

        # Identify blade elements on which to suppress shedding
        for i in 1:vlm.get_m(rotor)
            HS = vlm.getHorseshoe(rotor, i)
            CP = HS[5]

            if uns.vlm.norm(CP - vlm._get_O(rotor)) <= no_shedding_Rthreshold*R
                push!(omit_shedding, i)
            end
        end
    end

    # Case: sufficient time steps -> enable shedding
    if sim.nt == no_shedding_nstepsthreshold

        # Flag to stop suppressing
        omit_shedding .= -1

    end

    return false
end

################################################################
# build vehicle

#rotors
println("Generating rotors...")
rotors = vlm.Rotor[]
J = magVref/(rotorRPM[1]/60)/ROut[1]
for ri in 1:numberOfRotors
    df = DataFrame(CSV.File(data_path*"/rotors/"*bladeFile[ri]))

    cordDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[1]))
    pitchDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[2]))
    sweepDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[3]))
    heightDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[4]))
    airfoilDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[5]))

    sweepDistNew  = cosd(collective[ri])*sweepDist."y/R" - sind(collective[ri])*heightDist."z/R"
    heightDistNew = sind(collective[ri])*sweepDist."y/R" + cosd(collective[ri])*heightDist."z/R"

    
    cordDist = [cordDist."r/R" cordDist."c/R"]
    pitchDist = [pitchDist."r/R" pitchDist."twist (deg)"]
    sweepDist = [sweepDist."r/R" sweepDistNew]
    heightDist = [heightDist."r/R" heightDistNew+sind(coningAngle[ri])*heightDistNew]
    airfoil_contours::Vector{Tuple{Float64, String, String}} = [(airfoilDist."r/R"[1], airfoilDist."Contour file"[1], airfoilDist."Aero file"[1])]

    for i in 1:length(airfoilDist."r/R")-1
        println((airfoilDist."r/R"[i+1], airfoilDist."Contour file"[i+1], airfoilDist."Aero file"[i+1]))
        push!(airfoil_contours, (airfoilDist."r/R"[i+1], airfoilDist."Contour file"[i+1], airfoilDist."Aero file"[i+1]))
    end



    rotor = uns.generate_rotor(ROut[ri], RInner[ri], numberOfBlades[ri], cordDist, pitchDist, sweepDist, heightDist, airfoil_contours;
                                        pitch=collective[ri],
                                        n=n_rotor, CW=CWs[ri], blade_r=r_rotor,
                                        altReD=[rotorRPM[ri], J, mu/rho],
                                        xfoil=runXFoil,
                                        read_polar=vlm.ap.read_polar2,
                                        ncrit=ncrit,
                                        data_path=data_path,
                                        verbose=true,
                                        verbose_xfoil=true,
                                        plot_disc=true,
                                        ReD=ReD[ri],
                                        Matip=0,  #this is done such that compressibility is taken care of by the solver, another option is to take compressibility into account with the drag polars
                                        spline_s=0.0,
                                        );

    # Determine position along wing LE
    y = ypos[ri]
    x = xpos[ri]
    z = zpos[ri]

    # Account for angle of attack of wing



    # Translate rotor to its position along wing
    O = [x, y, z]                                       # New position
    Oaxis = uns.gt.rotation_matrix2(roll[ri], pitch[ri], yaw[ri])          # New orientation  Roll(x), pitch(y), yaw(z)
    vlm.setcoordsystem(rotor, O, Oaxis)

    push!(rotors, rotor)
end



#vehicle
println("Generating vehicle...")
system = vlm.WingSystem()
rotor_systems = (rotors, );
vlm_system = vlm.WingSystem()
wake_system = vlm.WingSystem() 

for (ri, rotor) in enumerate(rotors)
    vlm.addwing(system, "Rotor$(ri)", rotor)
    vlm.addwing(wake_system, "Rotor$(ri)", rotor)
end
#=
if VehicleType != uns.QVLMVehicle
    for (ri, rotor) in enumerate(rotors)
        vlm.addwing(wake_system, "Rotor$(ri)", rotor)
    end
end
=#
# Pitch vehicle to its angle of attack
Oaxis = uns.gt.rotation_matrix2(0, AOA, 0)         # New orientation
vlm.setcoordsystem(system, translateVehicle, Oaxis)


vehicle = VehicleType(   system;
                            vlm_system=vlm_system,
                            rotor_systems=rotor_systems,
                            wake_system=wake_system
                         );
############################manoeuvre definition####################################


maneuver = uns.KinematicManeuver(angles, RPMs, Vvehicle, anglevehicle)
############################simulation definition####################################
Vinit = Vref*Vvehicle(0)                    # Initial vehicle velocity
Winit = pi/180*(anglevehicle(1e-6) - anglevehicle(0))/(1e-6*ttot)  # Initial angular velocity

simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                    Vinit=Vinit, Winit=Winit)
#################################monitors and wake treatment###############################
#Wake treatment

wake_treatment1 = uns.remove_particles_strength(gammaCrit, gammaCritup)
wake_treatment2 = uns.remove_particles_sphere(domainLimits, wakeTreatmentStep, Xoff=domainCentre)

# Generate rotors monitor
monitors = uns.generate_monitor_rotors(rotors, J, rho, RPMref, nsteps;
                                            t_scale=RPMref/60,        # Scaling factor for time in plots
                                            t_lbl="Revolutions",   # Label for time axis
                                            save_path=bulkDataSavePath,
                                            run_name=caseName*"-rotors",
                                            figname="rotors monitor",
                                            )
# Concatenate monitors

extra_runtime_function = uns.concatenate(monitors, wake_treatment1, wake_treatment2)

#################################run simulation###############################
# run model
println("""
    Vref:   $(round(magVref, digits=1)) m/s
    RPM:    $(rotorRPM)
    Matip:  $(round.(Matip, digits=3))
    ReD:    $(round.(ReD, digits=0))
    Max particles: $(max_particles)
""")
println("Running simulation...")

uns.run_simulation(

    simulation,                    # Simulation object
    nsteps;                        # Total time steps in simulation

    # -------- SIMULATION OPTIONS -----------------------------------------
    Vinf            = Vinf,#Vinf, # Freestream velocity
    sound_spd       = c,# (m/s) speed of sound, apparently adding it here and in xfoil is double accounting for compressibility and shouldnt be done
    rho             = rho,            # (kg/m^3) air density
    mu              = mu,          # (Pa*s) air dynamic viscosity
    #tquit           = Inf,              # (s) force quit the simulation at this time
    #rand_RPM        = false,            # (experimental) randomize RPM fluctuations

    extra_runtime_function = extra_runtime_function,

    # -------- SOLVERS OPTIONS --------------------------------------------
    # Vortex particle method
    max_particles   = max_particles,         # Maximum number of particles
    #max_static_particles = nothing,     # Maximum number of static particles (use `nothing` to automatically estimate it)
    p_per_step      = p_per_step,                # Particle sheds per time step
    #vpm_formulation = vpm.rVPM,         # VPM formulation (`vpm.rVPM` or `vpm.cVPM`)
    #vpm_kernel      = vpm.gaussianerf,  # VPM kernel (`vpm.gaussianerf` or `vpm.winckelmans`)
    #vpm_UJ          = vpm.UJ_fmm,       # VPM particle-to-particle interaction scheme (`vpm.UJ_fmm` or `vpm.UJ_direct`)
    vpm_SFS         = vpm_SFS,          # VPM LES subfilter-scale model (`SFS_none`, `SFS_Cd_threelevel_nobackscatter`, `SFS_Cd_twolevel_nobackscatter`, or `SFS_Cs_nobackscatter`)
    vpm_integration = vpm_integration,  # VPM time integration scheme (`vpm.euler` or `vpm.rungekutta3`)
    #vpm_transposed  = true,             # VPM transposed stretching scheme
    #vpm_viscous     = vpm_viscous,      # VPM viscous diffusion scheme (`vpm.Inviscid()`, `vpm.CoreSpreading(nu, sgm0, zeta)`, or `vpm.ParticleStrengthExchange(nu)`)
    #vpm_fmm         = vpm.FMM(; p=4, ncrit=50, theta=0.4, phi=0.5), # VPM's FMM settings
    #vpm_relaxation  = vpm.pedrizzetti,  # VPM relaxation scheme (`vpm.norelaxation`, `vpm.correctedpedrizzetti`, or `vpm.pedrizzetti`)
    #vpm_surface     = true,             # Whether to include surfaces in the VPM through ASM/ALM

    # Actuator surface/line model (ASM/ALM): VLM and blade elements
    #vlm_vortexsheet = false,            # Whether to spread surface circulation as a vortex sheet in the VPM (turns ASM on; ALM if false)
    #vlm_vortexsheet_overlap     = 2.125,# Overlap of particles that make the vortex sheet
    #vlm_vortexsheet_distribution= uns.g_pressure, # Vorticity distribution of vortex sheet (`g_uniform`, `g_linear`, or `g_pressure`)
    #vlm_vortexsheet_sigma_tbv   = nothing, # Size of particles in trailing bound vortices (defaults to `sigma_vlm_surf` if not given)
    vlm_rlx         = vlm_rlx,               # VLM relaxation (>0.9 can cause divergence, <0.2 slows simulation too much, deactivated with <0)
    #vlm_init        = false,            # Initialize the first step with the VLM semi-infinite wake solution
    hubtiploss_correction = hubtiploss_correction, # Hub and tip loss correction of rotors (ignored in quasi-steady solver)

    # Wake shedding
    #wake_coupled        = true,         # Couple VPM wake -> VLM solution
    shed_unsteady       = shed_unsteady,         # Whether to shed vorticity from unsteady loading
    unsteady_shedcrit   = unsteady_shedcrit,         # Criterion for unsteady-loading shedding
    shed_starting       = shed_starting,        # Whether to shed starting vortex (only when `shed_unsteady=true`)
    #shed_boundarylayer  = false,        # (experimental) whether to shed vorticity from boundary layer of surfaces
    #boundarylayer_prescribedCd = 0.1,   # (experimental) prescribed Cd for boundary layer shedding used for wings
    #boundarylayer_d     = 0.0,          # (experimental) dipole width for boundary layer shedding
    #omit_shedding       = [],           # Indices of elements in `sim.vehicle.wake_system` on which omit shedding VPM particles

    # Regularization of solvers
    #sigma_vlm_solver    = -1,           # Regularization of VLM solver (internal VLM-on-VLM)
    sigma_vlm_surf      = 1,           # (REQUIRED!) Size of embedded particles in ASM/ALM wing surfaces (for VLM-on-VPM and VLM-on-Rotor)
    sigma_rotor_surf    = sigma_rotor_surf,           # (REQUIRED!) Size of embedded particles in ALM blade surfaces (for Rotor-on-VPM, Rotor-on-VLM, and Rotor-on-Rotor)
    sigmafactor_vpm     = lambda_vpm,          # Core overlap of wake particles
    sigmafactor_vpmonvlm = sigmafactor_vpmonvlm,           # (experimental) shrinks the particles by this factor when calculating VPM-on-VLM/Rotor induced velocities
    sigma_vpm_overwrite = sigma_vpm_overwrite,      # Overwrite core size of wake to this value (ignoring `sigmafactor_vpm`)

    # -------- RESTART OPTIONS --------------------------------------------
    #restart_vpmfile     = nothing,      # VPM restart file to restart simulation

    # -------- OUTPUT OPTIONS ---------------------------------------------
    save_path       = bulkDataSavePath, # Where to save simulation
    run_name        = caseName,# Suffix of output files
    #create_savepath = true,             # Whether to create `save_path`
    prompt          = false,             # Whether to prompt the user
    #verbose         = true,             # Enable verbose
    #v_lvl           = 0,                # Indentation level of verbose
    #verbose_nsteps  = 10,               # Verbose every this many steps
    #raisewarnings   = true,             # Whether to raise warnings
    #debug           = true,            # Output extra states for debugging
    #nsteps_save     = 1,                # Save vtks every this many steps
    #nsteps_restart  = -1,               # Save jlds every this many steps (restart files)
    #save_code       = bulkDataSavePath, # Copy the source code in this path to `save_path`
    #save_horseshoes = false,            # Whether to output VLM horseshoes in VTKs
    #save_static_particles = true,       # Whether to save ASM/ALM embedded particles
    save_wopwopin   = save_wopwopin,            # Generate input files for PSU-WOPWOP

)
println("Complete:)")
################################################################
# post processing

if plotResults


    plotly()

    function hanning(x)
        w = Vector(LinRange(-0.5pi, 0.5pi, length(x)))
        w = cos.(w).^2
        return x.*w
    end



    df = DataFrame(CSV.File("$bulkDataSavePath/$caseName-rotors_convergence.csv"))

    df.T_1 = df.CT_1.*rho*(rotorRPM[1]/60)^2*(ROut[1]*2)^4
    df.Q_1 = df.CQ_1.*rho*(rotorRPM[1]/60)^2*(ROut[1]*2)^5




    plt1 = plot(df.T, df.CT_1, label="rotor 1")
    title!(plt1, "CT = T/(ρn^2d^4)")
    xlabel!("time[s]")

    plt2 = plot(df.T, df.CQ_1, label="rotor 1")
    title!(plt2, "CQ = Q/(ρn^2d^5)")
    xlabel!("time[s]")

    y = fft(hanning(df.CT_1)).*2
    y = real(y).^2+imag(y).^2
    y = y.^0.5
    x = fftfreq(length(df.CT_1), nsteps_per_rev)
    y = y[x.>0.5]
    x = x[x.>0.5]
    y = y[x.<12]
    x = x[x.<12]
    plt3 = plot(x, y, label="rotor 1")
    title!(plt3, "fft(CT)")
    xlabel!("per rev")


    y = fft(hanning(df.CQ_1)).*2
    y = real(y).^2+imag(y).^2
    y = y.^0.5
    x = fftfreq(length(df.CQ_1), nsteps_per_rev)
    y = y[x.>0.5]
    x = x[x.>0.5]
    y = y[x.<12]
    x = x[x.<12]
    plt4 = plot(x, y, label="rotor 1")
    title!(plt4, "fft(CQ)")
    xlabel!("per rev")



    plt5 = plot(df.T, df.T_1, label="rotor 1")
    title!(plt5, "Trust")
    xlabel!("time[s]")

    plt6 = plot(df.T, df.Q_1, label="rotor 1")
    title!(plt6, "Torque")
    xlabel!("time[s]")

    plt7 = plot(df.T, df.CT_1./df.CQ_1, label="rotor 1")
    title!(plt7, "Ct/CQ")
    xlabel!("time[s]")

    

    for ri in 1:numberOfRotors-1
        CT = df[:,5+ri*4]
        CQ = df[:,6+ri*4]

        T = CT.*rho*(rotorRPM[ri]/60)^2*(ROut[ri]*2)^4
        Q = CQ*rho*(rotorRPM[ri]/60)^2*(ROut[ri]*2)^5

        plot!(plt1, df.T, CT, label="rotor $(ri+1)")

        plot!(plt2, df.T, CQ, label="rotor $(ri+1)")

        y = fft(hanning(CT)).*2
        y = real(y).^2+imag(y).^2
        y = y.^0.5
        x = fftfreq(length(CT), nsteps_per_rev)
        y = y[x.>0.5]
        x = x[x.>0.5]
        y = y[x.<12]
        x = x[x.<12]
        plot!(plt3, x, y, label="rotor $(ri+1)")

        y = fft(hanning(CQ)).*2
        y = real(y).^2+imag(y).^2
        y = y.^0.5
        x = fftfreq(length(CQ), nsteps_per_rev)
        y = y[x.>0.5]
        x = x[x.>0.5]
        y = y[x.<12]
        x = x[x.<12]
        plot!(plt4, x, y, label="rotor $(ri+1)")

        plot!(plt5, df.T, T, label="rotor $(ri+1)")

        plot!(plt6, df.T, Q, label="rotor $(ri+1)")
        plot!(plt7,df.T, CT./CQ, label="rotor $(ri+1)")

    end


    nrows = 4
    ncols = 2
    plt = plot(plt1, plt2,plt3, plt4,plt5,plt6,plt7, layout=(nrows,ncols), size=(ncols*700.0, nrows*600.0), legend=false)

    savefig(plt, bulkDataSavePath*"Results.html")
    display(plt)
    

end ##if

####################################################################
if runParaview
    println("Calling Paraview...")

    # Files to open in Paraview
    files = joinpath(bulkDataSavePath, caseName*"_pfield...xmf;")
    for (ri, blade) in enumerate(numberOfBlades)
        for bi in 1:blade
            global files *= caseName*"_Rotor$(ri)_Blade$(bi)_loft...vtk;"
        end
    end
    files *= caseName*"_Wing_vlm...vtk;"

    # Call Paraview
    run(`paraview --data=$(files)`)

end
