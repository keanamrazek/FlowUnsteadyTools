#=
Description:
case script for rotor only systems

Author: Keana Mrazek

=#
caseName = "test1"

bulkDataSavePath = caseName

data_path        = "/home/theaero/Dropbox/VPMData/database" #/mnt/c/Dropbox/VPMData/database/"

savePath = false

runParaview=false

plotResults=true

runXFoil = true

save_wopwopin = true

################################################################
## Case description
fidelityLevel = "Low"  # Low, Mid, High

includeViscous = "CoreSpreading" # "Inviscid", "CoreSpreading", "ParticleStrengthExchange"

magVinf         = 1e-8    
magVvehicle     = 10.0 

AOA = 0.0

sideSlip = 0.0

temperature = 15.0 #deg C ground temperature

altitude = 1600.0 #m

nrevs    = 1 # Number of revolutions in simulation

domainLimits = 11.0^2 #simulation domain radius

gammaCrit = 0.000005^2  #minimum vorticity

gammaCritup = 10.0^2 #maximum vorticity

wakeTreatmentStep = 20 #wake treatment interval

domainCentre = [5.0, 0.0, 2.0]

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
RPMref = 688.0                                # Reference RPM to scale maneuver by
################################################################
#Rotor system description

ROut            = [2.0] #m

RInner          = [0.5] #m 

coningAngle     = [0.0] #m 

numberOfBlades  = [3]

ypos            = [0.0]

xpos            = [0.0]

zpos            = [0.0] 

roll            = [0.0] # of rotor
pitch           = [0.0]
yaw             = [0.0]


rotorAngle      = [90.0]

CWs             = [true]

bladeFile      = ["qx6_blade.csv"]       # Blade geometry

collective      = [0.0]

rotorRPM        = [688]

numberOfRotors = length(ROut)

referenceSpan = 1.5

aspectRatio = 5

ncrit  = 9 # Turbulence criterion for XFOIL

################################################################
# imports

import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm        

################################################################
# Pre proccessing
if fidelityLevel=="High"
    n_rotor         = 50                  # Number of blade elements per blade
    r_rotor         = 1/10                # Geometric expansion of elements
    nsteps_per_rev  = 360
    p_per_step      = 2                   # Sheds per time step
    vpm_integration = vpm.rungekutta3
elseif fidelityLevel=="Mid"
    n_rotor         = 20                  # Number of blade elements per blade
    r_rotor         = 1/10                # Geometric expansion of elements
    nsteps_per_rev  = 72
    p_per_step      = 4                   # Sheds per time step
    vpm_integration = vpm.rungekutta3
elseif fidelityLevel=="Low"
    n_rotor         = 10                  # Number of blade elements per blade
    r_rotor         = 1/10                # Geometric expansion of elements
    nsteps_per_rev  = 36
    p_per_step      = 4                   # Sheds per time step
    vpm_integration = vpm.euler
end

# time parameters nsteps/nsteps_per_rev / (RPM/60) 
nsteps = nrevs*nsteps_per_rev
ttot            = nsteps/nsteps_per_rev / (rotorRPM[1]/60)
max_particles   = Int(round(sum(numberOfBlades)*n_rotor*p_per_step*nsteps*2+1, digits=0))#nrotors*((2*n_rotor+1)*B)*nsteps*p_per_step + 1 # Maximum number of particles
# density altitude
TAlt = temperature - 0.00659*altitude
PAlt = 101.29*((TAlt+273.15)/288.08)^5.256
rho = PAlt/(0.2869*(TAlt+273.15))

#viscousity
mu_ = 1.85508e-5 
T_  = 15+273.15
S = 113
mu = mu_*(TAlt/T_)^(3/2)*(T_+S)/(TAlt+273.15+ S)

#speed of sound
c = sqrt(1.4*286.9*(TAlt+273.15))

#
ReD             = ((2*pi*rotorRPM./60).*ROut* rho/mu * 2).*ROut      # Diameter-based rotor Reynolds number
Matip           = (2*pi*rotorRPM./60/c).*ROut   # Tip Mach number

println("""
    Vref:   $(round(magVref, digits=1)) m/s
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
sigma_rotor_surf= ROut[1]/10                      # Rotor-on-VPM smoothing radius
lambda_vpm      = 2.125                     # VPM core overlap
                                            # VPM smoothing radius
sigma_vpm_overwrite = lambda_vpm * 2*pi*ROut[1]/(nsteps_per_rev*p_per_step) #can set to "nothing" then, lambda_vpm will be used
sigmafactor_vpmonvlm= 1                     # Shrink particles by this factor when
                                            #  calculating VPM-on-VLM/Rotor induced velocities
#Viscousity dissapation scheme
#vpm.Inviscid()`, `vpm.CoreSpreading(nu, sgm0, zeta)`, or `vpm.ParticleStrengthExchange(nu)
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
hubtiploss_correction = vlm.hubtiploss_nocorrection # Hub and tip correction
################################################################
# fountain effect
suppress_fountain   = false                  # Toggle

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
J = 0
for ri in 1:numberOfRotors

    rotor = uns.generate_rotor(ROut[ri], RInner[ri], numberOfBlades[ri], bladeFile[ri];
                                        pitch=collective[ri],
                                        n=n_rotor, CW=CWs[ri], blade_r=r_rotor,
                                        altReD=[rotorRPM[ri], J, mu/rho],
                                        xfoil=runXFoil,
                                        ncrit=ncrit,
                                        data_path=data_path,
                                        verbose=true,
                                        verbose_xfoil=true,
                                        plot_disc=true,
                                        ReD=findmax(ReD)[1],
                                        Matip=0,

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
Oaxis = uns.gt.rotation_matrix2(0, -AOA, 0)         # New orientation
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
    vlm_rlx         = 0.5,               # VLM relaxation (>0.9 can cause divergence, <0.2 slows simulation too much, deactivated with <0)
    #vlm_init        = false,            # Initialize the first step with the VLM semi-infinite wake solution
    hubtiploss_correction = hubtiploss_correction, # Hub and tip loss correction of rotors (ignored in quasi-steady solver)

    # Wake shedding
    #wake_coupled        = true,         # Couple VPM wake -> VLM solution
    shed_unsteady       = true,         # Whether to shed vorticity from unsteady loading
    unsteady_shedcrit   = unsteady_shedcrit,         # Criterion for unsteady-loading shedding
    shed_starting       = true,        # Whether to shed starting vortex (only when `shed_unsteady=true`)
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


####################################################################
if runParaview
    println("Calling Paraview...")

    # Files to open in Paraview
    files = joinpath(savePath, caseName*"_pfield...xmf;")
    for (ri, blade) in enumerate(numberOfBlades)
        for bi in 1:blade
            global files *= caseName*"_Rotor$(ri)_Blade$(bi)_loft...vtk;"
        end
    end
    files *= caseName*"_Wing_vlm...vtk;"

    # Call Paraview
    run(`paraview --data=$(files)`)

end