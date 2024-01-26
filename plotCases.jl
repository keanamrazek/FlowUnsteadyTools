using Plots
using Plots.PlotMeasures
using CSV
using DataFrames
using FFTW

plotly()
function average(x)
    x = x./length(x)
    return sum(x)

end

bulkDataSavePath = "/home/aero/Desktop/tempVPMFilesKeana/jack"
n = 2500 #rpm
d = 0.420*2 # rotorDiameter
rho = 0.984096785411496 #air density

V = [0, 10, ]
Pitch= [18,20,22,24]#vec(collect(8:4:24))

CT = zeros((length(V), length(Pitch)))
CP = zeros((length(V), length(Pitch)))
CQ = zeros((length(V), length(Pitch)))
T = zeros((length(V), length(Pitch)))
P = zeros((length(V), length(Pitch)))
Q = zeros((length(V), length(Pitch)))
eta = zeros((length(V), length(Pitch)))
J = zeros((length(V), length(Pitch)))

for i in 1:length(V)
    for j in 1:length(Pitch)
        v=V[i]
        p=Pitch[j]
        df = DataFrame(CSV.File("/home/aero/Desktop/tempVPMFilesKeana/jack/V$(v)/VPMCaseP$(Int(p))/VPMCaseP$(Int(p))-rotors_convergence.csv"))
        df.T_1 = df.CT_1.*rho*(n/60)^2*d^4
        df.Q_1 = df.CQ_1.*rho*(n/60)^2*d^5
        df.P_1 = df.Q_1.*n*2pi/60
        df.CP_1 = df.P_1./(rho*(n/60)^3*d^5)
        J[i, j]  = v/(n*2pi/60)/d
        CT[i, j] = average(df.CT_1[end-90:end])
        CP[i, j] = average(df.CP_1[end-90:end])
        CQ[i, j] = average(df.CQ_1[end-90:end])
        T[i, j] = average(df.T_1[end-90:end])
        P[i, j] = average(df.P_1[end-90:end])
        Q[i, j] = average(df.Q_1[end-90:end])
        eta[i, j] = average(df.eta_1[end-90:end])
    end
end



plt1 = plot(legend=true)
ylabel!(plt1, "CT")
xlabel!(plt1, "Pitch [deg]")

plt2 = plot(legend=false)
ylabel!(plt2, "CQ")
xlabel!(plt2, "Pitch [deg]")

plt3 = plot(legend=false)
ylabel!(plt3, "CP")
xlabel!(plt3, "Pitch [deg]")

plt4 = plot(legend=false)
ylabel!(plt4, "CT/CP")
xlabel!(plt4, "Pitch [deg]")

plt5 = plot(legend=false)
ylabel!(plt5, "eta")
xlabel!(plt5, "Pitch [deg]")

plt6 = plot(legend=false)
ylabel!(plt6, "CT/CQ")
xlabel!(plt6, "Pitch [deg]")


for i in 1:length(V)

    j=J[i,1]
    plot!(plt1, Pitch, CT[i,:], label="AR: $(round(j*1e3)/1e3)")
    plot!(plt2, Pitch, CQ[i,:], label=nothing)
    plot!(plt3, Pitch, CP[i,:]./1000, label=nothing)
    plot!(plt4, Pitch, CT[i,:]./CP[i,:], label=nothing)
    plot!(plt5, Pitch, eta[i,:], label=nothing)
    plot!(plt6, Pitch, CT[i,:]./CQ[i,:], label=nothing)
end

ylims!(plt4, (2,4))
yticks!(plt4, 2:0.2:4)

ylims!(plt5, (0,2))
yticks!(plt5, 0:(4/20):2)

ylims!(plt6, (12,24))
yticks!(plt6, 12:0.5:24)

nrows = 3
ncols = 2
plt = plot(plt1, plt2,plt3,plt4, plt5,plt6,  layout=(nrows,ncols), size=(ncols*400.0, nrows*300.0), legend=true, left_margin = 30px, bottom_margin = 50px)
savefig(plt, bulkDataSavePath*"ResultsSummary.html")
display(plt)


plt1 = plot(legend=true)
ylabel!(plt1, "Thrust [N]")
xlabel!(plt1, "Velocity [m/s]")

plt2 = plot(legend=false)
ylabel!(plt2, "Torque [Nm]")
xlabel!(plt2, "Velocity [m/s]")

plt3 = plot(legend=false)
ylabel!(plt3, "Power [kW]")
xlabel!(plt3, "Velocity [m/s]")

plt4 = plot(legend=false)
ylabel!(plt4, "CT/CP")
xlabel!(plt4, "Velocity [m/s]")

plt5 = plot(legend=false)
ylabel!(plt5, "eta")
xlabel!(plt5, "Velocity [m/s]")

plt6 = plot(legend=false)
ylabel!(plt6, "CT/CQ")
xlabel!(plt6, "Velocity [m/s]")



for i in 1:length(Pitch)

    p=Pitch[i]
    plot!(plt1, V, T[:,i], label="$p deg pitch")
    plot!(plt2, V, Q[:,i],label=nothing)
    plot!(plt3, V, P[:,i]./1000,label=nothing)
    plot!(plt4, V, CT[:,i]./CP[:,i],label=nothing)
    plot!(plt5, V, eta[:,i],label=nothing)
    plot!(plt6, V, CT[:,i]./CQ[:,i],label=nothing)
end

ylims!(plt4, (2,4))
yticks!(plt4, 2:0.2:4)

ylims!(plt5, (0,2))
yticks!(plt5, 0:(4/20):2)
ylims!(plt6, (12,24))
yticks!(plt6, 12:0.5:24)


nrows = 6
ncols = 2
plt = plot(plt1, plt2,plt3,plt4,plt5,plt6,  layout=(nrows,ncols), size=(ncols*400.0, nrows*300.0), legend=true, left_margin = 30px, bottom_margin = 50px)
savefig(plt, bulkDataSavePath*"ResultsSummary2.html")
display(plt)



plt1 = plot(legend=true)
ylabel!(plt1, "CT")
xlabel!(plt1, "Advance Ratio")

plt2 = plot(legend=false)
ylabel!(plt2, "CQ")
xlabel!(plt2, "Advance Ratio")

plt3 = plot(legend=false)
ylabel!(plt3, "CP")
xlabel!(plt3, "Advance Ratio")

plt4 = plot(legend=false)
ylabel!(plt4, "CT/CP")
xlabel!(plt4, "Advance Ratio")

plt5 = plot(legend=false)
ylabel!(plt5, "eta")
xlabel!(plt5, "Advance Ratio")

plt6 = plot(legend=false)
ylabel!(plt6, "CT/CQ")
xlabel!(plt6, "Advance Ratio")


for i in 1:length(Pitch)

    p=Pitch[i]
    plot!(plt1, J[:,i], CT[:,i], label="$p deg pitch")
    plot!(plt2, J[:,i], CQ[:,i],label=nothing)
    plot!(plt3, J[:,i], CP[:,i]./1000,label=nothing)
    plot!(plt4, J[:,i], CT[:,i]./CP[:,i],label=nothing)
    plot!(plt5, J[:,i], eta[:,i],label=nothing)
    plot!(plt6, J[:,i], CT[:,i]./CQ[:,i],label=nothing)
end

ylims!(plt4, (2,4))
yticks!(plt4, 2:0.2:4)

ylims!(plt5, (0,2))
yticks!(plt5, 0:(4/20):2)

ylims!(plt6, (12,24))
yticks!(plt6, 12:0.5:24)

nrows = 6
ncols = 2
plt = plot(plt1, plt2,plt3,plt4,plt5,plt6,  layout=(nrows,ncols), size=(ncols*400.0, nrows*300.0), legend=true, left_margin = 30px, bottom_margin = 50px)
savefig(plt, bulkDataSavePath*"ResultsSummary3.html")
display(plt)
