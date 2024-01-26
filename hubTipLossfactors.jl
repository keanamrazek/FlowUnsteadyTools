using Plots
plotly()

Rtip = 2
Rhub = 0.5
B = 2

ri = LinRange(Rhub, Rtip, 1000)

hubtiploss_correction = ((1, 1, 1, 0.05), (2.6, 1, 1, 0.05))
function factor(r, theta)
    ((eh1, eh2, eh3, maxah), (et1, et2, et3, maxat)) = hubtiploss_correction

    fh = B/2*((r/Rhub)^eh1 - 1)^eh2/abs(sind(max(maxah, theta)))^eh3

    ft = B/2*((Rtip/r)^et1 - 1)^et2/abs(sind(max(maxat, theta)))^et3
    return 2/pi*acos(exp(-fh))*2/pi*acos(exp(-ft))

end

function prandtlFactor(r, theta)
    ((eh1, eh2, eh3, maxah), (et1, et2, et3, maxat)) = hubtiploss_correction

    ft = 1/2*((Rtip*B/r) - 1)/abs(sind(max(maxat, theta)))

    fh = 1/2*(B-(Rhub/r))/abs(sind(max(maxat, theta)))
    return 2/pi*acos(exp(-fh))*2/pi*acos(exp(-ft))

end

plt1 = plot(ri./Rtip, factor.(ri, 0), label="0 deg", legend=false)

for theta in 2:2:10
    plot!(ri./Rtip, factor.(ri, theta), label="$theta deg", legend=false)
end
title!(plt1, "lift loss factor")
xlabel!("r/R")
display(plt1)