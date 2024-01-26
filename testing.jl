using CSV
using DataFrames

data_path        = "/home/theaero/Dropbox/VPMData/database"
bladeFile      = ["qx6_blade.csv"]

df = DataFrame(CSV.File(data_path*"/rotors/"*bladeFile[1]))

cordDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[1]))
pitchDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[2]))
sweepDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[3]))
heightDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[4]))
airfoilDist = DataFrame(CSV.File(data_path*"/rotors/"*df.file[5]))
coningAngle = 3.0
col = 5.0
M1 = [cordDist."r/R" cordDist."c/R"]
M2 = [pitchDist."r/R" pitchDist."twist (deg)".+col]
M3 = [sweepDist."r/R" sweepDist."y/R"]
M4 = [heightDist."r/R" heightDist."z/R"+sind(coningAngle)*heightDist."r/R"]

M5 = [(airfoilDist."r/R"[1], airfoilDist."Contour file"[1], airfoilDist."Aero file"[1])]

for i in 1:length(airfoilDist."r/R")-1
    push!(M5, (airfoilDist."r/R"[i+1], airfoilDist."Contour file"[i+1], airfoilDist."Aero file"[i+1]))
end




