using CSV
using DataFrames
using Printf


polars = ["noble1.datRE1061000Mach62.csv", 
"noble2.datRE1184000Mach61.csv", 
"noble3.datRE1264000Mach59.csv", 
"noble4.datRE1318000Mach57.csv", 
"noble5.datRE1355000Mach55.csv", 
"noble6.datRE1379000Mach54.csv", 
"noble7.datRE1392000Mach52.csv", 
"noble8.datRE1398000Mach50.csv", 
"noble9.datRE1392000Mach47.csv", 
"noble10.datRE1367000Mach43.csv", 
"noble11.datRE1318000Mach39.csv", 
"noble12.datRE1197000Mach34.csv", 
"noble13.datRE1019000Mach29.csv", 
"noble14.datRE795000Mach24.csv", 
"noble15.datRE521000Mach18.csv"]

for j in 1:length(polars)
    polar = DataFrame(CSV.File("polars/"*polars[j]))

    text = """#
    ##=---------------------------------------------------------------------------*| 
    #|            #####        |                                                   | 
    #|       ####              |                                                   | 
    #|    ####    ###   ###    | Bronberg Dynamics                                 | 
    #|  #####     #  #  #  #   | Web: https://bronbergdynamics.co.za               | 
    #| #####      ###   #  #   | Author: Keana Mrazek                              | 
    #|  ####      #  #  #  #   | Email: Keana@bronbergdynamics.co.za               | 
    #|    ###     ###   ###    |                                                   | 
    #|      ##                 |                                                   | 
    #|*---------------------------------------------------------------------------=#
    !T AOACD
    !U deg\n"""

    for i in polar.Alpha
        text = text*@sprintf("%.6f\n",i)
    end
    text =text*"!T AOACL\n!U deg\n"
    for i in polar.Alpha
        text = text*@sprintf("%.6f\n",i)
    end
    text =text*"!T AOACM\n!U deg\n"
    for i in polar.Alpha
        text = text*@sprintf("%.6f\n",i)
    end
    text =text*"!T MACHCD\n!U deg\n"
    text =text*"0."*polars[j][end-5:end-4]*"\n"
    text =text*"!T MACHCL\n!U deg\n"
    text =text*"0."*polars[j][end-5:end-4]*"\n"
    text =text*"!T MACHCM\n!U deg\n"
    text =text*"0."*polars[j][end-5:end-4]*"\n"

    text =text*"!M CDTAB\n!U nd\n"
    for i in polar.Cd
        text = text*@sprintf("%.6f\n",i)
    end

    text =text*"!M CLTAB\n!U nd\n"
    for i in polar.Cl
        text = text*@sprintf("%.6f\n",i)
    end

    text =text*"!M CMTAB\n!U nd\n"
    for i in polar.Cm
        text = text*@sprintf("%.6f\n",i)
    end

    file = open("flightLab/$(polars[j][1:end-4]).tab", "w")
    write(file, text)
    close(file)
end