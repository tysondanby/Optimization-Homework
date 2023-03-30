include("GeneticMinimize.jl")
include("Eggshell.jl")
include("Visualization.jl")
using LinearAlgebra

numruns = 100
notminimum = 0.1
#Eggshell - initial test
xmin = [-5.0,-5.0]
xmax = [5.0,5.0]
xlims = [xmin,xmax]
tol = 3
n = 50
kcross = 2.0
kmutate = 1.0
pmutate = 0.03
psurvive = 0.33
kspread = 1.0

pophist, fminhist, itters =geneticminimize(eggshell,dummy,xlims,tol,n,kcross,kmutate,pmutate,psurvive,kspread)

println(fminhist[end])
#println(eggshell(pophist[end][1]))
println(pophist[end][1])
println("Itterations: $itters")
#plothist(pophist,eggshell,xlims)


gr()
#Sensitivity study 1
ns = 20:1:250
itterations = []
errors = []
global itters
Threads.@threads for i = 1:1:length(ns)
    global itters = 0.0
    minimums = 0.0
    for j = 1:1:numruns
        @label one
        local fminhist
        local pophist
        pophist, fminhist, nitters =geneticminimize(eggshell,dummy,xlims,tol,ns[i],kcross,kmutate,pmutate,psurvive,kspread)
        if norm(pophist[end][1]) > notminimum
            @goto one
        end
        itters = itters + nitters
        minimums = minimums + norm(pophist[end][1])
    end
    push!(itterations,itters/numruns)
    push!(errors,minimums/numruns)
end

p1 = plot(ns,itterations,xlabel = "Population Size",ylabel = "Itterations to Converge")
p2 = plot(ns,errors,xlabel = "Population Size",ylabel = "Absolute Error")


#Sensitivity study 2
psurvives = 0.05:0.01:0.75
itterations2 = []
errors2 = []
Threads.@threads for i = 1:1:length(psurvives)
    global itters = 0.0
    minimums = 0.0
    for j = 1:1:numruns
        @label two
        local fminhist
        local pophist
        pophist, fminhist, nitters =geneticminimize(eggshell,dummy,xlims,tol,n,kcross,kmutate,pmutate,psurvives[i],kspread)
        if norm(pophist[end][1]) > notminimum
            @goto two
        end
        itters = itters + nitters
        minimums = minimums + norm(pophist[end][1])
    end
    push!(itterations2,itters/numruns)
    push!(errors2,minimums/numruns)
end
#println(length(itterations2))
p3 = plot(psurvives,itterations2,xlabel = "Fraction of Population Retained in New Generation",ylabel = "Itterations to Converge")
p4 = plot(psurvives,errors2,xlabel = "Fraction of Population Retained in New Generation",ylabel = "Absolute Error")


#Sensitivity study 3
pmutates = 0.005:0.005:0.25
itterations3 = []
errors3 = []
Threads.@threads for i = 1:1:length(pmutates)
    global itters = 0.0
    minimums = 0.0
    for j = 1:1:numruns
        @label three
        local fminhist
        local pophist
        pophist, fminhist, nitters =geneticminimize(eggshell,dummy,xlims,tol,n,kcross,kmutate,pmutates[i],psurvive,kspread)
        if norm(pophist[end][1]) > notminimum
            @goto three
        end
        itters = itters + nitters
        minimums = minimums + norm(pophist[end][1])
    end
    push!(itterations3,itters/numruns)
    push!(errors3,minimums/numruns)
end

p5 = plot(pmutates,itterations3,xlabel = "Mutation Probability",ylabel = "Itterations to Converge")
p6 = plot(pmutates,errors3,xlabel = "Mutation Probability",ylabel = "Absolute Error")

#
n = 41
lastconvergence = 1.0
while lastconvergence >= 0.5
    global n = n - 1
    fail = 0
    #println("Line 114")
    Threads.@threads for i = 1:1:numruns
        history, ~, ~ =geneticminimize(eggshell,dummy,xlims,tol,n,kcross,kmutate,pmutate,psurvive,kspread)
        #println(history)
        if norm(history[end][1]) > notminimum
            fail = fail +1
        end
    end
    global lastconvergence = (numruns-fail)/numruns
end
println("Size to Fail 50% of time: $n")

#5.2
n = 50
xlims = [[-1.0,-1.0,-1.0],[2.0,2.0,2.0]]
pophist, fminhist, itters =geneticminimize(hartman,dummy,xlims,tol,n,kcross,kmutate,pmutate,psurvive,kspread)
println("Problem 5.2")
println(fminhist[end])
#println(hartman(pophist[end][1]))
println(pophist[end][1])
println("Itterations: $itters")

using Optim
println()
result = Optim.optimize(hartman,xlims[1],ParticleSwarm(; lower = xlims[1], upper = xlims[2], n_particles = n))