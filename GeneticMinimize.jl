include("GeneticAlgorithmSubroutines.jl")
function geneticminimize(f,g,xlims,tol,n,kcross,kmutate,pmutate,psurvive,kspread)
    itters = 0
    generationssteady = 0
    fminhist = [1e30]
    fmaxhist = [1e30]
    population=newpopulation(xlims[1],xlims[2], n,kspread)
    populationhist = [copy(population)]
    converged = false
    while converged==false
        population,fmin,fmax = geneticstep(population,f,g,kcross,kmutate,pmutate,psurvive)
        push!(fminhist,fmin)
        push!(fmaxhist,fmax)
        push!(populationhist,population)
        if fmin == fminhist[end-1]
            generationssteady =generationssteady +1
        else
            generationssteady = 0
        end
        if generationssteady >tol
            converged = true
        end
        itters = itters + 1
    end
    return populationhist, fminhist, itters
end