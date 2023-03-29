include("Genetics.jl")

# g should be less than or equal to 0
function geneticstep(population,f,g,kcross,kmutate,pmutate,psurvive)
    ninit =length(population)

    #EVALUATE f(x) AND g(x)
    populationF = zeros(ninit)
    feasible = ones(Bool,ninit)
    for i = 1:1:ninit
        populationF[i] = f(population[i])
        geval = g(population)
        for j = 1:1:length(geval)
            if geval[j] > 0.0
                feasible[i] = false
            end
        end
    end
    #Filter out the unfeasible
    feasiblepopulation,feasiblepopulationF=removeunfeasible(population,populationF,feasible)

    #ORDER FROM LOWEST TO HIGHEST
    order = sortperm(feasiblepopulationF)
    feasiblepopulation = feasiblepopulation[order]
    feasiblepopulationF = feasiblepopulationF[order]

    fmax = feasiblepopulationF[end]
    #Remove the least fit members of the population
    nkill = length(feasiblepopulation) - ninit*psurvive
    if nkill > 0 #Remove the highest value members
        for i = 1:1:nkill
            pop!(feasiblepopulation)
            pop!(feasiblepopulationF)
        end
    end

    #Repopulate and mutate
    repopulate!(feasiblepopulation, feasiblepopulationF,ninit,kcross)
    for i = 2:1:ninit #don't mutate the winner
        if rand() < pmutate
            feasiblepopulation[i] = mutate(feasiblepopulation[i],kmutate)
        end
    end
    #Return lowest function call and the new population
    return feasiblepopulation, feasiblepopulationF[1], fmax
end

function removeunfeasible(population,populationF,feasible)
    feasiblepopulation = []
    feasiblepopulationF = []
    for i = 1:1:length(population)
        if feasible[i] == true
            push!(feasiblepopulation,population[i])
            push!(feasiblepopulationF,populationF[i])
        end
    end
    return feasiblepopulation,feasiblepopulationF
end

#Population should come sorted from most to least fit.
function repopulate!(population,populationF,n,kcross)
    #only add children using push!()
    flow = populationF[1]
    fhigh = populationF[end]
    DF = 1.1*fhigh - 0.1*flow

    S = copy(populationF)
    for i in 1:1:length(population)
        fitness= (-populationF[i] + DF)/max(1,DF-flow)
        S[i] = fitness
    end
    S = S / sum(S)

    #print(length(population))
    #print(length(populationF))
    for i = 1:1:(n-length(population)) #number of offspring that need to be generated
        index1 = roulette(rand(),S)
        index2 = roulette(rand(),S)
        push!(population,crossover(population[index1],population[index2],kcross))
    end
    return newpopulation
end

function roulette(r,S)
    i = 0
    while r > 0.0
        i = i +1
        r=r-S[i]
    end
    return i
end