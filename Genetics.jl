
#k = (0,inf) higher k means child is more 
#likely a direct average of the parent
#Lower k means that the children inherit traits 
#more similar to one parent or the other.
# k = 1 gives a flat distribution
function crossover(x1,x2,k)
    crossoverparam1 = rand(length(x1))
    @. crossoverparam1 = distribute(crossoverparam1,k)
    crossoverparam2 = ones(length(x1)) - crossoverparam1
    return crossoverparam1.*x1 + crossoverparam2.*x2
end

# k is the maximum magnitude of a mutation
# as a fraction of x1
function mutate(x1,k)
    return x1 + k*((2*rand(length(x1))-ones(length(x1))).*x1)
end

function newpopulation(xmin,xmax, n,kspread)
    population = []
    for i = 1:1:n
        push!(population,crossover(xmin,xmax,kspread))
    end
    return population
end

function distribute(p,k)
    p = 2*(p -0.5)
    p = sign(p)*abs(p)^k
    return p/2 + 0.5
end