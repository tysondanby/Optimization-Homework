function eggshell(x)
    return 0.1*x[1]^2 +0.1*x[2]^2 -cos(3*x[1]) -cos(3*x[2])
end

function eggshell2d(a,b)
    x = [a,b]
    return 0.1*x[1]^2 +0.1*x[2]^2 -cos(3*x[1]) -cos(3*x[2])
end

function dummy(x)
    return zeros(length(x))
end