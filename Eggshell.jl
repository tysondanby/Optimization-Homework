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

function hartman(x)
    alpha = [1.0,1.2,3.0,3.2]
    A = [3.0 10.0 30.0;0.1 10.0 35.0;3.0 10.0 30.0;0.1 10.0 35.0]
    P = 1e-4*[3689.0 1170.0 2673.0;4699.0 4387.0 7470.0;1091.0 8732.0 5547.0;381.0 5743.0 8828.0]
    sum = 0
    for i = 1:1:4
        smallsum = 0
        for j = 1:1:3
            smallsum = smallsum - A[i,j]*(x[j]-P[i,j])^2
        end
        sum = sum - alpha[i]*exp(smallsum)
    end
    return sum
end