include("Eggshell.jl")
using Plots; pythonplot()
function plothist(pophist,f,xlims)
    frames = []
    anim = Animation()
    for i = 1:1:length(pophist)
        x = range(xlims[1][1], xlims[2][1], length=100)
        y = range(xlims[1][2], xlims[2][2], length=100)
        z = @. eggshell2d(x',y)
        contour(x, y, z)
        x1 =[]
        x2 =[]
        for x in pophist[i]
            push!(x1,x[1])
            push!(x2,x[2])
        end
        scatter!(x1,x2)
        frame(anim)
    end
    gif(anim, "anim_fps2.gif", fps = 2)
end