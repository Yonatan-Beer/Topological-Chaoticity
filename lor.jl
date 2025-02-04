using Plots
using Ripserer
using PersistenceDiagrams
using LinearAlgebra: norm

function lorenz(σ, r, b, dt)
    function lor(v)
        xp = σ.*(v[2] - v[1])
        yp = v[1].*(r - v[3]) - v[2]
        zp = v[1].*v[2] - b.*v[3]
        return(v .+ dt.*(xp,yp,zp))
    end
    return lor
end

function process(f,init,iters)
    cur = init
    out = Vector{Tuple{Float64,Float64,Float64}}(undef,iters)
    for i in 1:iters
        cur = f(cur)
        out[i] = cur
    end
    return out
end

function mag(r, dim=2)
    f = lorenz(10.0, r, 8/3, 0.001)
    x0 = (1.0,1.0,1.0)
    
    noise = 200
    iters = 100000+noise

    sample = process(f, x0, iters)[noise+1:20:end]
    diagram = ripserer(Alpha, sample; dim_max=2)
    form = PersistenceImage(diagram; size = 120)
    image = form(diagram[dim])
    return norm(image)
end

dt = 0.001
f = lorenz(10.0, 28.0, 8/3, dt)
iters = 100000+200
noise = 200

x0 = (1.0,1.0,1.0)
sample = process(f, x0, iters)[noise+1:1:end]

#x = [p[1] for p in sample]
#y = [p[2] for p in sample]
#z = [p[3] for p in sample]
#colors = range(0, 1, length=length(x)) 
#plot(x,y,z, linez= colors, c = :thermal, legend=false)

sample = sample[begin:20:end]
diagram = ripserer(Alpha, sample; dim_max=2)

form = PersistenceImage(diagram; size = 120)
image = form(diagram[2])
#heatmap(1:120, 1:120, image)

