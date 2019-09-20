using Statistics, Random
# version which generates shock internally
function SVmodel(θ, n, burnin)
    η = randn(n+burnin)
    ϵ = randn(n+burnin)
    SVmodel(θ, n, η, ϵ, false)
end    

# the dgp: simple discrete time stochastic volatility (SV) model
function SVmodel(θ, n, η, ϵ, savedata=false)
    α = 2.0*log(θ[1])
    ρ = θ[2]
    σ = θ[3]
    burnin = size(η,1) - n
    hlag = 0.0
    ys = zeros(n,1)
    for t = 1:burnin+n
        h = α + ρ*(hlag - α) + σ*η[t] # figure out type
        σt = exp(h/2.0)
        y = σt*ϵ[t]
        if t > burnin 
            ys[t-burnin] = y
        end    
        hlag = h
    end
    if savedata == true
        writedlm("svdata.txt", ys)
    end    
    ys
end
