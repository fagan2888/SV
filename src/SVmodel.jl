using Statistics, Random

# version which generates shock internally
function SVmodel(θ, n, burnin)
    ϕ = θ[1]
    ρ = θ[2]
    σ = θ[3]
    hlag = 0.0
    ys = zeros(n)
    for t = 1:burnin+n
        h = ρ*hlag + σ*randn()
        y = ϕ*exp(h/2.0)*randn()
        if t > burnin 
            ys[t-burnin] = y
        end    
        hlag = h
    end
    ys, ϕ*exp(hlag/2.0) # return the sample of returns, plus the final period volatility
end
# the dgp: simple discrete time stochastic volatility (SV) model
function SVmodel(θ, n, η, ϵ, savedata=false)
    ϕ = θ[1]
    ρ = θ[2]
    σ = θ[3]
    burnin = size(η,1) - n
    hlag = 0.0
    ys = zeros(n,1)
    for t = 1:burnin+n
        h = ρ*hlag + σ*η[t]
        σt = ϕ*exp(h/2.0)
        y = σt*ϵ[t]
        if t > burnin 
            ys[t-burnin] = y
        end    
        hlag = h
    end
    if savedata == true
        writedlm("svdata.txt", ys)
    end    
    ys, ϕ*exp(hlag/2.0) # return the sample of returns, plus the final period volatility
end
