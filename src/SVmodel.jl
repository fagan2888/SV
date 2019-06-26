# version which generates shock internally
function SVmodel(θ, n, burnin)
    shocks_u = randn(n+burnin,1)
    shocks_e = randn(n+burnin,1)
    SVmodel(θ, n, shocks_u, shocks_e, savedata=false)
end    

# the dgp: simple discrete time stochastic volatility (SV) model
function SVmodel(θ, n, shocks_u, shocks_e, savedata=false)
    σe = θ[1]
    ρ = θ[2]
    σu = θ[3]
    burnin = size(shocks_u,1) - n
    hlag = 0.0
    h = ρ.*hlag .+ σu.*shocks_u[1] # figure out type
    y = σe.*exp(h./2.0).*shocks_e[1]
    ys = zeros(n,1)
    for t = 1:burnin+n
        h = ρ.*hlag .+ σu.*shocks_u[t]
        y = σe.*exp(h./2.0).*shocks_e[t]
        if t > burnin 
            ys[t-burnin] = y
        end    
        hlag = h
    end
    if savedata == true
        writedlm("svdata.txt", ys)
    end    
    sqrt(n)*aux_stat(ys)
end
