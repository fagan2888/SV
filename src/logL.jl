# asymptotic Gaussian likelihood function of statistic
function logL(θ, m, n, shocks_u, shocks_e, withdet=true)
    S = size(shocks_u,2)
    k = size(m,1)
    ms = zeros(eltype(SVmodel(θ, n, shocks_u[:,1], shocks_e[:,1])), S, k)
    # this loop could be parallelized!
    Threads.@threads for s = 1:S
        ms[s,:] = SVmodel(θ, n, shocks_u[:,s], shocks_e[:,s])
    end
    mbar = mean(ms,dims=1)[:]
    Σ = cov(ms)
    x = (m .- mbar)
    logL = try
        if withdet
            logL = -0.5*log(det(Σ)) - 0.5*x'*inv(Σ)*x # for Bayesian
        else    
            logL = 0.5*x'*inv(Σ)*x # for classic indirect inference (note sign change)
        end    
    catch
        logL = -Inf
    end
end
