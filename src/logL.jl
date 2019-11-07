# asymptotic Gaussian likelihood function of statistic
function logL(θ, m, n, η, ϵ, withdet=true)
    S = size(η,2)
    k = size(m,1)
    ms = zeros(S, k)
    # this loop could be parallelized!
    Threads.@threads for s = 1:S
        ms[s,:] = sqrt(n)*aux_stat(SVmodel(θ, n, η[:,s], ϵ[:,s])[1]) # SVmodel returns last period vol, don't use it
    end
    mbar = mean(ms,dims=1)[:]
    if ~any(isnan.(mbar))
        Σ = cov(ms)
        x = (m .- mbar)
        lnL = try
            if withdet
                lnL = -0.5*log(det(Σ)) - 0.5*x'*inv(Σ)*x # for Bayesian
            else    
                lnL = 0.5*x'*inv(Σ)*x # for classic indirect inference (note sign change)
            end    
        catch
            lnL = -Inf
        end
     else
         lnL = -Inf
     end
     return lnL
end
