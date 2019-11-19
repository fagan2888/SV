# asymptotic Gaussian likelihood function of statistic
function logL(θ, m, n, η, ϵ, withdet=true)
    S = size(η,2)
    k = size(m,1)
    ms = zeros(S, k)
    Threads.@threads for s = 1:S
        y, junk = SVmodel(θ, n, η[:,s], ϵ[:,s])
        y = min.(y, 100.0)
        y = max.(y,-100.0)
        ms[s,:] = sqrt(n)*aux_stat(y)
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
