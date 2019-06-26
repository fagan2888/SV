# for GMM 
function SVmoments(m, n, θ, shocks_u, shocks_e)
    S = size(shocks_u, 2)
    ms = zeros(S,size(m,1))
    Threads.@threads for s=1:S
        ms[s,:] = SVmodel(θ, n, shocks_u[:,s], shocks_e[:,s])
    end
    ms .- m'
end
