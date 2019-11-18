# for GMM 
function SVmoments(m, n, θ, η, ϵ)
    S = size(ϵ, 2)
    ms = zeros(S,size(m,1))
    Threads.@threads for s=1:S
        ms[s,:] = sqrt(n)*aux_stat(SVmodel(θ, n, η[:,s], ϵ[:,s])[1])
    end
    ms .- m'
end
