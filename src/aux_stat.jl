
function aux_stat(y)
    y, m, s = stnorm(abs.(y))
    # look for evidence of volatility clusters
    mm = ma(y,5)
    mm = mm[5:end]
    clusters = 0.0
    try
        clusters = quantile(mm,0.75)-quantile(mm, 0.25)
    catch
        clusters = 1.0
    end    
    ϕ = HAR(y)
    vcat(m, s, clusters, ϕ)
end
