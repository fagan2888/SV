using Statistics
function aux_stat(y)
    α = sqrt(mean(y.^2.0))
    y = abs.(y)
    # look for evidence of volatility clusters, for ρ
    mm = ma(y,5)
    mm = mm[5:end]
    clusters = 0.0
    try
        clusters = quantile(mm,0.75)-quantile(mm, 0.25)
    catch
        clusters = 1.0
    end
    # HAR model, for all params
    vcat(α, clusters, HAR(y))
end
