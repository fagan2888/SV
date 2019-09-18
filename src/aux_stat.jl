using Statistics
function aux_stat(y)
    y = abs.(y)
    end    
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
    vcat(clusters, HAR(y))
end
