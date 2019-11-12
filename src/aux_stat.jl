using Statistics

function aux_stat(y)
    α = sqrt(mean(y.^2.0))
    y = abs.(y)
    m = mean(y)
    s = std(y)
    k = std(y.^2.0)
    # look for evidence of volatility clusters, for ρ
    mm = ma(y,5)
    mm = mm[5:end]
    clusters1 = quantile(mm,0.75)-quantile(mm, 0.25)
    clusters2 = quantile(mm,0.9)-quantile(mm, 0.1)
    mm = ma(y,10)
    mm = mm[10:end]
    clusters3 = quantile(mm,0.75)-quantile(mm, 0.25)
    clusters4 = quantile(mm,0.9)-quantile(mm, 0.1)
    # HAR model, for all params
    vcat(α, m, s, k, clusters1, clusters2, clusters3, clusters4, HAR(y))
end


#=
function RawMoment(θ, m)
    α = 2.0*log(θ[1])
    ρ = θ[2]
    σ = θ[3]
    exp(m*α/2.0 + (m^2.0)*σ/(1.0 - ρ^2.0)/4.0) # Shepphard GMM notes, page 384
end

function UnconditionalMoments(θ)
    α = 2.0*log(θ[1])
    ρ = θ[2]
    σ = θ[3]
    vcat(
        sqrt(2.0/pi)*RawMoment(θ,1),   # mean (abs(y))
        RawMoment(θ,2) # mean y^2
        )
        #2.0*sqrt(2.0/pi)*RawMoment(θ,3), # mean abs(y^3))
        #3.0*RawMoment(θ,4)               # mean y^4
        #)
end
=##=
function aux_stat(y, θ)
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
    vcat(α, clusters, HAR(y), UnconditionalMoments(θ))
end
=#
