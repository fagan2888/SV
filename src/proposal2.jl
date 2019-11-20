function proposal2(current, cholV, lb, ub)
    trial = copy(current)
    if rand() > 0.1
        trial += cholV'*randn(size(trial))
    else
        trial = lb + (ub - lb).*rand(size(lb))
    end
    return trial
end
