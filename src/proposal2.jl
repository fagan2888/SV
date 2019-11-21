function proposal2(trial, cholV, lb, ub)
    if rand() > 0.1
        trial += cholV'*randn(size(trial))
    else
        trial = lb + (ub - lb).*rand(size(lb,1))
    end
    return trial
end
