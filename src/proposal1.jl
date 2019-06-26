# uniform random walk, with bounds check
function proposal1(current, tuning, lb, ub)
    trial = copy(current)
    if rand() > 0.1
        i = rand(1:size(current,1))
        trial[i] = current[i] + tuning[i].*randn()
    else
        trial = lb + (ub - lb).*rand(size(lb))
    end    
    return trial
end
