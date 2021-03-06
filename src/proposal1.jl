# uniform random walk
function proposal1(current, tuning, lb, ub)
    trial = copy(current)
    if rand() > 0.1
        i = rand(1:size(trial,1))
        trial[i] += tuning[i].*randn()
    else
        trial = lb + (ub - lb).*rand(size(lb,1))
    end    
    return trial
end
