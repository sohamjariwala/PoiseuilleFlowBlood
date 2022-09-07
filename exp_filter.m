function sigma = exp_filter(m,order)

sigma = exp(log(eps)*([0:(m-1)]/(m-1)).^order);

end