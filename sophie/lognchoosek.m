function l = lognchoosek(n,k,flag)
if strcmp(flag,'int')
    l = sum(log(1:n)) - sum(log(1:(n-k))) - sum(log(1:k));
end
if strcmp(flag,'real')
    l = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
end