function cholUpdate(x,sgn)
    
d = size(x,1);
L = Matrix{Real}(I,d,d);

for k = 1:d
    r = sqrt(L[k, k]^2 + sgn*x[k]^2);
    c = r / L[k, k];
    s = x[k] / L[k, k];
    L[k, k] = r;
    if k < d
        L[(k+1):d, k] = (L[(k+1):d, k] + sgn * s * x[(k+1):d]) / c;
        x[(k+1):d] = c * x[(k+1):d] - s * L[(k+1):d, k];
    end
end
       
return L

end