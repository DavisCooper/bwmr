function bw_proj(X, Y, V_0, sgn, iters, U)

    m,_ = size(X);
    l = zeros(m,1);
    # l_ = zeros(m,1);
    N = ceil(iters/m); 
    
    U_i = Array{Matrix}(undef,0);
    count = 1;
    
    for i=1:N
        for k=1:m 
            U,a = bw_H_proj(U,X[k,:],Y[k],V_0,sgn[k],copy(l[k]));
            l[k] = l[k] - a;
            push!(U_i,U);
            count = count + 1;
        end
        
    # norm_sum = norm(l) + norm(l_);
    # if norm_sum == 0
    #     break
    # else
    #     conv = norm(l_ - l,1)/norm_sum;
    #     if conv < tol
    #         break
    #     end
    #     l_ = l;        
    # end
        
    end
        return U,l,U_i
end