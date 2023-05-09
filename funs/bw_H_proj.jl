function bw_H_proj(U,x,y,V_0,sgn,l)

    r = size(V_0,2);
    p = norm((V_0'*U)*(U'*x))^2;
    q = norm(U'*x)^2;

    c_2 = q^2;
    c_1 = -2*sgn*q;
    c_0 = 1 - p/y;

    desc = sqrt(c_1^2 - 4*c_2*c_0);    
    l_ = [
        (-c_1+real(desc))/(2*c_2),
        (-c_1-real(desc))/(2*c_2)
        ];

    if sgn == 1
        k = l_ .< 1/q;
    elseif sgn == -1
        k = l_ .> -1/q;
    end

    l_next = minimum(push!(l_[k],l));
    s = sqrt(abs(sgn*l_next/(1-sgn*l_next*q)));
    sgn = sign(sgn*l_next/(1-sgn*l_next*q));
    U_next = U*cholUpdate(s*U'*x,sgn)';
    
    return  U_next, l_next
end