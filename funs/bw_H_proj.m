 function [U,l] = bw_H_proj(U,x,y,V_0,sgn,l)

r = size(V_0,2);
p = norm((V_0'*U)*(U'*x))^2;
q = norm(U'*x)^2;

c_2 = q^2;
c_1 = -2*sgn*q;
c_0 = 1 - p/y;

desc = sqrt(c_1^2 - 4*c_2*c_0);     
l_(1) = (-c_1+real(desc))/(2*c_2);
l_(2) = (-c_1-real(desc))/(2*c_2);

if sgn == 1
    k = l_ < 1/q;
elseif sgn == -1
    k = l_ > -1/q;
end

l = min([l_(k),l]);
s = sqrt(abs(sgn*l/(1-sgn*l*q)));
sgn = sign(sgn*l/(1-sgn*l*q));

if sgn >= 0
    U = U*cholupdate(eye(r),s*U'*x,"+")';
else
    U = U*cholupdate(eye(r),s*U'*x,"-")';
end

end