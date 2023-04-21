 function [T,l] = bwT_H_proj(T,a,Q_0,b_0,sgn,l)

p = a'*T*Q_0*T*a;
q = a'*T*a;

c_2 = q^2;
c_1 = -2*sgn*q;
c_0 = 1 - p/b_0;

desc = sqrt(c_1^2 - 4*c_2*c_0);     
l_(1) = (-c_1+real(desc))/(2*c_2);
l_(2) = (-c_1-real(desc))/(2*c_2);

if sgn == 1
    k = l_ < 1/q;
elseif sgn == -1
    k = l_ > -1/q;
end

l = min([l_(k),l]);
T = T + sgn*l*T*a*a'*T/(1-sgn*l*q);

 end