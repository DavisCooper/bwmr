function [T,l,varargout] = bwT_proj(X, Y, Q_0, sgn, iters, T)

[m,n] = size(X);
l = zeros(m,1);
l_ = zeros(m,1);
N = ceil(iters/m); 

T_i = {};
count = 1;

for i=1:N
    for k=1:m 
        [T,a] = bwT_H_proj(T,X(k,:)',Q_0,Y(k),sgn(k),l(k));
        l(k) = l(k) - a;
        T_i{count} = T;
        count = count + 1;
%         Q_i{(i-1)*m+k} = sqrt_Q_0\Q*Q/sqrt_Q_0;
    end
    
%     norm_sum = norm(l) + norm(l_);
%     if norm_sum == 0
%         break
%     else
%         conv = norm(l_ - l,1)/norm_sum;
%         if conv < tol
%             break
%         end
%         l_ = l;
%         
%     end
    
end

varargout{1} = T_i;
end