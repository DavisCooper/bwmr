function [U,l,varargout] = bw_proj(X, Y, V_0, sgn, iters, U)

[n,d] = size(X);
l = zeros(n,1);
l_ = zeros(n,1);
N = ceil(iters/n); 

U_i = {};
l_i = {};

count = 1;
for i=1:N
    for k=1:n 
        [U,a] = bw_H_proj(U,X(k,:)',Y(k),V_0,sgn(k),l(k));
        l(k) = l(k) - a;
        U_i{count} = U;
        count = count + 1;
    end
    l_i{i} = l;

% *** not using convergence criteria ***
%     norm_sum = norm(l) + norm(l_);
%     if norm_sum == 0
%         break
%     else
%         conv = norm(l_ - l,1)/norm_sum;
%         if conv < tol
%             break
%         end
%         l_ = l;
%     end
    
end

varargout{1} = U_i;
varargout{2} = l_i;

end