using BenchmarkTools
using LinearAlgebra
using Plots
include("bw_proj.jl")
include("funs/bw_H_proj.jl")
include("funs/cholUpdate.jl")

n = 2;
p = 2;
m = 2;
maxItrs = 1e2;

# hypercube constraints 
A = randn(m,n);
A = A./[norm(A[i,:]) for i=1:size(A,1)];
A = [A;-A];
m = 2*m;

b_0 = 0.5*ones(m,1);
sgn = ones(m,1);

L_0 = randn(n,p);
Q_0 = L_0*L_0'; 

U0,_,_ = svd(Q_0);

U_proj,l,U_i = bw_proj(A,b_0,L_0,sgn,maxItrs,U0);
Q_proj =  U_proj*U_proj'*Q_0*U_proj*U_proj';

d = Array{Real}(undef,Int(maxItrs));
for i=1:Int(maxItrs)
    T_i = U_i[i]*U_i[i]';
    d[i] = max(real(tr(Q_proj + T_i*Q_0*T_i - 2*sqrt(Q_proj*T_i*Q_0*T_i))),1e-15);
end

plot([1:maxItrs],d)
plot!(yscale=:log10)