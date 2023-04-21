clear all
addpath('../util')
addpath('../opt')
warning('off','MATLAB:sqrtm:SingularMatrix')

n = 2;
p = 2;
m = 2;
maxItrs = 1e3;

% hypercube constraints 
A = randn(m,n);
A = A./vecnorm(A,2,2);
A = [A;-A];
m = 2*m;

b_0 = 0.5*ones(m,1);
sgn = ones(m,1);

L_0 = randn(n,p);
Q_0 = L_0*L_0';
% Q_0 = eye(n);
sqrt_Q_0 = sqrtm(Q_0);

tic
cvx_begin sdp quiet
    variable Q(n,n)
    variable b(m,1) 
    minimize(trace(Q) + trace(Q_0) - 2*trace_sqrtm(sqrt_Q_0*Q*sqrt_Q_0))
    subject to 
        for i = 1:m
            A(i,:)*Q*A(i,:)' <= b_0(i);
        end
        Q >= 0;
        b >= 0;
cvx_end
toc

Q_sol = Q;
b_sol = b;


[U0,~,~] = svd(Q_0);
U0 = U0(:,1:p);
tic
[U_proj,~,U_i,l_i] = bw_proj(A,b_0,L_0,sgn,maxItrs,U0);
Q_proj = U_proj*U_proj'*Q_0*U_proj*U_proj';
toc

tic
[T_proj,~] = bwT_proj(A,b_0,Q_0,sgn,maxItrs,U0*U0');
Q_Tproj = T_proj*Q_0*T_proj;
toc

if n == 2
    
    angles = linspace(0,2*pi,100);
   
    figure(2)
    ellipse = sqrtm(Q_0)*[ cos(angles) ; sin(angles) ];
    plot(ellipse(1,:),ellipse(2,:),'k-');
    
    hold on

    ellipse = sqrtm(Q_sol)*[ cos(angles) ; sin(angles) ];
    plot(ellipse(1,:),ellipse(2,:),'b-');

    ellipse = sqrtm(Q_proj)*[ cos(angles) ; sin(angles) ];
    plot(ellipse(1,:),ellipse(2,:),'r-');
    
%     v = con2vert(A,sqrt(b_0));
%     [k,~] = convhull(v);
%     v = v(k,:);
%     plot(polyshape(v(:,1),v(:,2)),'FaceAlpha',0.5,'EdgeColor','k','FaceColor','b');
    
    hold off
 
end

warning('on','MATLAB:sqrtm:SingularMatrix')
