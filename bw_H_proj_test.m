clear all
addpath('funs')
warning('off','MATLAB:sqrtm:SingularMatrix')

d = 2;
p = 2;
m = 1;
maxItrs = 1e3;

% hypercube constraints 
x = randn(d,1);
x = x/norm(x);
y = 0.5;
sgn = 1;

L_0 = randn(d,p);
Q_0 = L_0*L_0';

[U0,~,~] = svd(Q_0);
U0 = U0(:,1:p);

tic
[U_proj,l_proj] = bw_H_proj(U0,x,y,L_0,sgn,0);
Q_proj = U_proj*U_proj'*L_0*L_0'*U_proj*U_proj';
toc

tic
[T_proj,l_Tproj] = bwT_H_proj(U0*U0',x,Q_0,y,sgn,0);
Q_Tproj = T_proj*Q_0*T_proj;
toc

if d == 2
    
    x_max = -inf;
    x_min = inf;
    
    angles = linspace(0,2*pi,100);
   
    ellipse = real(sqrtm(Q_0))*[ cos(angles) ; sin(angles) ];
    x_min = min([ellipse(1,:),x_min]);
    x_max = max([ellipse(1,:),x_max]);
    plot(ellipse(1,:),ellipse(2,:),'k-');
    
    hold on

    ellipse = real(sqrtm(Q_Tproj))*[ cos(angles) ; sin(angles) ];
    x_min = min([ellipse(1,:),x_min]);
    x_max = max([ellipse(1,:),x_max]);
    plot(ellipse(1,:),ellipse(2,:),'g-');
    
    ellipse = real(sqrtm(Q_proj))*[ cos(angles) ; sin(angles) ];
    x_min = min([ellipse(1,:),x_min]);
    x_max = max([ellipse(1,:),x_max]);
    plot(ellipse(1,:),ellipse(2,:),'b-');
   
    margin_0 = (sqrt(y) - x(1)*[x_min,x_max])./x(2);
    plot([x_min,x_max],margin_0,'r-')

    hold off
 
end

warning('on','MATLAB:sqrtm:SingularMatrix')
