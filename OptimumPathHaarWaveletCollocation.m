clear all;
close all;



j = 2;
v = 1;
m = 2^(j+1);
pmat = 1/2; 
Hmat = 1;
for ind=1:j
    Hmat = Hm(ind);
    pmat = Pm(pmat,Hmat,2^(ind+1)); 
end





P_matrix = pmat;
H_matrix =Hm(3);




global nmc;


nmc = 2*2^(j+1);
Cx = [sym('c_x',[1 nmc])];
Cy = [sym('c_y',[1 nmc])];
Cu = [sym('c_u',[1 nmc])];



% x_dot = v*cos(u)
% y_dot = v*sin(u)
% (x-3)^2 + (y-3)^2 >= (0.5)^2


% Cx * hk = x_dot_dot
% Cx * P * hk = x_dot
% Cx * P * P *hk = x

% Cy * hk = y_dot_dot
% Cy * P * hk = y_dot
% Cy * P * P * hk = y_dot

% Cu * hk = u



% Cx * P * hk = v * cos(Cu * hk)
% Cy * P * hk = v * sin(Cu * hk) 



% x_dot_lhs = Cx * P_matrix * H_matrix;
% x_dot_rhs = v * cos(Cu * H_matrix);
% 
% y_dot_lhs = Cy * P_matrix * H_matrix;
% y_dot_rhs = v * sin(Cu * H_matrix);
% 
% 
% 
% x_lhs = Cx *P_matrix * P_matrix * H_matrix(:,nmc);
% x_rhs = 5;
% 
% y_lhs = Cy *P_matrix * P_matrix * H_matrix(:,nmc);
% y_rhs = 5;




% circle condition


%*((Cx * P_matrix * P_matrix *H_matrix)-3);


% c_temp = ((0.5)^2 -(((Cx * P_matrix * P_matrix *H_matrix)-3).^2 +((Cy * P_matrix * P_matrix *H_matrix)-3).^2)).';
% ceq_temp = [(x_dot_lhs - x_dot_rhs).';(y_dot_lhs - y_dot_rhs).'];

% Y = [Cx,Cy,Cu]
% circlecon(Y) = [c_temp,ceq_temp]
nonlcon = @circlecon;


% global c_temp,ceq_temp;
A = [];



% Aeq and beq   that is linear equality constraint

dumdum = ((P_matrix * P_matrix * H_matrix(:,nmc)).' );

linearAmat_x = zeros(1,3*nmc);

for z= 1:nmc
    linearAmat_x(z) = dumdum(z);
end


linearAmat_y = zeros(1,3*nmc);

for z=1:nmc
    linearAmat_y(nmc+z) = dumdum(z);
end

Aeq = [linearAmat_x;linearAmat_y];
Beq = [5;5];








AHHH = zeros(1,3*nmc);


A = [];
b = [];
LB = zeros(1,nmc*3);
UB = [];
% syms X;
% X = [Cx,Cy,Cu];


temp_integral = @(X) sum((X(1:nmc) * H_matrix).^2 +(X(nmc+1:2*nmc) *H_matrix).^2) ;
%ti = matlabFunction(temp_integral);

%chaddha = ones(1,nmc*3);
%temp_integral(chaddha)
options = optimoptions('fmincon','Display','iter-detailed','Algorithm','sqp');
result = fmincon(temp_integral,AHHH,A,b,Aeq,Beq,LB,UB,nonlcon,options)

for k = 1:nmc
        discrete_time(k) = (k-0.5)/nmc;
end

xt = result(1:nmc)*P_matrix*P_matrix*H_matrix;

yt = result(nmc+1:nmc*2)*P_matrix*P_matrix*H_matrix;


plot(discrete_time,xt);
plot(discrete_time,yt);

plot(xt,yt);
scatter(xt,yt);





function Hmatrix = Hm(J)
    syms x u y c t k;
    
    g = 10;
    y0 = 0;
    x0 = 0;
    
    nmc = 2^(J+1);
    
    discrete_time = zeros(1,nmc);
    
    for k = 1:nmc
        discrete_time(k) = (k-0.5)/nmc;
    end
    
    
    
    
    % x_dot =  sqrt(2*g*(y0 - y))*cos(u);
    % y_dot =  sqrt(2*g*(y0 - y))*sin(u);
    % 
    % 
    % Cx = [sym('c_x',[1 nmc])];
    % Cy = [sym('c_y',[1 nmc])];
    % Cu = [sym('c_u',[1 nmc])];
    
    
    
    
    
    
    Hmatrix = zeros(nmc,nmc);
    Hmatrix(:,1) = ones(nmc,1);
    
    for j = 0:J
    
        m = 2^j;
    
        for k =0:m-1
    
            i = m+k;
            
            fun(t) = piecewise((k/m)<= t <(k+0.5)/m,1,(k+0.5)/m<= t <(k+1)/m,-1,0);
            %fun(0.3)
            
    
            for q = 1:nmc
                
                time = discrete_time(q);
                %hmatrices(1,i,q) = fun(time);
                Hmatrix(q,i+1) = fun(time);
    
            end
        end
    end

    Hmatrix = Hmatrix.';
end



% rohan = Cx*(Hmatrix(1,:).')
function p = Pm(p2, H,m)
    p = zeros(m*2,m*2);
    p(1:m,1:m) = p2;
    p(1:m,m+1:end) = -1/(m)*H;
    p(m+1:end,1:m) = 1/m*inv(H);
    p(m+1:end,m+1:end) = 0;
end



function [c,ceq] = circlecon(X)

j = 2;
v = 1;
m = 2^(j+1);
pmat = 1/2; 
Hmat = 1;
for ind=1:j
    Hmat = Hm(ind);
    pmat = Pm(pmat,Hmat,2^(ind+1)); 
end





P_matrix = pmat;
H_matrix =Hm(3);







nmc = 2*2^(j+1);
Cx = [sym('c_x',[1 nmc])];
Cy = [sym('c_y',[1 nmc])];
Cu = [sym('c_u',[1 nmc])];


nmc = 16;
Cx = X(1:nmc);
Cy = X(nmc+1:nmc*2);
Cu = X(nmc*2+1:end);


x_dot_lhs = Cx * P_matrix * H_matrix;
x_dot_rhs = v * cos(Cu * H_matrix);
y_dot_lhs = Cy * P_matrix * H_matrix;
y_dot_rhs = v * sin(Cu * H_matrix);
x_lhs = Cx *P_matrix * P_matrix * H_matrix(:,nmc);
x_rhs = 5;
y_lhs = Cy *P_matrix * P_matrix * H_matrix(:,nmc);
y_rhs = 5;

c = ((0.5)^2 -(((Cx * P_matrix * P_matrix *H_matrix)-3).^2 +((Cy * P_matrix * P_matrix *H_matrix)-3).^2)).';
ceq = [(x_dot_lhs - x_dot_rhs).';(y_dot_lhs - y_dot_rhs).'];

% c = c_temp;     
% ceq = ceq_temp;
end

function stop = outfun(x,optimValues,state)
stop = toc>T;
end

