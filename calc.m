%% calculating the early and later capillary rise 
function [N, result_H_2step, H_HD, H_HD_later]=calc(T,param)
global D; 
global n
global Sy
global s;

n=param.n;
Sy=param.Sy;
s=1+1/n;
D =[param.rho, param.rho*param.R, -2*param.sigma, 2, param.rho*9.8*param.R, cos(param.theta_zero/180*3.141592), -param.X2*param.R/param.sigma/param.alpha, -param.alpha, -param.R*(1/param.k)^(1/param.n), -param.X1*param.k/param.sigma/(param.R^(param.n-1)),-param.X1*param.R/param.sigma];



%% early calculation
tspan = T;

% initial condition and the convergence option
y0  = [0.0001; 0.001; 3.141592/180*60;100]; % 89 deg is the initial contact angle to prevent singularity
yp0 = [0.001 ; 0; 0; 10];
options = odeset('RelTol',1e-3,'AbsTol',[1e-7 1e-7 1e-3 1e-5]);

% HB model with HD model ----------------------------
[y0_1,yp0_1] = decic(@HB,0,y0,[1 0 0 0],yp0,[]);
[~,H_HD] = ode15i(@HB,tspan,y0_1,yp0_1,options);

%% later calculation
% find the initial condition for the later regime

temp =(H_HD==real(H_HD));
temp = find(temp-1);
if isempty(temp)
    N=0;
else
    N=temp(1)-1; % this is the minimum index for real part of H_HD
end

H_HD=real(H_HD);

if N>1
    % to get the initial condition for the later regime
    d_H =diff(H_HD(1:N,:));
    d_T=diff(tspan(1:N));
    yp=d_H./d_T;
    
    %run the calculation ----------------------------
    [y0_later,yp0_later] = decic(@HB_later,T(N),H_HD(N,:)',[1 0 0 0],yp(N-1,:)',[]);
    [T_HD_later,H_HD_later] = ode15i(@HB_later,tspan(N:length(tspan)),y0_later,yp0_later,options);
    
    result_H_2step(1:N,:) = H_HD(1:N,:); % height and velocity 
    if N==length(tspan)-1
    result_H_2step(N+1,:)=H_HD_later(length(T_HD_later),:);
    else
        result_H_2step(N+1:N-1+length(T_HD_later),:)=H_HD_later(2:length(T_HD_later),:);
    end
        
elseif N<=1
   [y0_later,yp0_later] = decic(@HB_later,0,y0,[1 0 0 0],yp0,[]);
    [~,H_HD_later] = ode15i(@HB_later,tspan,y0_later,yp0_later,options);
    result_H_2step= H_HD_later(:,:);
end
result_H_2step(:,1:2)=result_H_2step(:,1:2).*1000;



end 

%% define HB-capillary rise with out contact angle dynamics 
function modelHB = HB(~,y,yp) % DCA is shear thinning only 
global D
global n
global Sy
global s

 
 modelHB =[ D(3)*cos(y(3))+D(4)*y(1)*y(4)+D(5)*y(1)
           yp(1)-y(2)
           y(3)^n*(D(6)-cos(y(3)))+D(10)*y(2)^n     %+D(11)*Sy yield stress effect how?
           y(2)+D(8)*y(4)+D(9)*(y(4))^(-3)*(y(4)-Sy)^s*(Sy^2/s+(y(4)-Sy)^2/(1/n+3)+2*Sy*(y(4)-Sy)/(1/n+2))];
       
end

%% define HB-capillary rise with out contact angle dynamics 
function modelHB_later = HB_later(~,y,yp) % DCA is slip only
global D


 
 modelHB_later =[ D(3)*cos(y(3))+D(4)*y(1)*y(4)+D(5)*y(1)
           yp(1)-y(2)
           cos(y(3))*(D(6)-cos(y(3)))/(1-sin(y(3)))+D(7)*y(2)
           y(2)+D(8)*y(4)];
       
end