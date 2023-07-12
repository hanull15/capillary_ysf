%% model  of LWRE + HB + Slip + DCA (shear thinnin to slip model)
% develope not complete yet! (2022.06.15.)
%% 2022.03.01. last update
% clear;
% theta e 활용해서 rise data fitting 해서 DCA 같이 그려주는 코드 
% 2022,02.25. HB model and slip is added
% 2022.03.01. 
% 2step_2_2 eanables the N=1 calculation 
clear result*
win=3;
%% Basic input parameter  : 물질마다 바꾸기
global n ;
global k;
global Sy ;

% k=0.7; %recent measurement
%  n= 0.4;
% Sy=1; %[Pa]

%  sigma=0.063;
% rho =1000.7;
% % theta_zero=75.5225; % in degree
% % R=0.002;
%   alpha=5.46e-5; %  the slip prefactor [vel/stress]
% % % 
% % % beta0=1;
% X1=28.669;
% X2=3.5;

%% other parameters given 
NN=1000; % 100까지 sweep
g = 9.8;
mu=0.001; % [Pa s] solvent viscosity 

%% variables caculation 

s=1/n;
sigma_costhetazero = sigma*cos(theta_zero/180*3.141592);
H_e=2*sigma_costhetazero/R/rho/g;
cos_theta_zero=sigma_costhetazero/sigma;
theta_zero_rad = acos(cos_theta_zero);
cos_theta_zero=acos(cos_theta_zero)/pi*180;
% H=2*sigmatheta/rho/R/g;

t_v=(3+s)*(2^(n+1)*k/(rho*g)^(n+1)/R^(2*n+1))^s*sigma_costhetazero;
result_param= [Sy,k,n,sigma,rho,theta_zero,R,alpha,X,t_v];
%% load "h vs t" data

% [fileName,pathName] = uigetfile('*.xlsx','Select the capiallry rise excel file','YYMMDD_#');
% if pathName == 0
%     error('No file is selected!');
% end
% 
% file_with_loca=[pathName,fileName]
% M = readmatrix(file_with_loca, 'Filetype','spreadsheet');
% 
% open(file_with_loca);
% pause;

T=logspace(log(t_v/10000)/log(10),log(t_v*1000)/log(10),71)'; % 7 decades before and after t_v ... 10 pts / dec
t=M(:,1);
H=M(:,2);

%1st element of T is not zero -> add (0,0)
if T(1)~=0
    T=cat(1,[0],T);
    t=cat(1,[0],t);
    H=cat(1,[0],H);
end

num=length(T);
H=H./1000; % invert from mm unit to m unit

% T=T(1:(length(T)-sum(isnan(T))),:);
% H=H(1:(length(H)-sum(isnan(H))),:);
% T=[0;0.001;0.0502288160000000;0.100457633000000;0.200915266000000;0.401830532000000;0.803661063000000;1;2;3;4;5;6;7;9;12;15;18;23;28;36;44;56;69;87;108;136;169;212;265;331;414;517;646;808;1010;1262;1578;1972;2465;3081;3852;4815;6019;10000;20000;40000;80000;160000];


%% Calculating the constants for govering equation of SCA

c1=(3+1/n)^n*2*k/(R^(n+1));
c2=rho*g;
% c3=2*k/(R^n)*beta0;
c4=2*sigma_costhetazero/R;

%% contstants for equation _ HB with HD model

global D; 
D =[rho, rho*R, -2*sigma, 2, rho*g*R, cos(theta_zero_rad), -X2*R/sigma/alpha, -alpha, -R/(1/n+1)/(k^(1/n)), -X1*k/sigma/(R^(n-1)),-X1*R/sigma];


%% contstants for equation _ HB with HD model

global E; 
E =[rho, rho*R, -2*sigma, 2, rho*g*R, cos(theta_zero_rad), -beta0*k/sigma/(R^(n-1)), -alpha, -R/(1/n+1)/(k^(1/n))];


%% 각 DCA model numerical solving with 1st fitting

tspan = T;
% initial condition and the convergence option
y0  = [0.00001; 0.00001; 3.141592/180*60;100]; % 89 deg is the initial contact angle not to prevent singularity
yp0 = [0.000001 ; 0; 0; 10];
options = odeset('RelTol',1e-3,'AbsTol',[1e-7 1e-7 1e-3 1e-5]);

% HB model with HD model ----------------------------
[y0_1,yp0_1] = decic(@HB,0,y0,[1 0 0 0],yp0,[]);
[T_HD,H_HD] = ode15i(@HB,tspan,y0_1,yp0_1,options);

%% later regime with HD model ----------------------------
% find the initial condition for the later regime
% N=length(H_HD);
temp =(H_HD==real(H_HD));
temp = find(temp-1);
if isempty(temp)
    N=0
else
    N=temp(1)-1 % this is the minimum index for real part of H_HD
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
    
    result_H_2step(1:N,1:2) = H_HD(1:N,1:2)*1000; % height and velocity 
    result_H_2step(N+1:length(tspan),1:2)=H_HD_later(2:length(tspan)-N+1,1:2)*1000;
    
elseif N<=1
   [y0_later,yp0_later] = decic(@HB_later,0,y0,[1 0 0 0],yp0,[]);
    [T_HD_later,H_HD_later] = ode15i(@HB_later,tspan,y0_later,yp0_later,options);
    result_H_2step= H_HD_later(:,1:2)*1000;
end

%SCA  ------------------------------------------------
[T_SCA,H_SCA] = ode45(@(t,h) ((c4-c2*h)/(c1*h))^s,T,0.001*0.0000000001);

%% result for phase diagram 
phase_x = Sy/2/sigma_costhetazero*R ;                   %Bc number
phase_y = (3+s)/R*(k*alpha)^s*(sigma/4/pi/mu)^(1-s);    %V slip to V shaer thinning

power=gradient(log(result_H_2step(:,1)))./gradient(log(T)); % find the local power
plateau=islocalmin(power);
plateau=power(plateau);
TF=(plateau-0.01)>0;
plateau=plateau(TF);
% TF=(plateau<n);% the minimum power and maximum power at plateau.
% plateau=plateau(TF);
group=1;

if isempty(plateau)
    if phase_x>=0.01
        group=2;
    end
else
    group=3;
end 



  
 
 %% calculate the relative area of deviation of plateau 

sum_p=sum(log(result_H_2step(2:num,1)/H_e/1000));                                         % height is normalized with H_e , final height
Sy_original=Sy;
Sy=0.001;                                                                    % almost zero 


% HB model with HD model ----------------------------
[y0_np,yp0_np] = decic(@HB,0,y0,[1 0 0 0],yp0,[]);
[T_HD_np,H_HD_np] = ode15i(@HB,tspan,y0_np,yp0_np,options); 
 
temp_np =(H_HD_np==real(H_HD_np));
temp_np = find(temp_np-1);
if isempty(temp_np) || (temp_np(1))==2                                                         % every H_HD is real for ho Sy case
    result_H_2step_np(:,1:2) = H_HD_np(:,1:2)*1000;
else
    N_np=temp_np(1)-1                                                       % this is the minimum index for real part of H_HD
    H_HD_np=real(H_HD_np);
    
    % to get the initial condition for the later regime
    d_H_np =diff(H_HD_np(1:N_np,:));
    d_T_np=diff(tspan(1:N_np));
    yp_np=d_H_np./d_T_np;
    
    %run the calculation ----------------------------
    [y0_np_later,yp0_np_later] = decic(@HB_later,T_HD_np(N_np),H_HD_np(N_np,:)',[1 0 0 0],yp_np(N_np-1,:)',[]);
    [T_HD_later_np,H_HD_later_np] = ode15i(@HB_later,tspan(N_np:length(tspan)),y0_np_later,yp0_np_later,options);
    
    result_H_2step_np(1:N_np,1:2) = H_HD_np(1:N_np,1:2)*1000;
    result_H_2step_np(N_np+1:length(tspan),1:2)=H_HD_later_np(2:length(tspan)-N_np+1,1:2)*1000;
    
end
num_np=length(result_H_2step_np(:,1));
T_np=T(1:num_np);

sum_np=sum(log(result_H_2step_np(2:num_np,1)/H_e/1000));                                         % height is normalized with H_e , final height
Sy=Sy_original;

 result_phase = [phase_x,phase_y,group, plateau, sum_p, sum_np, sum_np-sum_p]
%  result_param= [Sy,k,n,sigma,rho,theta_zero,R,alpha,X];

%% error estimation 전에 한번 보여주기
figure(win)
subplot(1,2,1)
plot(t,H*1000,'o')
hold on 
% plot(T_MK,H_MK(:,1)*1000)
plot(T_HD(1:N),H_HD(1:N,1)*1000,'LineWidth',1.5)
plot(T_HD_later,H_HD_later(:,1)*1000,'LineWidth',1.5)
plot(T_SCA,H_SCA(:,1)*1000,'--')
plot(T_np,result_H_2step_np(:,1),'-.')
title('Entire linear');
xlabel('time (sec)')
ylabel('height (mm)')
legend('RAW DATA','HD','HD later','SCA','thininng','Location','southeast')
% xlim([0.01 2000])
ylim([0.1 12])

subplot(1,2,2)
loglog(t,H*1000,'o')
hold on
loglog(T_HD(1:N),H_HD(1:N,1)*1000,'LineWidth',1.5)
loglog(T_HD_later,H_HD_later(:,1)*1000,'LineWidth',1.5)
loglog(T_SCA,H_SCA(:,1)*1000,'--')
loglog(T_np,result_H_2step_np(:,1),'-.')
title('Entire log-log');
xlabel('time (sec)')
ylabel('height (mm)')
legend('RAW DATA','HD','HD later','SCA','thininng','Location','southeast')
% xlim([0.01 10000])
ylim([0.1 100])


%% define HB-capillary rise with out contact angle dynamics 
function modelHB = HB(~,y,yp) % DCA is shear thinning only 
global D
global n
global Sy

 
 modelHB =[ D(3)*cos(y(3))+D(4)*y(1)*y(4)+D(5)*y(1)
           yp(1)-y(2)
           y(3)^n*(D(6)-cos(y(3)))+D(10)*y(2)^n     %+D(11)*Sy yield stress effect how?
           y(2)+D(8)*y(4)+D(9)*(y(4)-Sy)^(1/n+1)/y(4)*(1-2*(1-Sy/y(4))^2/(1/n+3)-2*(1-Sy/y(4))/(1/n+2)*Sy/y(4))];
       
end

%% define HB-capillary rise with out contact angle dynamics 
function modelHB_MK = HB_MK(~,y,yp) % C는 constant 벡터 
global E
global n
global Sy

 modelHB_MK =[ E(3)*cos(y(3))+E(4)*y(1)*y(4)+E(5)*y(1)
           yp(1)-y(2)
           E(6)-cos(y(3))+E(7)*y(2)^n+Sy %the stress at the meniscus follows HB model 
           y(2)+E(8)*y(4)+E(9)*(y(4)-Sy)^(1/n+1)/y(4)*(1-2*(1-Sy/y(4))^2/(1/n+3)-2*(1-Sy/y(4))/(1/n+2)*Sy/y(4))];
       
       
end

%% define HB-capillary rise with out contact angle dynamics 
function modelHB_later = HB_later(~,y,yp) % DCA is slip only
global D
global n
% global Sy

 
 modelHB_later =[ D(3)*cos(y(3))+D(4)*y(1)*y(4)+D(5)*y(1)
           yp(1)-y(2)
           cos(y(3))*(D(6)-cos(y(3)))/(1-sin(y(3)))+D(7)*y(2)
           y(2)+D(8)*y(4)];
       
end