%% capillary rise model ... LWRE + HB + Slip + DCA (shear thinnin -> slip model)
% 2022.03.01. last update
% written by hanull15@kaist.ac.kr

% For the first time user, just hit the 'Run' button. Please be sure calc.m
% and carbopol_0.12.mat file are in the same directory with this file.
% After running this file, you can see a window named as 'figure(1)'.

% "(optional)" sections, line 15~35 and 55~60, are for those who want to put
% their own data set. You can unlock those section simply by pressing ctl +
% t.

%% Initialization 
% clear 
% 
% % Load parameters and data for the caropol 0.12 wt. percent
% load('carbopol_0.12.mat')   

%%  (optional) Parameter setting : put your system variables manually
% % The default setting is Carbopol on glass capillary 

% global n ;
% global k;
% global Sy ;
% 
% % % H-B model rheology parameters
% k=3.5;              % viscosity at shear rate = 1 1/s [Pa s^n]
% n= 0.41;             % shear thinning power [1]
% Sy=3.5;               % yield stress [Pa]
% 
% sigma=0.0315;        % surface tension [N/m]
% rho =1000.8;        % density [kg/m^3]
% theta_zero=54; % equilibrium contact angle [degree]
% R=0.00085;           % capillary radius [m]
% l_s=1;            % slip layer thickness [nm]
% 
% % DCA constants
% X1=3.5;          % for the early regime 
% X2=0.6;             % for the later regime

%% other parameters given 
NN=1000; 
g = 9.8;            % gravitational constant
mu=0.001;           % [Pa s] solvent viscosity 

%% variables caculation 

alpha=l_s*10^(-9)/mu;      % slip prefactor [vel/stress], l_s = alpha*mu 
s=1/n;
sigma_costhetazero = sigma*cos(theta_zero/180*3.141592);
H_e=2*sigma_costhetazero/R/rho/g;
cos_theta_zero=sigma_costhetazero/sigma;
theta_zero_rad = acos(cos_theta_zero);
cos_theta_zero=acos(cos_theta_zero)/pi*180;

t_v=(3+s)*(2^(n+1)*k/(rho*g)^(n+1)/R^(2*n+1))^s*sigma_costhetazero;         % viscous regime time scale

%% make struct variable from pararmeters 
param= struct ('Sy',Sy, 'n', n, 'k',k, 'sigma', sigma, 'rho',rho,'theta_zero', theta_zero, 'R',R, 'alpha', alpha,'X1',X1, 'X2',X2,'t_v',t_v,'Sy_fit',[]);

%% (optional) Load your own rise curve 
% unlock is section when you want to load your own experimental capillary rise data

%t = load('capillaryRise_time.mat')
%H = load('capillaryRise_height.mat')

%% set time scale based on t_v

T=logspace(log(t_v/10000)/log(10),log(t_v*1000)/log(10),71)'; % 7 decades before and after t_v ... 10 pts / dec


%1st element of T is not zero -> add (0,0)
if T(1)~=0
    T=cat(1,[0],T);
end

num=length(T);

%% Calculating the constants for govering equation of SCA

c1=(3+1/n)^n*2*k/(R^(n+1));
c2=rho*g;
c4=2*sigma_costhetazero/R;

%% Numerical solver : with yield stress 
[N, result_H_2step, H_HD, H_HD_later]=calc(T, param);
sum_p=sum(log(result_H_2step(2:num,1)/H_e/1000));                                           % height is normalized with H_e , final height

%% Numerical solver : without yield stress 

% no DCA  no slip no yield stress
[T_SCA,H_SCA] = ode45(@(t,h) ((c4-c2*h)/(c1*h))^s,T,0.001*0.0000000001);


% only without yield stress 
Sy_original=Sy;
param.Sy=0.001;                                                                                % Sy ~ zero 
[N_np, result_H_2step_np, H_HD_np, H_HD_later_np]=calc(T,param);

num_np=length(result_H_2step_np(:,1));
T_np=T(1:num_np);

sum_np=sum(log(result_H_2step_np(2:num_np,1)/H_e/1000));                                      % height is normalized with H_e , final height
Sy=Sy_original;

result_phase = [ sum_p, sum_np, sum_np-sum_p]                                                 % Calculate delta 


%% error estimation 전에 한번 보여주기
figure(1)
subplot(1,2,1)
plot(t,H*1000,'o')
hold on 
% plot(T_MK,H_MK(:,1)*1000)
plot(T(1:N),H_HD(1:N,1)*1000,'LineWidth',1.5)
plot(T(N:length(T)),H_HD_later(:,1)*1000,'LineWidth',1.5)
plot(T_SCA,H_SCA(:,1)*1000,'--')
plot(T_np,result_H_2step_np(:,1),'-.')
title('Entire linear');
xlabel('time (sec)')
ylabel('height (mm)')
legend('RAW DATA','HD','HD later','SCA','thininng','Location','southeast')
% xlim([0.01 2000])
% ylim([0.1 6])

subplot(1,2,2)
loglog(t,H*1000,'o')
hold on
loglog(T(1:N),H_HD(1:N,1)*1000,'LineWidth',1.5)
loglog(T(N:length(T)),H_HD_later(:,1)*1000,'LineWidth',1.5)
loglog(T_SCA,H_SCA(:,1)*1000,'--')
loglog(T_np,result_H_2step_np(:,1),'-.')
title('Entire log-log');
xlabel('time (sec)')
ylabel('height (mm)')
legend('RAW DATA','HD','HD later','SCA','thininng','Location','southeast')
% xlim([0.01 10000])
% ylim([0.1 10])

  
 