%% model  of LWRE + HB + Slip + DCA (shear thinnin to slip model)
%% 2022.08.22. last update

clear result*
win=3;
step=0.1;
%% 10 basic input parameters : material dependent

k=10.9; %recent measurement
n= 0.36;
Sy=10.2; %[Pa]

sigma=0.0345;
rho =1001.2;
theta_zero=14.46; % in degree
R=0.00085;
alpha=1e-4; %  the slip prefactor [vel/stress]

X1=2.3;
X2=0.039;


%% other parameters given 
% NN=1000; % 100까지 sweep
g = 9.8;
mu=0.001; % [Pa s] solvent viscosity 

%% variables caculation 
s=1/n;
sigma_costhetazero = sigma*cos(theta_zero/180*3.141592);
H_e=2*sigma_costhetazero/R/rho/g;
cos_theta_zero=sigma_costhetazero/sigma;
theta_zero_rad = acos(cos_theta_zero);
cos_theta_zero=acos(cos_theta_zero)/pi*180;
t_v=(3+s)*(2^(n+1)*k/(rho*g)^(n+1)/R^(2*n+1))^s*sigma_costhetazero;

%% make struct variable from pararmeters 
param= struct ('Sy',Sy, 'n', n, 'k',k, 'sigma', sigma, 'rho',rho,'theta_zero', theta_zero, 'R',R, 'alpha', alpha,'X1',X1, 'X2',X2,'t_v',t_v,'Sy_fit',[]);

%% load experimental data
% % 
% [fileName,pathName] = uigetfile('*.txt','Select the capiallry rise excel file','YYMMDD_#');
% if pathName == 0
%     error('No file is selected!');
% end
% 
% file_with_loca=[pathName,fileName]
% M = readmatrix(file_with_loca);


T=logspace(log(t_v/10000)/log(10),log(t_v*1000)/log(10),71)'; % 7 decades before and after t_v ... 10 pts / dec
t=M(:,1);
H=M(:,2);

%1st element of T is not zero -> add (0,0)
if t(1)~=0
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

%% 각 DCA model numerical solving with 1st fitting
[N, result_H_2step, H_HD, H_HD_later]=calc(T, param);

%SCA  ------------------------------------------------
[T_SCA,H_SCA] = ode45(@(t,h) ((c4-c2*h)/(c1*h))^s,T,0.001*0.0000000001);

 
 %% calculate the relative area of deviation of plateau 

sum_p=sum(log(result_H_2step(2:num,1)/H_e/1000));                                         % height is normalized with H_e , final height
Sy_original=Sy;
param.Sy=0.01;                                                                    % almost zero 


% HB model with HD model ----------------------------
[N_np, result_H_2step_np, H_HD_np, H_HD_later_np]=calc(T,param);

num_np=length(result_H_2step_np(:,1));
T_np=T(1:num_np);

sum_np=sum(log(result_H_2step_np(2:num_np,1)/H_e/1000));                                         % height is normalized with H_e , final height
Sy=Sy_original;

 result_phase = [ sum_p, sum_np, sum_np-sum_p]


%% error estimation 전에 한번 보여주기
figure(win)
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
ylim([0.1 12])

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
ylim([0.1 100])

%% yield stress fitting (22.08.07.) 
% chagne the ode solver in two steps
% calculate squared relative error varying the yield stress! 

param.Sy=Sy_original;
round_fit = round(Sy_original/step*5); % fitting range is until the five times of the original Sy
n_data=length(t);
param.Sy=step;        %starting Sy is about 0
error=zeros(round_fit,2); % first is Sy and second is error
turn_fit=1;
while turn_fit<=round_fit
    [~, H_fit,~, ~]=calc(t,param);
%     error_this = ((H(:,1)-H_fit(:,1)./1000)./H(:,1)).^2; % squared relative error for each time step
    error_this = ((H(:,1)-H_fit(:,1)./1000)).^2; % squared relative error for each time step
    error(turn_fit,1)=param.Sy;
    error(turn_fit,2)=sum(error_this(2:length(error_this)));
    param.Sy=param.Sy+step;
    turn_fit=turn_fit+1;
end
[fit_res,fit_ind]=min(error(:,2));
Sy_fit=error(fit_ind,1);
param.Sy_fit=Sy_fit;

param.Sy=Sy_fit;
[~, H_fit, ~,~]=calc(t, param);
param.Sy=Sy_original;
    
%% After error estimation 
figure(win+2)
subplot(1,2,1)
plot(t,H*1000,'o')
hold on 
plot(t,H_fit(:,1),'-','Linewidth',1.5)
plot(T_SCA,H_SCA(:,1)*1000,'--')
plot(T_np,result_H_2step_np(:,1),'-.')
title('Entire linear');
xlabel('time (sec)')
ylabel('height (mm)')
legend('RAW DATA','FIT','SCA','thininng','Location','southeast')
% xlim([0.01 2000])
ylim([0.1 12])

subplot(1,2,2)
loglog(t,H*1000,'o')
hold on
plot(t,H_fit(:,1),'-','Linewidth',1.5)
loglog(T_SCA,H_SCA(:,1)*1000,'--')
loglog(T_np,result_H_2step_np(:,1),'-.')
title('Entire log-log');
xlabel('time (sec)')
ylabel('height (mm)')
legend('RAW DATA','FIT','SCA','thininng','Location','southeast')
% xlim([0.01 10000])
ylim([0.1 100])

