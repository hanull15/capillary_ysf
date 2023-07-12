%% making phase diagram 
% read txt sheet compromised with parameters 
% claculate all the parameters for each parameter set 
% make the contour map using Hy/He and Vs/Vthinning 

g = 9.8;
mu=0.001; % [Pa s] solvent viscosity 
%% Bring
[fileName,pathName] = uigetfile('*.txt','Select the capiallry rise excel file','YYMMDD_#');
if pathName == 0
    error('No file is selected!');
end

file_with_loca=[pathName,fileName]
M = readmatrix(file_with_loca);

[num, m] = size(M); %n is the number of the phase diagram data points. 

if m~=10 
    error('the number of columns must be 10 !');
end
if num==0
    error('There is no data')
end

phase = zeros(num,3);

%% calculate necessay values. 

i=1;
while i<=num
param_p= struct ('Sy',M(i,1), 'n', M(i,2), 'k',M(i,3), 'sigma', M(i,4), 'rho',M(i,5),'theta_zero', M(i,6), 'R',M(i,7), 'alpha', M(i,8),'X1',M(i,9), 'X2',M(i,10),'t_v',[]);
s=1/param_p.n;
sigma_costhetazero = param_p.sigma*cos(param_p.theta_zero/180*3.141592);
H_e=2*sigma_costhetazero/param_p.R/param_p.rho/g;

t_v=(3+s)*(2^(param_p.n+1)*param_p.k/(param_p.rho*g)^(param_p.n+1)/param_p.R^(2*param_p.n+1))^s*sigma_costhetazero;
param_p.t_v=t_v;

T=logspace(log(t_v/10000)/log(10),log(t_v*1000)/log(10),71)'; % 7 decades around t_v: this term must be dependent on the parameters. 
[N, result_H_2step, H_HD, H_HD_later]=calc(T,param_p);
num_p=length(result_H_2step(:,1));
sum_p=sum(log(result_H_2step(2:num_p,1)/H_e/1000));   


Sy_original=param_p.Sy;
param_p.Sy=0.001; 
% HB model with HD model ----------------------------
[N_np, result_H_2step_np, H_HD_np, H_HD_later_np]=calc(T,param_p);


num_np=length(result_H_2step_np(:,1));
% T_np=T(1:num_np);

sum_np=sum(log(result_H_2step_np(2:num_np,1)/H_e/1000));                                         % height is normalized with H_e , final height
param_p.Sy=Sy_original;

phase(i,1)=  1/(1+param_p.Sy*2/param_p.rho/g/param_p.R);               % x axis
phase(i,2)=  8*param_p.alpha*param_p.Sy/param_p.R*(param_p.k/param_p.Sy)^s;               % y axis
phase(i,3)= sum_np-sum_p;

i=i+1;
end

%% interpolation

% study the range of the region 
min_x = round(log10(min(phase(:,1)))-0.5);
max_x = round(log10(max(phase(:,1)))+0.5);
min_y = round(log10(min(phase(:,2)))-0.5);
max_y = round(log10(max(phase(:,2)))+0.5);

phase_log(:,1:2)=log10(phase(:,1:2));

[Xq,Yq] = meshgrid(min_x:0.05:max_x,min_y:0.05:max_y);
Vq = griddata(phase_log(:,1),phase_log(:,2),phase(:,3),Xq,Yq,'linear'); %%%% ERROR LINE -----------------------------------------------------------------------------

% Vq = interp2(phase_log(:,1),phase_log(:,2),phase(:,3),5);
%% plot interpolated figure 
figure(100); surf(Xq,Yq,Vq); % axis is in log10 unit
hold on;
plot3(phase_log(:,1),phase_log(:,2),phase(:,3),'o');