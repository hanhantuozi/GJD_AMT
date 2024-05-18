function main_p2
clc
clear all
close all

global m_slv J_se J_re J_ce J_p R_s R_r R_p R_c r_g P1 P2...
    c_slv theta_g N N_h k ig kesi K_con D_con mu_con Jx1 Jx2...
    Jx3 Jx4 Jx Lrs1 Lrs2 Lrs3 Lrs4 Lcp1 Lcp2 Lcp3 Lcp4 K_d...
    data_save sum_e last_e

%% --------- geometry ---------
R_s = 36.76e-3;   % 太阳轮
R_r = R_s*1.9;    % 齿圈
R_p = (R_r-R_s)/2;% 行星
R_c = (R_r+R_s)/2;% 桁架
r_g = R_r; 

P1 = 0.02-0.002;
P2 = 0.02;

c_slv = 1;
theta_g = 0.70;
N = 4;
N_h = 25;

k = R_r/R_s;%一档变速器传动比
ig = 20; % 末端减速比 

%% --------- mass --------- 
m_slv = 2.543;
J_se = 0.122; % 待定，应该加上转子和输出轴的总的转动惯量？？
J_re = 0.112;
J_ce = 1.528;
J_p = 0.02;

% for N stage
Jx1 = J_re+(R_r^2)/4/(R_c^2)*J_ce+N*(R_r^2)/4/(R_p^2)*J_p;
Jx2 = R_r*R_s/4/(R_c^2)*J_ce-N*R_r*R_s/4/(R_p^2)*J_p;
Jx3 = Jx2;
Jx4 = J_se+(R_s^2)/4/(R_c^2)*J_ce+N*(R_s^2)/4/(R_p^2)*J_p;
Jx = [Jx1 Jx2;Jx3 Jx4];
Lrs1 = Jx4/(Jx1*Jx4-Jx2^2);
Lrs2 = -Jx2/(Jx1*Jx4-Jx2^2);
Lrs3 = Lrs2;
Lrs4 = Jx1/(Jx1*Jx4-Jx2^2);
Lcp1 = (k*Lrs1+Lrs2)/(k+1);
Lcp2 = (k*Lrs2+Lrs4)/(k+1);
Lcp3 = (Lrs2-k*Lrs1)/(k-1);
Lcp4 = (Lrs4-k*Lrs2)/(k-1);
K_d = cos(theta_g)^2*r_g^2*(1/J_se+1/J_ce)+sin(theta_g)^2/m_slv;


%% --------- collision ----------
kesi = 0.3; %泊松恢复系数
K_con = 1.07e7;
D_con = 1e2;
mu_con = 0.3;

sum_e = 0;
last_e = 0;
data_save=[];



s0=[0	5.21323383394755	5.23598775598298	0	104.264676678951	104.719755119660];
t0=0; %Simulation starting time for phase 1
t1f=0.05; %Simulation final time for phase 1

%%---------------------------phase 1 first free fly------------------------

%%------------------------- phase 2 -------------------------
p=2
s20=s0;
t20=0;
t2f=0.6;

options=@events2; %Create an optionsvariable
% [t1,s1]=ode45(@phase1,[t10 t1f],s10,options);
[t2,s2]=runge_kutta4(@phase_2,s20,1e-2,t20,t2f,options);
function [value,isterminal,direction]=events2(t2,s2)
value=s2(1)-0.4; %Stops when s1(11)=1e-3
isterminal=0; %Stop after the first event (=0 to get all the events)
direction=-1; % No matter which direction (+ -> - or - -> +)
% direction=-1;
end

%%-------------------------phase 3 Neutral  -------------------------



%------------------------------Plots Set----------------------------------
co = [0 0 1;
      1 0 0;
      0 0.5 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25];
set(groot,'defaultAxesColorOrder',co);
set(0,'DefaultLineLineWidth',1.5);
%------------------------------Plots----------------------------------
save ./data/test_p2

figure(1);
subplot(3,1,1); %
plot(t2,s2(:,1)); hold on; 
text(t2(end),s2(end,1)*30/pi,'$$\ x_{slv} $$','Interpreter','latex','FontSize',14);
ylabel('$$\ (m) $$','FontSize',12,'Interpreter','latex');

subplot(3,1,2);%
plot(t2,s2(:,2)); hold on;
text(t2(end),s2(end,2),'$$\ \theta_{slv} $$','FontSize',14,'Interpreter','latex');

plot(t2,s2(:,3)); hold on;
text(t2(end),s2(end,3),'$$\ \theta_{sun} $$','FontSize',14,'Interpreter','latex');
ylabel('$$\ (rad) $$','FontSize',12,'Interpreter','latex');

subplot(3,1,3);%
plot(t2,s2(:,5)); hold on;
text(t2(end),s2(end,5),'$$\ \omega_{slv} $$','FontSize',14,'Interpreter','latex');

plot(t2,s2(:,6)); hold on;
text(t2(end),s2(end,6),'$$\ \omega_{sun} $$','FontSize',14,'Interpreter','latex');

% % figure(1);
% subplot(3,1,1); %
% plot(t1,s1(:,1)); hold on; 
% text(t1(end),s1(end,1)*30/pi,'$$\ x_{slv} $$','Interpreter','latex','FontSize',14);
% ylabel('$$\ (m) $$','FontSize',12,'Interpreter','latex');
% 
% subplot(3,1,2);%
% plot(t1,s1(:,2)); hold on;
% text(t1(end),s1(end,2),'$$\ \theta_{slv} $$','FontSize',14,'Interpreter','latex');
% 
% plot(t1,s1(:,3)); hold on;
% text(t1(end),s1(end,3),'$$\ \theta_{sun} $$','FontSize',14,'Interpreter','latex');
% ylabel('$$\ (rad) $$','FontSize',12,'Interpreter','latex');
% 
% subplot(3,1,3);%
% plot(t1,s1(:,5)); hold on;
% text(t1(end),s1(end,5),'$$\ \omega_{slv} $$','FontSize',14,'Interpreter','latex');
% 
% plot(t1,s1(:,6)); hold on;
% text(t1(end),s1(end,6),'$$\ \omega_{sun} $$','FontSize',14,'Interpreter','latex');
figure(2);
subplot(3,1,1);plot(data_save(:,1),data_save(:,2));title("F");
% subplot(3,1,2);plot(data_save(:,1),data_save(:,6));title("colis");
% subplot(3,1,3);plot(data_save(:,1),data_save(:,5));title("v");
% data_save=[data_save;t,delta,F_slv,s(1),s(4)];
% 
ylabel('$$\ (rad/s) $$','FontSize',12,'Interpreter','latex');

xlabel('$$\ Time\,(s) $$','FontSize',12,'Interpreter','latex');
end