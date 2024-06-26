 function main_p4
clc
clear all
close all

global m_slv J_se J_re J_ce J_p R_s R_r R_p R_c r_g P1 P2...
    c_slv theta_g N N_h k ig kesi K_con D_con mu_con Jx1 Jx2...
    Jx3 Jx4 Jx Lrs1 Lrs4 Lrs4 Lrs4 Lcp1 Lcp2 Lcp3 Lcp4 K_d...
    data_save sum_e last_e

%% --------- geometry ---------
R_s = 36.76e-3;   % 太阳轮
R_r = R_s*1.9;    % 齿圈
R_p = (R_r-R_s)/2;% 行星
R_c = (R_r+R_s)/2;% 桁架
r_g = R_r; 

P1 = 0.02-0.002;   % ->  结合齿顶点平行
P2 = 0.02; % -> 结合齿完全接合

c_slv = 1;         % 横向阻尼系数
theta_g = 0.70; % 结合齿角度
N = 4;       % 行星轮齿数
N_h = 25; % 接合套齿数

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
Lrs4 = -Jx2/(Jx1*Jx4-Jx2^2);
Lrs4 = Lrs4;
Lrs4 = Jx1/(Jx1*Jx4-Jx2^2);
Lcp1 = (k*Lrs1+Lrs4)/(k+1);
Lcp2 = (k*Lrs4+Lrs4)/(k+1);
Lcp3 = (Lrs4-k*Lrs1)/(k-1);
Lcp4 = (Lrs4-k*Lrs4)/(k-1);
K_d = cos(theta_g)^2*r_g^2*(1/J_se+1/J_ce)+sin(theta_g)^2/m_slv;


%% --------- collision ----------
kesi = 0.3; %泊松恢复系数
K_con = 3.0e5 ;
D_con = 0e1;
mu_con = 0.3;

sum_e = 0;
last_e = 0;




%%---------------------------phase 1 first free fly------------------------

%%------------------------- phase 2  Shift in low gear-------------------------

%%-------------------------phase 3 Neutral  gear  -------------------------


%%-------------------------phase 4 Collision in High gear  -------------------------
p=4
 % s = [ x_slv theta_slv theta_sun v_slv omega_slv omega_sun ]
s40=[0.0200002343998499	71.6718926536225	71.6718926536225	1.93230694596753	156.616304774300	113.616304774300];
t40=0;
t4f=0.05;


% data_save=[data_save;t,delta,ddelta,F_slv,colis,N_cons*N_h,f_cons,T_r,T_s];
data_save = [ 0,0,0,0,0,0,0,0,0];


options=@events4; %Create an optionsvariable
% [t1,s1]=ode45(@phase1,[t10 t1f],s10,options);
[t4,s4]=runge_kutta4(@phase_4,s40,1e-5,t40,t4f,options);

function [value,isterminal,direction]=events4(t4,s4)
value=(s4(6)-s4(5)); %Stops when s1(11)=1e-3
isterminal=1; %Stop after the first event (=0 to get all the events)
direction=0; % No matter which direction (+ -> - or - -> +)
% direction=-1;
end


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
save ./data/test_p4

figure(1);
subplot(4,1,1); %
plot(t4,s4(:,1)); hold on; 
text(t4(end),s4(end,1)*30/pi,'$$\ x_{slv} $$','Interpreter','latex','FontSize',14);
ylabel('$$\ (m) $$','FontSize',12,'Interpreter','latex');

figure(1);
subplot(4,1,2); %
plot(t4,s4(:,4)); hold on; 
text(t4(end),s4(end,1)*30/pi,'$$\ v_{slv} $$','Interpreter','latex','FontSize',14);
ylabel('$$\ (m/s) $$','FontSize',12,'Interpreter','latex');


subplot(4,1,3);%
plot(t4,s4(:,2)); hold on;
text(t4(end),s4(end,2),'$$\ \theta_{slv} $$','FontSize',14,'Interpreter','latex');

plot(t4,s4(:,3)); hold on;
text(t4(end),s4(end,3),'$$\ \theta_{sun} $$','FontSize',14,'Interpreter','latex');
ylabel('$$\ (rad) $$','FontSize',12,'Interpreter','latex');

subplot(4,1,4);%
plot(t4,s4(:,5)); hold on;
text(t4(end),s4(end,5),'$$\ \omega_{slv} $$','FontSize',14,'Interpreter','latex');

plot(t4,s4(:,6)); hold on;
text(t4(end),s4(end,6),'$$\ \omega_{sun} $$','FontSize',14,'Interpreter','latex');


% data_save=[data_save;t,delta,ddelta,F_slv,colis,(f_con*sin(theta_g)+F_con*sin(theta_g))*N_h,T_r,T_c];

figure(2);
subplot(4,2,1);plot(data_save(:,1),data_save(:,2));title("delta");
ylabel('$$\ (dm) $$','FontSize',12,'Interpreter','latex');
subplot(4,2,6);plot(data_save(:,1),data_save(:,5));title("colis");
% subplot(3,1,3);plot(data_save(:,1),data_save(:,5));title("v");
subplot(4,2,3);plot(data_save(:,1),data_save(:,4));title("F_{slv}");
ylabel('$$\ (m) $$','FontSize',12,'Interpreter','latex');
subplot(4,2,2);plot(data_save(:,1),data_save(:,3));title("ddelta");
subplot(4,2,4);plot(data_save(:,1),data_save(:,6));title("N_{cons}");
subplot(4,2,5);plot(data_save(:,1),data_save(:,7));title("f_{cons}");
subplot(4,2,7);plot(data_save(:,1),data_save(:,8));title("T_r");
% subplot(3,1,3);plot(data_save(:,1),data_save(:,5));title("v");
subplot(4,2,8);plot(data_save(:,1),data_save(:,9));title("T_s");
xlabel('$$\ Time\,(s) $$','FontSize',12,'Interpreter','latex');


figure(3);
subplot(2,1,1);
plot(t4,mod(s4(:,2)-s4(:,3),pi*2/N_h)); hold on;
text(t4(end),s4(end,6),'$$\ \delta \theta$$','FontSize',14,'Interpreter','latex');

subplot(2,1,2);
plot(t4,s4(:,5)-s4(:,6)); hold on;
text(t4(end),s4(end,6),'$$\ \delta \omega $$','FontSize',14,'Interpreter','latex');
% figure(4);
% plot(t4,J_re.*s4(:,5).*s4(:,5)+J_se.*s4(:,6).*s4(:,6)); hold on;
% text(t4(end),s4(end,6),'$$\ \omega $$','FontSize',14,'Interpreter','latex');

 end