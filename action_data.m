m_slv = 0.543;
J_se = 0.022; % 待定，应该加上转子和输出轴的总的转动惯量？？
J_re = 0.012;
J_ce = 1.528;
J_p = 0.01;

R_s = 36.76e-3;% mm->m
R_r = R_s*1.9;
R_p = (R_r-R_s)/2;
R_c = (R_r+R_s)/2;
r_g = R_r;

ig = 20;
c_slv = 1;
theta_g = 0.70;
N = 4;
N_h = 20;
kesi = 0.3; %泊松恢复系数

Jx1 = J_re+(R_r^2)/4/(R_c^2)*J_ce+N*(R_r^2)/4/(R_p^2)*J_p;
Jx2 = R_r*R_s/4/(R_c^2)*J_ce-N*R_r*R_s/4/(R_p^2)*J_p;
Jx3 = Jx2;
Jx4 = J_se+(R_s^2)/4/(R_c^2)*J_ce+N*(R_s^2)/4/(R_p^2)*J_p;
Jx = [Jx1 Jx2;Jx3 Jx4];


temp = Jx1*Jx4-Jx2^2;

Lrs1 = Jx4/temp;
Lrs2 = -Jx2/temp;
Lrs3 = Lrs2;
Lrs4 = Jx1/temp;

k = R_r/R_s; %一档变速器传动比

Lcp1 = (k*Lrs1+Lrs2)/(k+1);
Lcp2 = (k*Lrs2+Lrs4)/(k+1);
Lcp3 = (Lrs2-k*Lrs1)/(k-1);
Lcp4 = (Lrs4-k*Lrs2)/(k-1);

K_d = cos(theta_g)^2*r_g^2*(1/J_se+1/J_ce)+sin(theta_g)^2/m_slv;
kcon = 1.07e10;  %考虑齿轮耦合振动的换挡过程非线性动力学分析 Sui
mu_con = 0.3;

A1 = [-c_slv/m_slv 0 0 ;0 0 0 ; 0 0 0];
A = [A1];

B = inv(Jx)*[1 0 R_r/2/R_c;0 1 R_s/2/R_c];

B1 = [-1/m_slv 0 ;0 B(1,2); 0 B(2,2)];
B2 = [0 0 ;B(1,1) B(1,3); B(2,1) B(2,3)];
B_u = [B1];
B_w = [B2];

T_r = 0;

K1 = sin(theta_g)*(1+kesi)/m_slv/K_d;
K2 = cos(theta_g)*(1+kesi)*r_g/(J_re+J_ce)/K_d;
K3 = cos(theta_g)*(1+kesi)*r_g/(J_se)/K_d;
a_11 = 1-K1*sin(theta_g);
a_12 = K1*cos(theta_g)*r_g;
a_13 = -K1*cos(theta_g)*r_g;
a_21 = -K2*sin(theta_g);
a_22 = 1-K2*cos(theta_g)*r_g;
a_23 = K2*cos(theta_g)*r_g;
a_31 = K3*sin(theta_g);
a_32 = K3*cos(theta_g)*r_g;
a_33 = 1-K3*cos(theta_g)*r_g;
M = [a_11 a_12 a_13;a_21 a_22 a_23; a_31 a_32 a_33];
M1 = [a_11 -a_12 -a_13;-a_21 a_22 a_23; -a_31 a_32 a_33];

P_0 = 0.0005;
P_1 = 0.002;
P_2 = 0.004;
P_3 = 0.013;

        A_1 = [ zeros(3,3)];
        B1_1 = [0 0;0 0;0 (1+k)^2/(J_se*(1+k)^2+J_ce)];
        Bu_1 = [B1_1];
        B2_1 = [0 0;0 0;0 -(1+k)/(J_se*(1+k)^2+J_ce)];
        Bw_1 = [B2_1];

        A_2 = A;
        B1_2 = [-1/m_slv 0;0 0;0 (1+k)^2/(J_se*(1+k)^2+J_ce)];
        Bu_2 = [B1_2];
        B2_2 = [0 0;0 0;0 -(1+k)/(J_se*(1+k)^2+J_ce)];
        Bw_2 = [B2_2];


        A_3 = A;
        Bu_3 = B_u;  
        Bw_3 = B_w;

        A_4 = A;
        B1_4 = [-1/m_slv 0;0 0;0 1/(J_se+J_re+J_ce)];
        Bu_4 = [B1_4];
        B2_4 = [0 0;0 0;0 -1/(J_se+J_re+J_ce)];
        Bw_4 = [B2_4];

        A_5 = [zeros(3,3)];
        B1_5 = [0 0;0 0;0 1/(J_se+J_re+J_ce)];
        Bu_5 = [B1_5];
        B2_5 = [0 0;0 0;0 -1/(J_se+J_re+J_ce)];
        Bw_5 = [B2_5];

