function S3=phase_3(t,s)
% %Motion equations phase 1: First free fly
global m_slv J_se J_re J_ce J_p R_s R_r R_p R_c r_g P1 P2...
    c_slv theta_g N N_h k ig kesi K_con D_con mu_con Jx1 Jx2...
    Jx3 Jx4 Jx Lrs1 Lrs2 Lrs3 Lrs4 Lcp1 Lcp2 Lcp3 Lcp4 K_d...
    data_save sum_e sum_e2 last_e
%% --------------------------Variables----------------------------
% S1=[
% 1)x_slv,
% 2)theta_ring,
% 3)theta_sun,
% 4)dx_slv,
% 5)dtheta_ring,
% 6)dtheta_sun,
%% --------------------------Dynamics-----------------------------R_s = 36.76e-3;% mm->m
r = R_s*1.9;

theta_g = 55*pi/180;

A1 = [-c_slv/m_slv 0 0 ;0 0 0 ; 0 0 0];
A = [zeros(3,3) eye(3);
    zeros(3,3) A1];

B = inv(Jx)*[1 0 R_r/2/R_c;0 1 R_s/2/R_c];
B1 = [1/m_slv 0 ;0 B(1,2); 0 B(2,2)];
B2 = [0 0 ;B(1,1) B(1,3); B(2,1) B(2,3)];
B_u = [zeros(3,2);B1];
B_w = [zeros(3,2);B2];
%% 碰撞发生条件
 F_con =0;
 f_con =0;


delta = 0;
ddelta = 0;
T_r = 0;
T_c = 0;

% theta_gr = 0;
% 如果在升档状态，theta_gr = thete_sun = x(3);
% 如果在降档状态，theta_gr = 0;


%% PID
% ref_p3 = 0.015;
% last_e=ref_p3-s(1);
% sum_e = sum_e+last_e;
% kp = 200000;
% kd = 1000;
% F_slv = PIDController(last_e,sum_e,kp,kd);
% function u=PIDController(e ,sum,kp,kd)
%     u=kp*e+kd*0;
% end

 F_slv =650;

%% collision heppen dynamics
colis = 0;
delta = 0;
F_con = 0;
f_con = 0;

if P1<s(1)&&s(1)<P2
    theta_gr = s(3);
    theta_sg = mod(-(s(2)-theta_gr),2*pi/N_h);
    dtheta = pi/N_h*(s(1)-P1)/0.004*(s(1)>=P1&&s(1)<=P2);

    if (theta_sg <= dtheta)&&theta_sg<pi/N_h&&theta_sg>0                                      %发生于接合齿圈下倒角
        colis = 1;
    elseif (theta_sg+dtheta>=2*pi/N_h)&&theta_sg>=pi/N_h&&theta_sg<2*pi/N_h %发生于接合齿圈上倒角
        colis = -1;
    else 
        colis= 0;
    end

        delta = ((dtheta-theta_sg)*sin(theta_g))*(colis==1) + ((dtheta+theta_sg-2*pi/N_h)*sin(theta_g))*(colis==-1);
        ddelta = (s(4)*sin(theta_g)+colis*(s(2)-s(3))*r*cos(theta_g))*(colis~=0);
        F_con = K_con*delta+D_con*ddelta;
        f_con = mu_con*F_con*(ddelta>0)+mu_con*F_con*(ddelta<0);
%         colis
%         F_con
        
end   
    F_slv = F_slv-(colis~=0)*(N_h*(f_con*cos(theta_g)+F_con*sin(theta_g)));
    % F_slv
    % s(1)
    % s(4)
    % delta
    T_r =  - colis*(f_con*sin(theta_g)+F_con*sin(theta_g))*N_h;         % slv == ring 认为向下为＋，认为slv > sun
    T_s =  0+ colis*(f_con*sin(theta_g)+F_con*sin(theta_g))*N_h;
    T_c = 0;

 Bu_3 = B_u;  
 Bw_3 = B_w;

S3 = A*s + Bu_3*[F_slv;T_s]+Bw_3*[T_r;T_c]+[0;0;0;-10*0.02*(s(4)>0)+10*0.02*(s(4)<0);0;0];

% data_save=[data_save;t,delta,F_slv,N_h*(f_con*cos(theta_g)+F_con*sin(theta_g)),delta,colis];
data_save=[data_save;t,delta,ddelta,F_slv,colis,(f_con*sin(theta_g)+F_con*sin(theta_g))*N_h,T_r,T_s];
end