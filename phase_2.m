function S2=phase_2(t,s)
% %Motion equations phase 1: First free fly
global m_slv J_se J_re J_ce J_p R_s R_r R_p R_c r_g P1 P2...
    c_slv theta_g N N_h k ig kesi K_con D_con mu_con Jx1 Jx2...
    Jx3 Jx4 Jx Lrs1 Lrs2 Lrs3 Lrs4 Lcp1 Lcp2 Lcp3 Lcp4 K_d...
    data_save sum_e
%% --------------------------Variables----------------------------
% S1=[
% 1)x_slv,
% 2)theta_ring,
% 3)theta_sun,
% 4)dx_slv,
% 5)dtheta_ring,
% 6)dtheta_sun,
%% --------------------------Dynamics-----------------------------
 A1 = [-c_slv/m_slv 0 0 ;0 0 0 ; 0 0 0];
 A = [zeros(3,3) eye(3);
      zeros(3,3) A1]; 
 A_2 = A;
 B1_2 = [-1/m_slv 0;0 0;0 (1+k)^2/(J_se*(1+k)^2+J_ce)];
 Bu_2 = [zeros(3,2);B1_2];
 B2_2 = [0 0;0 0;0 -(1+k)/(J_se*(1+k)^2+J_ce)];
 Bw_2 = [zeros(3,2);B2_2];

ref_p2=0.018; %参考输入
u=0; %控制量存储空间预设

e=ref_p2-s(1);
sum_e = sum_e+e;


F = PIDController(e,sum_e);
% F
% s(1)
function u=PIDController(e ,sum)
    u=-3000*e-0*sum;
end


S2 = A_2*s + Bu_2*[F;0]+Bw_2*[0;0]+[0;0;0;-10*0.02*(s(4)>0)+10*0.02*(s(4)<0);0;0];
%%
% data_save=[data_save;t];
% data_save=[data_save;t,Ffork,Fslspx,Fslspy,Fxsrsl,Nsr_str,h,Tsi_sr,Tstr_sr,tempd,S,fc,stage,Fyisl,Tsl_i];
end