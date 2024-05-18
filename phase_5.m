function S5=phase_5(t,s)
% %Motion equations phase 1: First free fly
global m_slv J_se J_re J_ce J_p R_s R_r R_p R_c r_g P1 P2...
    c_slv theta_g N N_h k ig kesi K_con D_con mu_con Jx1 Jx2...
    Jx3 Jx4 Jx Lrs1 Lrs2 Lrs3 Lrs4 Lcp1 Lcp2 Lcp3 Lcp4 K_d...
    data_save
%% --------------------------Variables----------------------------
% S1=[
% 1)x_slv,
% 2)theta_ring,
% 3)theta_sun,
% 4)dx_slv,
% 5)dtheta_ring,
% 6)dtheta_sun,
%% --------------------------Dynamics-----------------------------

       A = [-c_slv/m_slv 0 0 ;0 0 0 ; 0 0 0];
       A_5= [zeros(3,3) eye(3);
               zeros(3,3) A]; 
        B1_5 = [-1/m_slv 0;0 0;0 1/(J_se+J_re+J_ce)];
        Bu_5 = [zeros(3,2);B1_5];
        B2_5 = [0 0;0 0;0 -1/(J_se+J_re+J_ce)];
        Bw_5 = [zeros(3,2);B2_5];

%% dx = A_3*x + Bu_3*[F_slv;T_s]+Bw_3*[T_r;T_c];

S5 = A_5*s + Bu_5*[0;0]+Bw_5*[0;0];
%%
% data_save=[data_save;s,t];
end