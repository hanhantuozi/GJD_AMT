function S4 = phase_4(t,s)
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
x_slv = s(1);
theta_slv = s(2);
theta_sun = s(3);

dx_slv = s(4);
omega_slv = s(5);
omega_sun = s(6);


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
N_cons =0;
f_cons =0;
% K_cons =5e5; 
% D_cons = 2.5e3;
K_cons =3e6; 
D_cons =3e4;

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

 F_slv =250;

%% collision heppen dynamics
colis = 0;
delta = 0;
ddelta = 0;
N_cons = 0;
f_cons = 0;
mu_cons = 0.002;
% 计算单个齿间碰撞

%  % ----------- 错误思路：通过相对位移判断
% if x_slv>P2
%         delta =  R_r*mod(theta_slv -theta_sun ,2*pi/N_h)-R_r*2*pi/N_h/2; % 此处dtheta乘R_r得到相对位移，用于计算碰撞线性力
%         ddelta =  R_r*(s(5)-s(6));                                         % 此处domega乘R_r得到相对速度，用于计算碰撞阻尼力
%         if  delta<0
%             colis = 1 ;% 如果接合套速度大于太阳轮速度那么一定在接合套转动方向下齿面，太阳轮上齿面，进行碰撞
%             % 太阳轮本来规定向下为正，左侧顺时针，此时colis==1表示太阳轮加速，接合套降速
%         elseif  delta>0
%             colis = -1; % 如果接合套速度大于太阳轮速度那么一定在接合套转动方向上齿面，太阳轮下齿面，进行碰撞
%             % 太阳轮本来规定向下为正，左侧顺时针，此时colis==-1表示太阳轮降速，接合套加速
%         else 
%             colis = 0;
%         end
% end  

% ------- 错误思路： 通过相对速度判断是否产生碰撞， 这样会一直产生跳变碰撞力
if x_slv>P2
        delta =  R_r*mod(theta_slv -theta_sun ,2*pi/N_h); % 此处dtheta乘R_r得到相对位移，用于计算碰撞线性力
        ddelta =  R_r*(s(5)-s(6));                                         % 此处domega乘R_r得到相对速度，用于计算碰撞阻尼力
        if delta< R_r*pi/N_h&&delta>0 %s(5)-s(6)>  0 &&
            colis = 1 ;% 如果接合套速度大于太阳轮速度那么一定在接合套转动方向下齿面，太阳轮上齿面，进行碰撞
            % 太阳轮本来规定向下为正，左侧顺时针，此时colis==1表示太阳轮加速，接合套降速
                % if (s(5)-s(6)>0

        elseif delta>R_r*pi/N_h  %s(5)-s(6)<0&&
            colis = -1; % 如果接合套速度大于太阳轮速度那么一定在接合套转动方向上齿面，太阳轮下齿面，进行碰撞
            % 太阳轮本来规定向下为正，左侧顺时针，此时colis==-1表示太阳轮降速，接合套加速
        else 
            colis = 0;
        end
end  

    N_cons = (K_cons*delta+D_cons*ddelta);%*( abs(ddelta)>1e-8 )
    f_cons = - (abs(mu_cons*N_cons)*(s(4)>0)); % 只要速度向内侧则方向向外，当不存在速度时静摩擦力先不考虑，
% 接合套收到的总的齿的摩擦力 需要乘N_h
    F_slv = F_slv+N_h*(f_cons);
    % 转矩的方向，先按接合套往负方向减速定
    T_r =   -colis*N_cons*N_h*R_s;         % slv == ring 认为向下为＋，认为slv > sun , 符号仅表示colision =1时T_r 是使接合套减速的，反之加速
    T_s =    colis*N_cons*N_h*R_s;         
    T_c = 0;

 Bu_3 = B_u;  
 Bw_3 = B_w;

S4 = A*s + Bu_3*[F_slv;T_s]+Bw_3*[T_r;T_c]+[0;0;0;-10*0.02*(s(4)>0)+10*0.02*(s(4)<0);0;0];

% data_save=[data_save;t,delta,F_slv,N_h*(f_con*cos(theta_g)+F_con*sin(theta_g)),delta,colis];
data_save=[data_save;t,delta,ddelta,F_slv,colis,N_cons*N_h,f_cons,T_r,T_s];
end