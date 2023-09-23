clear;
clc;
s=tf('s');
states = ["u", "v", "w", "p", "q", "r", "phi", "theta", "psi", "x", "y", "z"];
inputs = ["phi_dd", "theta_dd", "psi_dd", "a_x", "a_y", "a_z"];
outputs = ["phi", "theta", "psi", "x", "y", "z"];
z = zeros(3);
one3 = eye(3);

w1 = 2.5/(s+0.001);
w2  = tf([105 63], [1 16 64]);

W1 = w1*eye(6);
W2 = w2*eye(6);

% Jacobi Linear approximation

% Ac = [z one3 z z;
%       z z z z;
%       z z z one3;
%       z z z z];
% 
% Bc = [z z;
%      one3 z;
%      z z;
%      z one3];
% 
% Cc = [one3 z z z;
%      z z one3 z];
% 
% Dc = 0;


% sys1 = ss(Ac,Bc,Cc,Dc, "statename", states, "inputname", inputs, "outputname", outputs);
% [K1,CL1,GAM1,INFO1] = ncfsyn(sys1,W1,W2);


% Robust feedback linearization

m = 0.5;
g = 9.81;
l = 0.23;
Ixx = 6.8*10^-3;
Iyy = 5.3*10^-3;
Izz = 1.7*10^-3;
Iyz = 3.1*10^-4;
kf = 1.97*10^-5;
kd = 2.88*10^-7;

Ar = [0 0 0 0 0 0 0 g 0 0 0 0;
      0 0 0 0 0 0 -g 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0;
      z z z z;
      z one3 z z;
      one3 z z z];

Br_r1 = [0 (-0.5*sqrt(3)*kf/m) (0.5*sqrt(3)*kf/m) 0 0 0];
Br_r2 = [(kf/m) (-0.5*kf/m) (-0.5*kf/m) 0 0 0];
Br_r3 = [0 0 0 (kf/m) (kf/m) (kf/m)];
Br_r4 = [0 (0.5*sqrt(3)*kd/Ixx) (-0.5*sqrt(3)*kd/Ixx) 0 (0.5*sqrt(3)*l*kf/Ixx) (-0.5*sqrt(3)*l*kf/Ixx)];
Br_r5 = [(Iyz*kf*l - Izz*kd) (Iyz*kf*l + 0.5*Izz*kd) (Iyz*kf*l + 0.5*Izz*kd) -(Izz*kf*l + Iyz*kd) (0.5*Izz*kf*l - Iyz*kd) (0.5*Izz*kf*l - Iyz*kd)]/(Iyy*Izz - Iyz^2);
Br_r6 = [(Iyy*kf*l - Iyz*kd) (Iyy*kf*l + 0.5*Iyz*kd) (Iyy*kf*l + 0.5*Iyz*kd) -(Iyz*kf*l + Iyy*kd) (0.5*Iyz*kf*l - Iyy*kd) (0.5*Iyz*kf*l - Iyy*kd)]/(Iyy*Izz - Iyz^2);
Br_low = [z z;
          z z];
Br = [Br_r1; Br_r2; Br_r3; Br_r4; Br_r5; Br_r6; Br_low];

Cr = [z z one3 z;
      z z z one3];

Dr = 0;

sys2 = ss(Ar,Br,Cr,Dr, "statename", states, "inputname", inputs, "outputname", outputs);
[K2,CL2,GAM2,INFO2] = ncfsyn(sys2,W1,W2);




