clc
clear all
close all

a_12 = 0.1;    %in
a_34 = 0.15;   %in
St_13 = 0.0034;   %in/tooth
St_24 = 0.0057;   %in/tooth
R=(3/4)/2;  %cutter radius

Le_datasheet = 0.393701;   %cutting edge effective length (in)
Re_datasheet = 0.031496;   %corner radius (in)
Bs_datasheet = 0.047244;   %wiper edge length (in)

Tc_1 = 19*0.5*0.25*12;   %max torque for cut 1 (in-lbs)
Tc_2 = 19*0.5*0.31*12;   %max torque for cut 2 (in-lbs)
Tc_3 = 19*0.5*0.33*12;   %max torque for cut 3 (in-lbs)
Tc_4 = 19*0.5*0.42*12;   %max torque for cut 4 (in-lbs)

Le_12 = 2*Bs_datasheet+(a_12-Re_datasheet)+((2*Re_datasheet*pi())/4);   %in
Le_34 = 2*Bs_datasheet+(a_34-Re_datasheet)+((2*Re_datasheet*pi())/4);  %in

syms k1 k2
Torque_1=(2*R*a_12*St_13*k1)/pi() + (R*k2*Le_12);
Torque_2=(2*R*a_12*St_24*k1)/pi() + (R*k2*Le_12);
Torque_3=(2*R*a_34*St_13*k1)/pi() + (R*k2*Le_34);
Torque_4=(2*R*a_34*St_24*k1)/pi() + (R*k2*Le_34);

eqn1 = Tc_1==vpa(Torque_1);    %solving first system of linear equation
eqn2 = Tc_2==vpa(Torque_2);
[A,B]=equationsToMatrix([eqn1 eqn2],[k1 k2]);
X1=linsolve(A,B);

eqn3 = Tc_3==vpa(Torque_3);    %solving second system of linear equations
eqn4 = Tc_4==vpa(Torque_4);
[A,B]=equationsToMatrix([eqn3 eqn4],[k1 k2]);
X2=linsolve(A,B);

[k1;k2]==(X1+X2)/2