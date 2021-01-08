% Matlab code to determine the maximum force on the human after collision

clc
clear all
close all

%constants 
rho = 68; % kg/m^3 density of Aluminium-SiC foam, shield material
rho_foam = 9; % density of foam in kg/m^3
t = 0.03; % thickness of shells in m
v0 = 140/3.6; % initial speed in m/s
mh = 100; % weight of human in kg
Telem   = 3.14E3; % Toughness of insulation board in J/m^2, breaking elements
Tshield = (133+234)/2; % Toughness of Aluminium-SiC foam, shield material
g = 9.81; % gravitational acceleration in m/s^2
Cfoam = (1E5 + 2e5)/2; % Compressive modulus of the foam in Pa
rho_elem = (160+480)/2; % density of insulation board in kg/m^3, breaking elements
epsilon = 0.87; % Maximum compression of foam
rd = 0.2; % size reduction of consecutive shells
W = 2.85; % width of the first shell
Fs = (1.5E6 + 5.5E6)/2; % Flexural strength of the breaking element

% volume of shells
V1 = W*2*t + 1*2*t + 1*2*t;
V2 = (W-rd)*(2-rd)*t + (1+0.1)*(2-rd)*t + (1+0.1)*(2-rd)*t;
V3 = (W-rd*2)*(2-rd*2)*t + (1+0.1*2)*(2-rd*2)*t + (1+0.1*2)*(2-rd*2)*t;
V4 = (W-rd*3)*(2-rd*3)*t + (1+0.1*3)*(2-rd*3)*t + (1+0.1*3)*(2-rd*3)*t;
V5 = (W-rd*4)*(2-rd*4)*t + (1+0.1*4)*(2-rd*4)*t + (1+0.1*4)*(2-rd*4)*t;
V6 = (W-rd*5)*(2-rd*5)*t + (1+0.1*5)*(2-rd*5)*t + (1+0.1*5)*(2-rd*5)*t;
V7 = (W-rd*6)*(2-rd*6)*t + (1+0.1*6)*(2-rd*6)*t + (1+0.1*6)*(2-rd*6)*t;
V8 = (W-rd*7)*(2-rd*7)*t + (1+0.1*7)*(2-rd*7)*t + (1+0.1*7)*(2-rd*7)*t;
V9 = (W-rd*8)*(2-rd*8)*t + (1+0.1*8)*(2-rd*8)*t + (1+0.1*8)*(2-rd*8)*t;

%thickness of the foam
t_foam = [2.8; 2.2; 1.1; 0.9; 0.7; 0.4;0.3;0.2;0.1]; %thickness of the foam in m

% volume of the foam
V1_f = W*2*t_foam(1);
V2_f = (W-rd)*(2-rd)*t_foam(2);
V3_f = (W-rd*2)*(2-rd*2)*t_foam(3);
V4_f = (W-rd*3)*(2-rd*3)*t_foam(4);
V5_f = (W-rd*4)*(2-rd*4)*t_foam(5);
V6_f = (W-rd*5)*(2-rd*5)*t_foam(6);
V7_f = (W-rd*6)*(2-rd*6)*t_foam(7);
V8_f = (W-rd*7)*(2-rd*7)*t_foam(8);
V9_f = (W-rd*8)*(2-rd*8)*t_foam(9);

% Area of foam
Afoam = [W*2;(W-rd)*(2-rd);(W-rd*2)*(2-rd*2);(W-rd*3)*(2-rd*3);(W-rd*4)*(2-rd*4);(W-rd*5)*(2-rd*5)
    (W-rd*6)*(2-rd*6);(W-rd*7)*(2-rd*7);(W-rd*8)*(2-rd*8)];

% mass of shells
m_s1 = rho*V1;
m_s2 = rho*V2;
m_s3 = rho*V3;
m_s4 = rho*V4;
m_s5 = rho*V5;
m_s6 = rho*V6;
m_s7 = rho*V7;
m_s8 = rho*V8;
m_s9 = rho*V9;

% mass of foam
mf_1 = rho_foam*V1_f;
mf_2 = rho_foam*V2_f;
mf_3 = rho_foam*V3_f;
mf_4 = rho_foam*V4_f;
mf_5 = rho_foam*V5_f;
mf_6 = rho_foam*V6_f;
mf_7 = rho_foam*V7_f;
mf_8 = rho_foam*V8_f;
mf_9 = rho_foam*V9_f;
m_s = [m_s1+mf_1;m_s2+mf_2;m_s3+mf_3;m_s4+mf_4;m_s5+mf_5;m_s6+mf_6;m_s7+mf_7;m_s8+mf_8;m_s9+mf_9];


%m_s = [m_s1;m_s2;m_s3;m_s4;m_s5;m_s6;m_s7;m_s8];
% volume of breaking elements 
Delem = 0.02; % diameter of breaking element in meters
A_elem = (Delem/2)^2*pi;
V_elem = A_elem*0.3; %

% mass of breaking elements
m_elem = rho_elem*V_elem;

% Breaking energy of element
Ebreak = Telem*A_elem;


% Energy absorbed by the foam
  k1 = Cfoam*V1_f/(t_foam(1)^2);
% k2 = Cfoam*V2_f/(t_foam(2)^2);
% k3 = Cfoam*V3_f/(t_foam(3)^2);
% k4 = Cfoam*V4_f/(t_foam(4)^2);
% k5 = Cfoam*V5_f/(t_foam(5)^2);
% k6 = Cfoam*V6_f/(t_foam(6)^2);
% k7 = Cfoam*V7_f/(t_foam(7)^2);
% k8 = Cfoam*V8_f/(t_foam(8)^2);
% k9 = Cfoam*V9_f/(t_foam(9)^2);
% kfoam = [k1;k2;k3;k4;k5;k6;k7;k8;k9];
% 
% Efoam = [0.5*k1*(t_foam*epsilon)^2; 0.5*k2*(t_foam*epsilon)^2;0.5*k3*(t_foam*epsilon)^2;0.5*k4*(t_foam*epsilon)^2
%      0.5*k4*(t_foam*epsilon)^2;0.5*k5*(t_foam*epsilon)^2;0.5*k6*(t_foam*epsilon)^2;0.5*k7*(t_foam*epsilon)^2;
%      0.5*k8*(t_foam*epsilon)^2]

%calculate wind force on first shell
rho_air = 1.225; %density of air in kg/m^3;
A_shell1 = 2.85*2; % surface area of first shell;
Cd = 0.3; % drag coefficient of first shell;
F_w = 0.5*Cd*A_shell1*rho_air*v0^2; % drag force on first shell

% calculate remaining compression length of first foam layer after drag
Lw = F_w/k1; 
t_foam(1) = t_foam(1) - Lw;

% ordering of elements per shield
N = [4 4 4 4 5 6 7 8 9];
Ebreak_vec = Ebreak*N;
% calculate critical force
Fcrit = 40*mh*g;

% creating vectors to store data
Fmax = []; 
vh = [];
E_foam = [];

% Force needed to break through layers
F_break = Fs*A_elem;


% calculate maximum force after each stage
for i = 1:9
    ufoam =  N(i)*F_break/(Cfoam*Afoam(i)/t_foam(i));
% ufoam =  N(i)*F_break/(kfoam(i));
    Efoam = 0.5*(Cfoam*Afoam(i)/t_foam(i))*ufoam^2
    v = sqrt((mh*v0^2 - 2*N(i)*Ebreak - 2*Efoam)/(mh+sum(m_s(1:i)))); %calculate new velocity
    vmean = (v0 + v)/2; % mean velocity
    t_impact = t_foam(i)*epsilon/vmean; %impact time
    Favg = mh*(v0-v)/t_impact; % average force
    
    if(v<=0)
        v = 0;
        Favg = 0;
    end
    vh = [vh;v];
    Fmax = [Fmax;Favg*2]; 
    v0 = v;
    E_foam = [E_foam;Efoam]
end

% plotting
figure
plot(Fmax)
hold on
plot(Fcrit*ones(9),'r')
ylabel('Force on human in N')
xlabel('# shell')
legend('Inertia force','Critcal force')
title('Force on human')

v0 = 140/3.6;
vh = [v0;vh];
figure 
plot(0:9,vh)
xlabel('# shell')
ylabel('velocity in m/s')
title('Velocity of human after each shell')
