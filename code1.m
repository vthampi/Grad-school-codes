function dydt = newmodel(t,x,Kappa,w)

%x(1) range - [7500 , 15000] kg C /ha
%x(2) range - [1000 , 7500]  kg C /ha
%x(3) range - dunno kg N /ha

dx = [0; 0; 0; 0];
w_n = 39.6;        %range [30,40] kg C /kg N]

w_e = 27.1;        %range [20,30] kg C /kg N]

v_n = 0.0211;        %range [0.02,0.5] ha /kgC /yr

v_e = 0.0877;        %range [0.02,0.5] ha /kgC /yr

mu_n = 0.14;       %range [0.0,0.2] /yr

mu_e = 0.81;       %range [0.0,1.0] /yr

m_n = 0.17;        %range [0.0,0.5] /yr

m_e = 0.32;        %range [0.0,0.5] /yr

delta_n = 0.06;    %range [0.003,0.7]

delta_e = 0.06;    %range [0.003,0.7]

k = 0.4;          %range [0.0,0.1] /yr

I = 2;             %range [0,15] kg N/ ha

h = 3890;          %range [500, 10000] kg C/ha     [only used for coexistence]

p = 3;          %range {1,2,3}     [only used for coexistence]

alpha_ne = 19.00;   %range [0.0,100.0]

alpha_en = 0.05;   %range [0.0,100.0]

%g_n = (I * v_n * w_n)/((v_n * delta_n * x(1) + k) * m_n) - mu_n/m_n; %only for dominant

%g_e = (I * v_e * w_e)/((v_e * delta_e * x(2) + k) * m_e) - mu_e/m_e; %only for dominant

g_n = x(1)^p/(h^p + x(1)^p + x(2)^p);  %[only used for coexistence]

g_e = x(2)^p/(h^p + x(1)^p + x(2)^p);  %[only used for coexistence]

f_n = g_n + alpha_ne * g_e;

f_e = alpha_en * g_n + g_e;

rho = 3; % 3 is baseline. efficacy of activists in reducing the Nitro input


%Kappa = 1.5; %rescaled sampling rate 1.5 is baseline

%w = 0.1; % baseline is  0.1. cost of supporting runoff reduction program, relative to utility penalty of tolerating invasive species

epsilon = 0.001; % 0.001 is baseline.

dx(1) = x(1) * (w_n * v_n * x(3) - mu_n - m_n * f_n);
dx(2) = x(2) * (w_e * v_e * x(3) - mu_e - m_e * f_e);
dx(3) = I - k * x(3) - x(1) * (v_n * x(3) - (mu_n + m_n * f_n) * (1 - delta_n)/w_n) - x(2) * (v_e * x(3) - ( mu_e + m_e *f_e) * ( 1 - delta_e)/w_e) + rho * (1 - x(4));
dx(4) = Kappa * x(4) * (1 - x(4)) * (-w + x(2)/(x(1) + x(2)))  + epsilon * (1 - 2 * x(4));
% if x(4) < 0.05
%     dx(4) = 0.05 + x(4);
% elseif x(4) > 0.95
%     dx(4) = 0.95 - x(4);
% else
%     dx(4) = Kappa * x(4) * (1 - x(4)) * (-w + g_e);
% end
dydt=dx;
