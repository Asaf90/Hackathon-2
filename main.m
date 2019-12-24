%%Wireless Communication Networks Hackathon 2
close all
clear variables
clc

%% Phase 1: Ray Matricec - ABCD Matrix
disp('-----Phase 1-----')
%1
f1 = 50;
f2 = 25;
dist = 40;
beam_in = [10 0]';
M1 = [1 0; -1/f1 1];
M2 = [1 dist; 0 1];
M3 = [1 0; -1/f2 1];

M = M3*M2*M1 %#ok<NOPTS>
beam_out = M*beam_in %#ok<NOPTS>


%2
syms d
syms n2

M4 = [1 0; 0 1/n2];     %Assuming n1 = 1
M5 = [1 d; 0 1];
M6 = [1 0; 0 n2];
M_slab = M6*M5*M4 %#ok<NOPTS>

%3
f = f1;                 %We chose f = f1 = 50mm;
delta_f = f * 1e-4;
dist = 2*f + delta_f;
M_lens = [1 0; -1/f 1];
M_dist = [1 dist; 0 1];
M_comp = M_lens * M_dist;
n = 2:50;
for i = n
    M_n = M_lens * M_dist * M_comp^(i-2) * M_lens;
    beam_out(:,i-1) = M_n * beam_in;
end

figure()
plot(n, beam_out(1,:), '*')
grid on
title('Exit Beam Height as a Function of N')
xlabel('Number of lenses')
ylabel('Exit Beam Height [mm]')

figure()
plot(n, beam_out(2,:), '*')
grid on
title('Exit Beam Angle as a Function of N')
xlabel('Number of lenses')
ylabel('Exit Beam Angle [deg]')
ylim([-0.3 0.3])

%% Phase 2: Spatial Propagation
disp('-----Phase 2-----')
%1
r = -5:0.01:5;
I0 = 1;             %Assume peak intensity is 1.
w = 20e-1;          %Assume w = 20cm
I = I0*exp(-2*r.^2 / w.^2);
figure()
plot(r,I)
title('Gaussian Profile Beam')
xlabel('Radial distance [m]')
ylabel('Intensity [W/m^2]')

%2
w0 = 20e-1;             %Assume w0 = 20cm
lambda = 1064e-9;       %Assume wavelength is 1064nm
z = 0:0.1:1e3;
w_z = w0 * sqrt(1+((lambda * z)/(pi * w0^2)).^2);

figure()
plot(z, w_z)
title('Propagation Distance')
xlabel('Distance z [m]')
ylabel('Spot Size [m]')

%3
w = 1;
lambda = 1e-9;

FWHM = sqrt(-0.5*log(0.5))  %#ok<NOPTS>
D_135 = 1                   %#ok<NOPTS>
z0 = w^2 * pi * lambda      %#ok<NASGU,NOPTS>

%4
w0 = 1e-2;
n = 1;
theta_beam = lambda / (pi * w * n);
z0 = w0 / theta_beam;

z = 1e3;
w = w0 * z/z0               %#ok<NOPTS>

%5
theta = -pi/2:pi/200:pi/2; 
I_lamb = cos(abs(theta));
P_lamb = cumtrapz(I_lamb);
figure()
plot(theta,P_lamb);
title('Accumulated Power of Lambertian Light Srouce')
xlabel('\theta')
ylabel('Accumulated Power')
%% Phase 3: Beam Expander and Concentrator
disp('-----Phase 3-----')
%1
Din = 1e-2;
Dout = 15e-2;
f_num = 1;
f1 = Din * f_num            %#ok<NASGU,NOPTS>
f2 = Dout * f_num           %#ok<NASGU,NOPTS>

%2
f1 = 10e-2;
f2 = 100e-2;
lambda = 1.55e-6;

theta_in = lambda / Din;
theta_out = (f1/f2) * theta_in    %#ok<NASGU,NOPTS>

%Assuming Gaussian Beam:
w0 = lambda / (pi* theta_out);
z0 = w0 / theta_out;
z = 1e3;

w = w0 * z/z0                       %#ok<NASGU,NOPTS>
%% Phase 4: The Communication Link


%% Phase 5: The Communication Link
disp('-----Phase 5-----')
%1
r = 0:1e-3:10;10e-6;
z = 0:0.01:100;
con = 2e5 * exp(-r/10^-6);
cross = 10^-6 * exp(-r/10^-6);
gamma = trapz(con.*cross);

tau = exp(-gamma*z);
figure()
plot(z,tau)
title('Attenuation of Beam Through Aerosols')
xlabel('Range [m]')
ylabel('Attenuation')

%2
P = 1013.25;
temp = 300;
n0 = 1 + 7.8e-5 * P/temp            %#ok<NASGU,NOPTS>

%3          TODO
I0 = 1;
mx = 1;
sigmax = 1;
th = 1;

%4
TSI = 1.365e3                        %#ok<NASGU,NOPTS>
%From Wikipeida
%% Phase 6: Filter
disp('-----Phase 6-----')
n = 1.5;
d = 1e-6;

theta_i = 0:pi/200:pi/2;

theta_t = asin(sin(theta_i)/n);

lambda = 2*d*n*cos(theta_t);
figure()
plot(theta_i, lambda)
title('Transmission of Interference Filter')
xlabel('\theta_i [rad]')
ylabel('Wavelength [m]')
%% Phase 7: Link Budget
disp('-----Phase 7-----')
%TODO
range = 10e3;           %range = 10km
tx_db = 30;             %tx power = 30dB
apper = 0.01;           %tx and rx apperture = 0.01
opt_amp = 30;           %optical amplifier = 30dB
NF = 5;                 %NF = 5dB
lambda = 1.55e-6;       %wavelength = 1.55um
qe = 0.8;               %PIN diode quantum efficiency
tx_point_err = 10e-3;   %tx pointing error = 10mrad


S = 10;                 %snowfall rate 10 mm/h
a = 0.0001023*lambda*1e9 + 3.7855466;   %a for wet snow
b = 0.72;                               %b for wet snow
gamma_snow = a*S^b;     %wet snow attenuation

P_r = 0 + 2*opt_amp - 20*log10(apper) - 10*log10((lambda/4*pi*range)^2) - NF;



%% Phase 8: Detector
disp('-----Phase 8-----')
%1
disp('Q1: 12 bits')

k_eff = 0.01;
sigma2_th = 10e-14;
R = 1;
B = 1e9;
P_R = 1e-6;

M = 0:0.0001:0.01;
SNR = (M*R*P_R).^2 ./ (2*k_eff*R*P_R*B*M.^2.5 + sigma2_th);
figure()
plot(M, SNR);
title('SNR as a function of Gain')
xlabel('Gain')
ylabel('SNR')
[max_SNR, max_IDX] = max(SNR);
max_M = M(max_IDX)              %#ok<NASGU,NOPTS>


%% Phase 9: Ideal Receiver
disp('-----Phase 9-----')
th_1 = 0:1:50;
% x = 1:length(th_1);
% th_2 = 0:1:50;
% y = 1:length(th_1);
% th_3 = 0:1:50;
% z = 1:length(th_1);
lambdaVec = [1 3 5 7];
Q_gen=@(x,lambda) ((lambda.^x).*exp(-lambda)./factorial(x));
min_location = [0,0,0];
min=1;
for x=1:length(th_1)
    for y=x:length(th_1)
        for z=y:length(th_1)
            seg0 = 0:x; 
            seg1 = x:y;
            seg2 = y:z;
            seg3 = z:length(th_1);
            
            correct1 = Q_gen(seg0,lambdaVec(1));
            val1 = sum(correct1);
            correct2 = Q_gen(seg1,lambdaVec(2));
            val2 = sum(correct2);
            correct3 = Q_gen(seg2,lambdaVec(3));
            val3 = sum(correct3);
            correct4 = Q_gen(seg3,lambdaVec(4));
            val4 = sum(correct4);
            val = 1-0.25*(val1+val2+val3+val4);
            if (val<min)
                min = val;
                min_location=[x,y,z];
            end
        end
    end
end
val                         %#ok<NASGU,NOPTS>
[x,y,z]                     %#ok<NASGU,NOPTS>