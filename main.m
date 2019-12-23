%%Wireless Communication Networks Hackathon 2
close all
clear variables
clc

%% Phase 1: Ray Matricec - ABCD Matrix
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
%% Phase 3: Beam Expander and Concentrator
%1

%% Phase 4: The Communication Link


%% Phase 5: The Communication Link


%% Phase 6: Filter


%% Phase 7: Link Budget


%% Phase 8: Detector


%% Phase 9: Ideal Receiver
