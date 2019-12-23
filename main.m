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


%% Phase 3: Beam Expander and Concentrator


%% Phase 4: The Communication Link


%% Phase 5: The Communication Link


%% Phase 6: Filter


%% Phase 7: Link Budget


%% Phase 8: Detector


%% Phase 9: Ideal Receiver
