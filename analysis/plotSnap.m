%close all
clear all
clc

%% Generate initial conditions


%% parameters
N       = 64;                                      % resolution
beta    = 0.06;  0.06;0.03;                         % gravity
Lbox    = beta^-1;

b = 1;
G = 1;

snap = 1;10;

simDir = '../output/';

%% read snapshot
snapFile = [simDir 'N' num2str(N) 'box' num2str(Lbox) 'B' num2str(beta) 'b' num2str(b) 'G' num2str(G) '/snap_' sprintf('%03d',snap) '.h5' ];

a = h5read(snapFile,'/a')

psi = h5read(snapFile,'/psiRe') + 1.i * h5read(snapFile,'/psiIm');



%% plot snapshot

figure;
imagesc(log10(mean(abs(psi).^2,3)))
axis square
colorbar
caxis([-.01 .01])
%caxis([-1 1])




