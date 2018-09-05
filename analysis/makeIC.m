close all
clear all
clc

%% Generate initial conditions


%% parameters
N       = 64;128;                      % resolution
beta    = 0.06;  0.06;0.03;                         % gravity
Lbox    = beta^-1;

a = 1;

simDir = '../output/';


%% setup IC
rng(42);

kmax = sqrt(2)/2; % /sqrt(a) ?? fastest growing


% spacing & grid
dx = Lbox / N;

% fourier space variables
klin = ((-N/2:N/2-1)') * (2*pi/Lbox);
kmaxRes = max(klin);
kminRes = min(klin(klin>0));
log10([kminRes kmax kmaxRes])


[kx, ky, kz] = meshgrid(klin, klin, klin);
kSq = fftshift(kx.^2 + ky.^2 + kz.^2);
clear klin;
clear kx;
clear ky;
clear kz;
[~,sid] = sort(kSq(:));


%% IC -- construct in fourier space accoridng to Eq (27) of our paper
variance = 0.5;

rng(42);
psi = zeros(size(kSq));
psi(sid)= exp(1.i*2*pi*rand(size(kSq(:))));  % initialize with random phases
rng(42);
psi(sid)= psi(sid).* randn(size(kSq(:)));
psi = psi .* ( sqrt(variance).* exp( -kSq/(2*(kmax)^2)) ) * (sqrt(2*pi)/Lbox*N)^3;
clear sid


psi = ifftn(psi);

psi = (1 + psi);

% check average density
rhobar = sum(abs(psi(:)).^2) * dx^3 / Lbox^3
assert(abs(rhobar-1) < 1e-2, 'more than 1 percent fluctuations in IC!');

% domain
xlin = Lbox * linspace(0,1,N+1)';
xlin = xlin(1:N) + dx/2;


%% Save IC file

snapFile = [simDir 'ic' num2str(N) 'box' num2str(Lbox) 'B' num2str(beta)  '.h5'];
if exist(snapFile,'file')
    delete(snapFile)
end
h5create(snapFile,'/a',size(a));
h5create(snapFile,'/psiRe',size(psi));
h5create(snapFile,'/psiIm',size(psi));
h5write(snapFile,'/a',a);
h5write(snapFile,'/psiRe',real(psi));
h5write(snapFile,'/psiIm',imag(psi));



figure;
imagesc(mean(log10(abs(psi).^2),3))
caxis([-.01 .01])


