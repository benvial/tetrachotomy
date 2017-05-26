
clear all
close all

% initialisation of variables
global om pol n nf dist h fg lambg eps2 theta_inc res op or einf
fg=[];lambg=[];

c=299792458; % light speed

% definition of rect0 (the rectangular restriction of complex plane where to search for the poles):
% M0 (X0+i*Y0) : lower left corner
% M1 (X1+i*Y1) : upper right corner

X0=3e14;X1=8e14;
Y0=-1e14;Y1=0e14;
rect0.x0=X0;rect0.y0=Y0;
rect0.x1=X1;rect0.y1=Y1;

theta0=0; % "angle" (there is actually no angle because no incident wave,
% it is just a dummy variable to define alpha0)
alpha0=1e14/c*sin(pi/180*theta0); % quasiperiodicity coefficient (transverse wavenumber)

% number of points for plots
nx=250;ny=250;
omegare=linspace(X0,X1,nx);
omegaim=linspace(Y0,Y1,ny);
lambda=2*pi*c./omegare;

% plot figures or not
figs=1;
%
[X,Y]=meshgrid(omegare,omegaim);
O_c=X+1i*Y; %complex omega

% permitivity of dispersive material

op=4e14;
or=5e14-1e13*1i;
einf=7;
ep=einf-op^2./(omegare.^2-or^2);

figure
plot(omegare,real(ep))
hold on
plot(omegare,imag(ep),'r--')
legend('Re','Im');xlabel('Re \omega');ylabel('\epsilon')
ep_c=einf-op^2./(O_c.^2-or^2);
figure;fh=gca;pcolor(X,Y,log10(abs(ep_c)));shading interp;hold on;colorbar
xlabel('Re \omega');ylabel('Im \omega');title('log|\epsilon|')

% f: function we search for poles (function handle)
%(here function trans: complex transmission coefficient of a multilayer stack)

f=inline(strcat('trans(Z,',num2str(alpha0),',0)'),'Z');
% plots of f on the real axis
figure, hold on
plot(omegare,real(f(omegare)))
plot(omegare,imag(f(omegare)),'r--')
legend('Re','Im');xlabel('Re \omega');ylabel('f')

figure
plot(omegare,abs(f(omegare)).^2,'r')
xlabel('Re \omega');ylabel('|f|^2')

% plot of ln(|f|) in complex plane
t_c=f(O_c);
figure;fh=gca;pcolor(X,Y,log10(abs(t_c)));shading interp;hold on;colorbar
xlabel('Re \omega');ylabel('Im \omega');title('log|f|')
set(gca,'Clim',[0 2])

%%

%%


% tolerance on I0
tolre0=(X1-X0)/10000;
tolim0=(Y1-Y0)/10000;
% tolerance on I1
tolre1=tolre0^2;
tolim1=tolre0^2;

% relative tolerance on ratios R21 et R10
tolreR=1e-6;
tolimR=1e-6;

tol_conv=1e-3;  %relatic tolerance for the convergence of contour integrals In

pol=[]; %initialisation : the pol variable contains the results of found poles
res=[]; %initialisation : the res variable contains the results of found residues

n=0;    % number of calls of the pole_2D function
nf=0;   % number of calls of the function f

figure
trace=1; % plot rectangular sections while searching

% tetrachotomy
tic
pole_2D(f,rect0,tolre0,tolim0,tolre1,tolim1,tolreR,tolimR,tol_conv,trace) ;
pol

res

n
nf
fprintf('%i pole(s) found\n',length(pol));
t_calcul=toc


plot(real(pol),imag(pol),'MarkerSize',8,'Marker','x','LineWidth',1,'LineStyle','none');
hold on
xlabel('Re \omega');ylabel('Im \omega');

plot(fh,real(pol),imag(pol),'MarkerSize',8,'Marker','x','LineWidth',1,'LineStyle','none');
hold on
