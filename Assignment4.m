clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','docked')
VR3 = linspace(0.1,10,10);
nx =  20;
ny =  10;
Vo = 0.1;
f = zeros(1,nx*ny);
Gmap = sparse(ny*nx);
M = zeros(nx,ny);
l = 200e-9;
w = 100e-9;
sigma = ones(nx,ny);
T = 300;
q = 1.602e-19;
n_elec = 1e19; %m^-2 = 1e15 cm^-2   (per cm squared)
kb = 1.38064852*10^-23;
kt = 1.38064852*10^-23*T;
m = (9.10938356*10^-31)*0.26;
vth = sqrt(2*kt/m);
numberofelectrons = 10000;
xp = zeros(numberofelectrons,1);
X = rand(numberofelectrons,1)*l;
Y = rand(numberofelectrons,1)*w;
yp = zeros(numberofelectrons,1);
Nt = 500;
tmn = 0.2e-12;
dt = (1/200)*(l/vth);
x = linspace(0,l,nx);
y = linspace(0,w,ny);
dx = l/(nx-1);
dy = w/(ny-1);
%boxes for particle trajectory plot
boxes = {};
boxes{1}.X = [0.8 0.4]*1e-7;
boxes{1}.Y = [0.6 0.4]*1e-7;
boxes{2}.X = [0.8 0.4]*1e-7;
boxes{2}.Y = [0 0.4]*1e-7;

%loop to make sure no particles are generated inside the boxes
it = (X < 1.2e-7) & (X > 0.8e-7) & (Y > 0.6e-7);
ib = (X < 1.2e-7) & (X > 0.8e-7) & (Y < 0.4e-7);
iXY = it | ib;
while(sum(iXY) > 0)
    X(iXY) = rand(1,sum(iXY))*l;
    Y(iXY) = rand(1,sum(iXY))*w;
    it = X < 1.2e-7 & X > 0.8e-7 & Y > 0.6e-7;
    ib = X < 1.2e-7 & X > 0.8e-7 & Y < 0.4e-7;
    iXY = it | ib;
end



%Mapping for conductivity
 for z = 1:length(VR3)
for i = 1:nx
   for j = 1:ny
    if x(i)> l*0.3 & x(i)< l*0.65 & (y(j)< w*0.3 | y(j)> w*0.6)
        sigma(i,j) = 10e-2;
    end
   end
end

    
%creation of G matrix
for i = 1:nx 
     for j = 1:ny
        n = j +(i-1)*ny;
        nxm = j + (i-2)*ny;
        nxp = j + (i)*ny;
        nyp = (j+1) + (i-1)*ny; 
        nym = (j-1) +(i-1)*ny;
        
        if i == 1
            Gmap(n,n) = sigma(i,j);
            f(n) = VR3(z);
            
        elseif i == nx
            Gmap(n,n) = sigma(i,j);   
            
        elseif j == 1 %bottom   
            Gmap(n,nxp) = sigma(i+1,j);
            Gmap(n,nxm) = sigma(i-1,j);
            Gmap(n,nyp) = sigma(i,j+1);
            Gmap(n,n) = -(Gmap(n,nxp)+Gmap(n,nxm)+Gmap(n,nyp));
            
        elseif j == ny %top
            Gmap(n,nxp) = sigma(i+1,j);
            Gmap(n,nxm) = sigma(i-1,j);
            Gmap(n,nym) = sigma(i,j-1);
            Gmap(n,n) = -(Gmap(n,nxp)+Gmap(n,nxm)+Gmap(n,nym));
          
        else 
            Gmap(n,nxp) = sigma(i+1,j);
            Gmap(n,nxm) = sigma(i-1,j);
            Gmap(n,nym) = sigma(i,j-1);
            Gmap(n,nyp) = sigma(i,j+1);
            Gmap(n,n) = -(Gmap(n,nxp)+Gmap(n,nxm)+Gmap(n,nym)+Gmap(n,nyp));
            end           
        end
end   

% creating V matrix  
V1 = Gmap\f';

    for i = 1:nx
        for j = 1:ny
            n = j + (i-1)*ny;
            M(i,j) = V1(n);
        end
       %surf(M) %,'linestyle','none')  
    end
% Electron, Force and acceleration calculation    
[Ey,Ex] = gradient(M);
Ex = -Ex/dx;
Ey = -Ey/dy;
Fx = (Ex)*(1.602e-19);
Fy = (Ey)*(1.602e-19);
ax = Fx/m; 
ay= Fy/m; 

%maxwell-boltzmann distribution
vx = randn(numberofelectrons,1) *vth/sqrt(2);
vy = randn(numberofelectrons,1) *vth/sqrt(2);
vx_thermalized = randn(numberofelectrons,1) *vth/sqrt(2);
vy_thermalized = randn(numberofelectrons,1) *vth/sqrt(2);



%For loop for plotting particle trajectory
for t =1:Nt  
    %position, velocity and acceleration calculations usining discretization
    xp = X;
    yp = Y;
   
    xbin = discretize(X,nx);
    ybin = discretize(Y,ny);
    
    axp = zeros(numberofelectrons,1);
    ayp = zeros(numberofelectrons,1);
    
    for i = 1:numberofelectrons
        axp(i) = ax(xbin(i),ybin(i));
        ayp(i) = ay(xbin(i),ybin(i));
    end

    X = X + vx*dt+(1/2)*(axp)*(dt^2);
    Y = Y + vy*dt+(1/2)*(ayp)*(dt^2);
    
    vx = vx + axp*dt;
    vy = vy + ayp*dt;

    %scattering probabilty
    Pscat = 1-exp(-dt/tmn);
    scat = Pscat > rand(numberofelectrons,1);
    vx(scat) = vx_thermalized(scat);
    vy(scat) = vy_thermalized(scat);
      
    %Boundary conditions 
    %inside top and bottom box
    InTopBox = X < 1.2e-7 & X > 0.8e-7 & Y > 0.6e-7;
    InBotBox = X < 1.2e-7 & X > 0.8e-7 & Y < 0.4e-7;
    
    %right and left of the boxes
    reflectxtopright = xp > 1.2e-7 & yp > 0.6e-7;
    reflextxtopleft = xp < 0.8e-7 & yp > 0.6e-7;
    reflectxbotright = xp > 1.2e-7 & yp < 0.4e-7;
    reflextxbotleft = xp < 0.8e-7 & yp < 0.4e-7;
    reflectmiddle = xp < 1.2e-7 & xp > 0.8e-7;

    %reverse xspeed for topbox 
    vx(reflectxtopright & InTopBox) = -vx(reflectxtopright & InTopBox);
    vx(reflextxtopleft & InTopBox) = -vx(reflextxtopleft & InTopBox);
   
    %reverse xspeed for bottom box
    vx(reflectxbotright & InBotBox) = -vx(reflectxbotright & InBotBox);
    vx(reflextxbotleft & InBotBox) = -vx(reflextxbotleft & InBotBox);
    
    %reverse x for bottom box
    X(reflectxbotright & InBotBox) = xp(reflectxbotright & InBotBox);
    X(reflextxbotleft & InBotBox) = xp(reflextxbotleft & InBotBox);
    
    %reverse x for top box
    X(reflectxtopright & InTopBox) = xp(reflectxtopright & InTopBox);
    X(reflextxtopleft & InTopBox) = xp(reflextxtopleft & InTopBox);
    
    %reverse y for top and bottom box
    Y(reflectxbotright & InBotBox) = yp(reflectxbotright & InBotBox);
    Y(reflextxbotleft & InBotBox) = yp(reflextxbotleft & InBotBox);
    Y(reflectxtopright & InTopBox) = yp(reflectxtopright & InTopBox);
    Y(reflextxtopleft & InTopBox) = yp(reflextxtopleft & InTopBox);
    
    %reverse yspeed and y for middle 
    vy(reflectmiddle & (InTopBox|InBotBox)) = -vy(reflectmiddle & (InTopBox|InBotBox));
    Y(reflectmiddle & (InTopBox|InBotBox)) = yp(reflectmiddle & (InTopBox|InBotBox));

    %periodic boundary condition
    xp(X>l) = xp(X>l)-l;
    X(X>l) = X(X>l)-l;
    xp(X<0) = xp(X<0)+l;
    X(X<0) = X(X<0)+l;
     

   %particle hitting top and bottom
   vy(Y>w) = -vy(Y>w);
   vy(Y<0) = -vy(Y<0);
    
   
    v_x = mean(vx);
                            % 3D: C elec/m^3 m/s = A/m^2
                            %     I = J*Area
    Jx(z) = q*n_elec*v_x;   % 2D: C elec/m^2 m/s = A/m
    I = Jx*w;           %     I = J*height
    % n_particle = num of particles/L*W
    % simple example 3D:
    % n_elec = 100 elec/unit vol
    % n_particles = 10 particles/unit vol
    % elec/particle = 10
end
end


%linear fit for R3 resistance
p1 = polyfit(I,VR3,1);


% Time Domain:
% v1 = vin;
% (v2-v1)*G1 + v2*G2 +C*d(v2-v1)/dt + IL = 0
% (dIL/dt)*L= (v3-v2);
% IL = v3*G3
% v4 = -a*v3*G3; 
% (v5-v4)*G4 + v5*G5 = 0

% Frequency Domain:
% v1 = vin; 
% (v2-v1)*G1 + v2*G2 +C*jw(v2-v1) + IL = 0
% jwIL*L = (v2-v3);
% IL = v3*G3 ;
% v4 = -a*v3*G3;
% (v5-v4)*G4 + v5*G5 = 0

% unknowns:
% v1
% v2
% v3
% v4
% v5
% IL

%variables and constants for circuit simulation
G1 = 1/1;
G2 = 1/2;
G3 = 1/abs(p1(1));%update resistor value with p1 value
G4 = 1/0.1;
G5 = 1/1000;
v1 = 0;
v2 = 0;
v3 = 0;
v4 = 0;
v5 = 0;
IL = 0;
a = 100;
Cap = 0.25;
L = 0.2;
n = 6;
step = 200;
G = zeros(n,n);
C = zeros(n,n);
V = zeros(n,1);
F = zeros(n,1);
V_5 = zeros(step,1);
Vin = linspace(-10,10,step);
V_3 = zeros(step,1);
V_1 = zeros(step,1);
freq = linspace(0,100,step);



%Vector for unknowns
V(1,1) = v1;
V(2,1) = v2;
V(3,1) = v3;
V(4,1) = v4;
V(5,1) = v5;
V(6,1) = IL;

%G matrix
G(1,1) = 1;
G(2,1) = -G1;
G(2,2) = G1 + G2;
G(2,6) = 1;
G(3,2) = 1;
G(3,3) = -1;
G(4,3) = G3;
G(4,6) = 1;
G(5,3) = -G3*a;
G(5,4) = 1;
G(6,4) = -G4;
G(6,5) = G4 + G5;

%Capacitance Matrix
C(2,1) = -Cap;
C(2,2) = Cap;
C(3,6) = L;

%F matrix 
F(1,1) = v1;

%DC Sweep over 0.1 to 10
for i = 1:step
   F(1,1) = Vin(i);
    V = G\F;
    V_5(i) = V(5,1);
    V_3(i) = V(3,1);
end

subplot(2,2,1)
plot(Vin,V_5)
hold on
plot(Vin,V_3)
title('DC')
legend('Vo','V3')
ylabel('Voltage(V)')
xlabel('Voltage(V)')


%AC simulation 
 for i = 1:step
    F(1,1) = Vin(i);
    A = G + 1j*2*pi*freq(i)*C;
    V = A\F;
    V_5(i) = V(5,1);
    V_1(i) = V(1,1);
 end
 
 
subplot(2,2,2)
plot(freq,(V_5/V_1))
title('Gain')
ylabel('Gain(V/V)')
xlabel('Frequency')

%Varying Capacitance
Capdistribution = (randn(step) *0.05)+0.25;
w =  linspace(0,pi,step);
gain = zeros(step,1);
for i = 1:step
   F(1,1) = Vin(i);
   C(2,1) = -Capdistribution(i);
   C(2,2) = Capdistribution(i);
   A = G + 1j*2*pi*w(i)*C;
   V = A\F;
   V_5(i) = V(5,1);
   V_1(i) = V(1,1);
   gain(i) = (abs(20*log(V_5(i)/V_1(i))));
end
subplot(2,2,3)
histogram(Capdistribution)
title('C')
ylabel('Distribution')
xlabel('Capacitance(F)')

subplot(2,2,4)
histogram(gain)
title('Gain Histogram')
ylabel('Distribution')
xlabel('Gain(|20log(vout/vin)|)')







