%clear all; close all; clc

%Now add 0% rock fraction to the ice shell to make melt posotively bouyant
rho_shell=0.92;
%Now with variable fluid viscosity with eta_l_ref_defined at T=176
itmax  = 7e5;
itplot = 7e2;
nz = 128;      % number of cell centers

%Saved time steps
Phi_saved=zeros(nz+2,1+(itmax/itplot));
T_saved=Phi_saved;
h_saved=zeros(1,1+(itmax/itplot));
t_saved=h_saved;
C_saved=Phi_saved;

% set physical parameter values
Kref=1;  %Nomalized reference permeability
xiref=1; %Normalized reference compaction viscosity
Lref=sqrt(Kref*xiref);   %Normalized length scale = 1

nk  = 3;   % permeability phi exponent
nxi = -1;  % compaction viscosity exponent
kmin=1e-8; %Minimum permability for stability

Vsed=0;  % Growth rate of ice shell determined through thermal evolution

Hinit=1e-1; %Non dimensional initial thickness of ice shell
H=Hinit;
phic=0.5;  %Melt fraction at base of layer 

dy = 1/nz;  %
y  = -0.5*dy:dy:1+0.5*dy;  %Y goes from 0 to 1, is non changing length as opposed to z
ycf=0:dy:1;  % cell faces; ycf(1)=0 corresponds to sedimenting boundary
CFL     = 0.125;  %Stability parameter

% -- Anonymous Functions --
Sf  = @(S) (S(1:end-1)+S(2:end))/2;
phif  = @(phi) (phi(1:end-1)+phi(2:end))/2;
% dirichlet = @(v) [2*vb - v(2) v(2:end-1) v(end-1)-2*vb];
neumann = @(v) [v(2) v(2:end-1) v(end-1)];
boundary = @(v) [phic v(2:end-1)  -v(end-1)];
bc.Snz=0;
bc.S0=0;
boundS = @(v) [(2 * bc.S0 - v(2)) v(2:nz + 1) (2 * bc.Snz - v(nz+1))];

% -- Initial condition -
% Initialize accumulators
t  = 0;
it = 1;
itp = 1;
phi = phic*ones(1,nz + 2);

%Material properties with real units
%Thermal
Q=33000;
R=8.31;
Ts=38; %Eutectic temperature of ammonia
Tm=230;
Cp=1600;
kk=3;
rho=920;
kappa=kk/(rho*Cp);
Lh=330000; %Latent heat

%Add in effective specific heat due to melt fraction
Cp_eff=zeros(1,nz+2)+Cp;

%Melt Migration constants
k0=1e-9; %reference permability
eta=1e14; %reference bulk viscosity
etal=1e-3;%viscosity of fluid
drho=80; %bouyancy ratio
g=0.86; %gravity

%Non-Dimensional conversions,
Lref_real=sqrt(k0*eta/etal);  %length
Sref_real=drho*g*k0/etal; %velocity
tref_real=Lref_real/Sref_real; %time

D=H*Lref_real; %Ice Shell thickness in m

kappa_n0=kappa/(Lref_real*Sref_real); %Non-dimensional thermal diffusivity

%Temparature initial condition
%Conductive Profile
T=linspace(Tm,Ts,nz+2);

%Calculate liquidous composition as a function of temperature
B=53.8; A=650;
C1=T-273.1;
Xl=(-B+sqrt(B^2-(4*A.*C1)))./(2*A); %Ammonia concentration at liquidous
C0=0.1; %BULK AMMONIA CONCENTRATION
C=C0+zeros(1,nz+2); %Allow ammonia concentration to vary with depth and time
%Calculate fluid density
a1=0.0017; b1=0.000013;
rho_l=1.0+(a1.*C1)+b1.*(C1).^2;

%Initial melt fraction defined from thermal equilibrium
phi=C0./Xl;
phieq=phi;
phi(phi>0.5)=0.5;
[junk,i_solidus]=min(abs(T-176));
if T(i_solidus)>176
    ip=0;
else
    ip=1;
end
i_solidus=i_solidus-ip;
phi(i_solidus:i_solidus+9)=linspace(phi(i_solidus),0.0001,10);
phi(i_solidus+9:end)=0.0001;
phi = boundary(phi);
captured=0;
Phi_saved(:,1)=phi;
C_saved(:,1)=C;

%Restart functions
restart=0;
fid = fopen('output.dat','wt');
if restart>0
 load restart.mat  phi Vsed H time
 t=time(end);
 Hinit=H;
end
time=zeros(itmax,1);

gamma_f=zeros(1,nz+2);
gamma_f_crazy=gamma_f;

dt=1e-4;


% -- Time integration --
while it <= itmax    
    
    % 1) Redefine z
    dz=H*dy; %These are constant with time but let's leave it for now
    lz=H;
    z = -0.5*dz:dz:lz+0.5*dz; %cell centers
           
    % 2) Identify solidus region
    [junk,i_solidus]=min(abs(T-176));
    if T(i_solidus)>176
       ip=0;
    else
       ip=1;
    end
    i_solidus=i_solidus-ip;  %Puts it and the warm side of the solidus
           
     % 3) Calculate phieq from previous T and C values
     C1=T-273.1;
     Xleq=(-B+sqrt(B^2-(4*A.*C1)))./(2*A);
     dXldT=-(B^2-4*A*C1).^-0.5;
     
     %Determine dPhidT
     dPhidT=(-C./(Xl.^2)).*dXldT; 
     dPhidT(T<184)=C(i_solidus-1)./Xl(i_solidus-1)./8; %Spreads solidus front out over 8 degrees
     dPhidT(T<176)=0;
     
     phieq_new=C./Xleq;
     %Adjust phieq near solidus region
     phieq_new(i_solidus:i_solidus+9)=linspace(phi(i_solidus),0.0001,10);
     phieq_new(i_solidus+9:end)=0.0001;
        
    % 4) Calculate Permeabilities and viscosities based on old T and phi
%    Calculate velocity field at cell faces
%    phi    = neumann(phi);
    phi    = boundary(phi);
    %Permeability
    k(1:nz+2,1)=Kref*(abs(phi)).^nk;
    
    %Create effective permeability which includes fluid viscosity and
    %density
    %k=k.*exp((T-273.1)./10.55);
    rho_l=1.0+(a1.*C1)+b1.*(C1).^2;
    rho_l(rho_l>1.05)=1.05;
    k(1:nz+2,1)=Kref*((abs(phi)).^nk).*exp((T-273.1)./10.55).*((rho_l-rho_shell)./0.08);
    kzero= k<kmin;    
    k(kzero)=kmin;
    
    %Compaction viscosity 
    xiref=exp((Q/R).*((1./T)-(1/273.1)));    %account for temperature effects    
    xi(1:nz+2,1)= xiref.*((abs(phi)+1d-6).^nxi).*(exp(-22.*phi)); %From Arakawa

    
    %Now redefine Lref as a vector modified by the density contrast
    %See equations showing how S is non-dimensionalized and solved
    Lref=1./sqrt(((rho_l-rho_shell)./0.08));
    
    %xi(1:nz+2,1)= xiref.*((abs(phi)+1d-8)).^nxi;
    xi(xi>10^8)=10^8;
    xi(xi<10^-6)=10^-6;
        
    % 5) Calculate segregation flux
    [S] = Helmholtz_Sv222noah(bc,k,xi,Lref,nz,lz);  %assumes drho*g=1
    %[S] = Helmholtz_Sv222(bc,k,xi,Lref,nz,lz); 

    Sref=ones(1,(nz+2));
    Sref(T<176)=0; %Now flow below solidus
    S=Sref.*S;
    
    p =-((xi(2:nz+1)'+xi(3:nz+2)').*S(3:nz+2) ...
        +(xi(2:nz+1)'-xi(3:nz+2)').*S(2:nz+1) ...
        -(xi(1:nz)'+xi(2:nz+1)').*S(1:nz))/(4*dz);
    p=[p(1) p(1:end) p(end)];
    
    phicf=phif(phi);
    Scf=Sf(S);
    
    % 6) Calculate how quickly boundary  moves
    Vsed=kk*(T(3)-T(4))/(dz*Lref_real)/(Lh*rho);
    Vsed=Vsed/Sref_real;
    VH=Vsed+Scf(1);
%    VH=0;
    Vlcf=VH+(1-phicf).*Scf./(phicf+1e-8);
    VY=(Vlcf-VH*ycf)/H; %Remember Y is non-dimensional depth, %VH*ycf is psuedo advective term
    
    % 7) Determine timestep
    dtcfl = CFL * dz / max(abs(Vlcf));
    dt=dtcfl;
    %Time step needs to be minimum between thermal and melt migration stability
    dt=min(0.2*abs(dz^2/max(kappa_n0)),dt); 
    
    % 8) Calculate new temperature 
    dTdt=(gradient(gradient(T)).*kappa_n0./(dz^2))./(1+(Lh/Cp).*dPhidT);
    T=T+(dt.*dTdt);  %gamma is freezing production rate           
    T(T>Tm)=Tm;
    T(1)=Tm; T(end)=Ts; %Reinforce boundary conditions
    
    % 8) Calculate Crazy Gamma
    term1=[VY,0].*H.*gradient(T)./dz;
    gamma_f_crazy=(-phi./Xl).*dXldT.*(dTdt+term1);
    gamma_f_crazy(i_solidus:end)=0;
             
    % 9) Calculate change in phi and ammonia concentration
    dphidt = [0 VanLeer1D(phi,VY,dy,dt) 0];
    dCdt = [0 VanLeer1D(phi.*Xl,VY,dy,dt) 0];%+(Xl.*dphidt);
    %Ensure no artificial enrichment at solidus boundary
    dCdt(i_solidus-3:end)=0;
   
    % 10) Calculate new phi and C
    phi = phi + dt * (dphidt)+dt*gamma_f_crazy;
    phi(phi>0.5)=0.5;
    phi(T<176)=0.0001;
    C=C+dt*(dCdt);
    C(C<0)=0;

    % 11) Calculate new concentration in the liquid
    Xl=C./phi;
    Xl(Xl>0.36)=0.36;
      
     % 12) Interpolate vectors due to ice shell
       x=linspace(1-(VH*dt/dz),nz+2,nz+2);
       C=interp1([1-(VH*dt/dz),(1:nz+2)],[C0,C],x);
       phi=interp1([1-(VH*dt/dz),(1:nz+2)],[phieq(1),phi],x);
       T=interp1([1-(VH*dt/dz),(1:nz+2)],[T(1),T],x);
     
    % 13) Check for erros
    if isnan(phi)
        error('NaN')
    end
    if (min(T)<(Ts-1))
        error('Negative Temp')
    end
    if C>1
        error('Composition gone bad')
    end
    %Add phi minimum
    phi(phi<0.0001)=0.0001;
    phi=boundary(phi);
    
    
    % 14) Update accumulators
    t  = t + dt;
    it = it + 1;
     H=H+dt*VH;
    
    %Prints out timestep, ice thickness (m), and time (years)
    if (it/100)==round(it/100)
        fprintf('it=%5i   %10.3e   %10.3e\n',it,H*Lref_real,t*tref_real/(pi*10^7))
    end
   
    % plotting
    %phimax=max(phi);
     if round(it/itplot)==it/itplot
        
        itp=itp+1;
        Phi_saved(:,itp)=phi;
        T_saved(:,itp)=T;
        t_saved(itp)=t;
        h_saved(itp)=H;
        C_saved(:,itp)=C;
    end
end

csvwrite(['Phi.Noah.csv'],Phi_saved);

%H=height(end);
save restart.mat phi 
 fclose(fid);
 
figure
subplot(1,3,2),plot(Phi_saved(2:nz+1,1),linspace(Lref_real*Hinit/1000,0,nz),'k','Linewidth',3)
hold on
set(gca,'Fontsize',24,'Ydir','reverse')
xlabel('\phi')
ylabel('Depth (km)')
subplot(1,3,1),plot(T_saved(2:nz+1,2),linspace(Lref_real*Hinit/1000,0,nz),'k','Linewidth',3)
hold on
set(gca,'Fontsize',24,'Ydir','reverse')
xlabel('T')
ylabel('Depth (km)')
subplot(1,3,3),plot(C_saved(2:nz+1,1),linspace(Lref_real*Hinit/1000,0,nz),'k','Linewidth',3)
hold on
set(gca,'Fontsize',24,'Ydir','reverse')
xlabel('\chi_{NH3}')
ylabel('Depth (km)')
for i=1:10  
    subplot(1,3,2), plot(Phi_saved(2:nz+1,i*100),linspace(h_saved(i*100)*Lref_real/1000,0,nz),'Linewidth',2)
    subplot(1,3,1), plot(T_saved(2:nz+1,i*100),linspace(h_saved(i*100)*Lref_real/1000,0,nz),'Linewidth',2)
    subplot(1,3,3), plot(C_saved(2:nz+1,i*100),linspace(h_saved(i*100)*Lref_real/1000,0,nz),'Linewidth',2)
    pause(0.5)
end


