% Compute principal stresses from NSIDC-0670 (20-year average velocities)
% across all Greenland, and their associated uncertainties
% kristin poinar
% univ buffalo
% march 26, 2024
%
%

clear variables

%% Setup
% Load velocity data, set domain and flow law parameters

% Use existing grid strain rates
% load allGreenlandStrainRates_NSIDC0670_3kmSmoothing.mat

% Calculate principal strain rates using Savitzky-Golay filter
dx = 250;
dy = 250;
L = 1e3;

% Domain: Full Greenland
% xlim = [-645125 859875-dx];
% ylim = [-3370125 -640125-dy];
xlim = [-620000 845000-dx];
ylim = [-3370125 -825000-dy];

% Domain: small Helheim
% xlim = [235000 330000];
% ylim = [-2605000 -2550000];

% GPS stations
gps = m_shaperead('~/Library/CloudStorage/Box-Box/gly-kpoinar/gly-kpoinar-projects/FTW-FAM/gis/July2023_imagea-GPS-locations');
gps = cell2mat(gps.ncst);

% Read velocity
vels = readNSIDC_0670(xlim, ylim);
[ep,e1,e3,shear,de] = principalstrainrates_savitzkygolay(vels.vx,vels.vy,dx,dy,vels.ex,vels.ey,L);
d_e1 = de;
d_e3 = de;
clear de

% Calculate strain rates: start with effective strain rate
epseff = sqrt(0.5*(e1.^2 + e3.^2)+shear.^2);

% Set flow law parameter (A) and flow law exponent (n)
n = 3;
T = -10;
A = AofT(T);
% Convert A from Pa-3s-1 to Pa-3 yr-1
A = A * 31557600;

% Calculate principal stresses from principal strain rates and the above
sigma1 = 1/A^(1/n) * epseff.^((1-n)/n) .* e1;
sigma3 = 1/A^(1/n) * epseff.^((1-n)/n) .* e3;



%% Error propagation


% Begin uncertainty calculations: Get errors in effective strain rate and A
%d_epseff = 1/sqrt(2) * d_e1 ./ sqrt(e1.^2+e3.^2);
%d_epseff = d_e1 .* (e1.^2 + e3.^2).^0.25;
% update d_epseff May 2024
h = e1/sqrt(2) * 1./(sqrt(e1.^2+e3.^2+2*shear.^2));
l = shear ./(sqrt(0.5*e1.^2 + 0.5*e3.^2 + shear.^2)); 
d_epseff = d_e1 .* sqrt(h.^2 .* (1+(e3./e1).^2) + l.^2);
d_A = (AofT(-5)-AofT(-15)) * 31557600;
% d_e3 and d_e1 are already delivered from the strain rate calculations


% Derivatives

dsigma1depseff = (1-n)/n * A^(-1/n) * e1 .* epseff.^((1-2*n)/n);
dsigma3depseff = (1-n)/n * A^(-1/n) * e3 .* epseff.^((1-2*n)/n);

dsigma1de1 = A^(-1/n) * epseff.^((1-n)/n);
dsigma3de3 = A^(-1/n) * epseff.^((1-n)/n);

dsigma1dA = 1/n * A^(-1/n-1) * epseff.^((1-n)/n) .* e1;   % 5/24/24 Fixed from ... epseff.^((1-n)/n)
dsigma3dA = 1/n * A^(-1/n-1) * epseff.^((1-n)/n) .* e3;   % 5/24/24 Fixed from ... epseff.^((1-n)/n)

% Errors from each source: uncertainties times derivatives
err1_from_epseff = abs(dsigma1depseff .* d_epseff);
err1_from_e1 = abs(dsigma1de1 .* d_e1);
err1_from_A = abs(dsigma1dA * d_A);

err3_from_epseff = dsigma3depseff .* d_epseff;
err3_from_e3 = dsigma3de3 .* d_e3;
err3_from_A = dsigma3dA * d_A;

% Add in quadrature to get full error in sigma1, sigma3
d_sigma1 = sqrt(err1_from_epseff.^2 + err1_from_e1.^2 + err1_from_A.^2);
d_sigma3 = sqrt(err3_from_epseff.^2 + err3_from_e3.^2 + err3_from_A.^2);

%% Plot the three error sources
figure(3); clf; 
av(1) = subplot(1,3,1); hold on;
pcolor(vels.x/1e3,vels.y/1e3,err1_from_epseff/1e3); shading flat; colorbar;
str = ('\delta\sigma_1 from \dot{\epsilon_E}}');
str = ('$\delta\sigma_1$ (kPa) from eps-eff');
t = title(str,'interpreter','latex');
plot(gps(:,1)/1e3,gps(:,2)/1e3,'.y')
clim([0 100])

av(2) = subplot(1,3,2); hold on;
pcolor(vels.x/1e3,vels.y/1e3,err1_from_e1/1e3); shading flat; colorbar;
title('$\delta\sigma_1$ (kPa) from eps1')
plot(gps(:,1)/1e3,gps(:,2)/1e3,'.y')
clim([0 100])

av(3) = subplot(1,3,3); hold on;
pcolor(vels.x/1e3,vels.y/1e3,err1_from_A/1e3); shading flat; colorbar;
title('$\delta\sigma_1$  (kPa) from A')
plot(gps(:,1)/1e3,gps(:,2)/1e3,'.y')
clim([0 100])
linkaxes(av);

%% Plot the raw uncertainty sources
figure(2); clf;
aw(1) = subplot(1,3,1); hold on;
pcolor(vels.x/1e3,vels.y/1e3,d_epseff); shading flat; colorbar;
title('$\delta \epsilon_E$ (per sec)');%,'interpreter','latex');
plot(gps(:,1)/1e3,gps(:,2)/1e3,'.y')
c1 = nanmedian(d_epseff(:));
caxis([0 10*c1])

aw(2) = subplot(1,3,2); hold on;
pcolor(vels.x/1e3,vels.y/1e3,d_e1); shading flat; colorbar;
title('$\delta \epsilon_1$ (per sec)')
plot(gps(:,1)/1e3,gps(:,2)/1e3,'.y')
c2 = nanmedian(d_e1(:));
caxis([0 10*c2])

aw(3) = subplot(1,3,3); hold on;
pcolor(vels.x/1e3,vels.y/1e3,d_A*ones(size(e1))); shading flat; colorbar;
title('$\delta A$  (Pa-3 yr-1)')
plot(gps(:,1)/1e3,gps(:,2)/1e3,'.y')
caxis([0 d_A])
linkaxes(aw);

%% Plot the total error (relative to the stress itself)
figure(13); clf; 
ax(1) = subplot(1,2,1); 
pcolor(vels.x/1e3,vels.y/1e3,d_sigma1./sigma1); 
shading flat; colorbar; clim([-1 1]*10);
colormap(flipud(cbrewer('div','RdBu',30)))
hold on 
% plot(264.5,-2583.6,'*k')
% plot(268.4,-2583.2,'*k')
plot(gps(:,1)/1e3,gps(:,2)/1e3,'*y')
title('Relative error in $\sigma_1$')

ax(3) = subplot(1,2,2); 
pcolor(vels.x/1e3,vels.y/1e3,d_sigma3./sigma3); 
shading flat; colorbar; clim([-1 1]*10);
linkaxes(ax)
hold on
% plot(264.5,-2583.6,'*k')
% plot(268.4,-2583.2,'*k')
plot(gps(:,1)/1e3,gps(:,2)/1e3,'*y')
title('Relative error in $\sigma_3$')

%
%% Write stress data to netcdf
% Create the ncfile
ncname = 'GreenlandWinterPrincipalStresses_Compressional_0670.nc';
system(['rm ' ncname])  

% Write each variable in turn
% note that the dimensions are not 'row' and 'column' or 'r' and 'c'
% they are 'x' and 'y' and 't' (if I had t)
%
nccreate(ncname,'x', 'Dimensions',{'x', length(vels.x)});
    ncwrite(ncname,'x',vels.x);    
    ncwriteatt(ncname,'x','description','Polar Sterographic x')
    ncwriteatt(ncname,'x','xnits','meters')
    ncwriteatt(ncname,'x','long_name', 'east-west coordinate polar stereographic projection (EPSG:3413) ');

nccreate(ncname,'y', 'Dimensions',{'y', length(vels.y)})
    ncwrite(ncname,'y',vels.y');
    ncwriteatt(ncname,'y','description','Polar Sterographic y')
    ncwriteatt(ncname,'y','units','meters')
    ncwriteatt(ncname,'y','long_name', 'north-south coordinate polar stereographic projection (EPSG:3413) ');

% Write the derived variables
% Note that I have to transpose all my sigma_1, d_sigma_1, etc.
% in order to have their x, y be correctly oriented for QGIS to read them
% in the right place
% Have not figured out how to attach EPSG:3413 projection to these,
% as of now I have to assign the projection manually in QGIS when I open
% this netcdf file to view there.
nccreate(ncname,'sigma_1','Dimensions', {'x', length(vels.x), 'y', length(vels.y)})
   ncwrite(ncname,'sigma_1',sigma1');
   ncwriteatt(ncname,'sigma_1','description','More-compressional principal stress')    
   ncwriteatt(ncname,'sigma_1','units','Pascal')

nccreate(ncname,'d_sigma_1','Dimensions', {'x', length(vels.x), 'y', length(vels.y)})
    ncwrite(ncname,'d_sigma_1',d_sigma1');    
    ncwriteatt(ncname,'d_sigma_1','description','Uncertainty in sigma_1 principal stresses')    
       ncwriteatt(ncname,'d_sigma_1','units','Pascal')



%% Create the ncfile
ncname = 'GreenlandWinterPrincipalStresses_Extensional_0670.nc';
system(['rm ' ncname])  

% Write each variable in turn
% note that the dimensions are not 'row' and 'column' or 'r' and 'c'
% they are 'x' and 'y' and 't' (if I had t)
%
nccreate(ncname,'x', 'Dimensions',{'x', length(vels.x)});
    ncwrite(ncname,'x',vels.x);    
    ncwriteatt(ncname,'x','description','Polar Sterographic x')
    ncwriteatt(ncname,'x','xnits','meters')
    ncwriteatt(ncname,'x','long_name', 'east-west coordinate polar stereographic projection (EPSG:3413) ');

nccreate(ncname,'y', 'Dimensions',{'y', length(vels.y)})
    ncwrite(ncname,'y',vels.y');
    ncwriteatt(ncname,'y','description','Polar Sterographic y')
    ncwriteatt(ncname,'y','units','meters')
    ncwriteatt(ncname,'y','long_name', 'north-south coordinate polar stereographic projection (EPSG:3413) ');

nccreate(ncname,'sigma_3','Dimensions', {'x', length(vels.x), 'y', length(vels.y)})
   ncwrite(ncname,'sigma_3',sigma3');
   ncwriteatt(ncname,'sigma_3','description','More-extensional principal stress')
   ncwriteatt(ncname,'sigma_3','units','Pascal')

nccreate(ncname,'d_sigma_3','Dimensions', {'x', length(vels.x), 'y', length(vels.y)})
    ncwrite(ncname,'d_sigma_3',d_sigma3');    
    ncwriteatt(ncname,'d_sigma_3','description','Uncertainty in sigma_3 principal stresses (Pa)')   
    ncwriteatt(ncname,'d_sigma_3','units','Pascal')


%% Write strain rate data to netcdf
% Create the ncfile
ncname = 'GreenlandWinterPrincipalStrainRates_Compressional_0670.nc';
system(['rm ' ncname])  

% Write each variable in turn
% note that the dimensions are not 'row' and 'column' or 'r' and 'c'
% they are 'x' and 'y' and 't' (if I had t)
%
nccreate(ncname,'x', 'Dimensions',{'x', length(vels.x)});
    ncwrite(ncname,'x',vels.x);    
    ncwriteatt(ncname,'x','description','Polar Sterographic x')
    ncwriteatt(ncname,'x','xnits','meters')
    ncwriteatt(ncname,'x','long_name', 'east-west coordinate polar stereographic projection (EPSG:3413) ');

nccreate(ncname,'y', 'Dimensions',{'y', length(vels.y)})
    ncwrite(ncname,'y',vels.y');
    ncwriteatt(ncname,'y','description','Polar Sterographic y')
    ncwriteatt(ncname,'y','units','meters')
    ncwriteatt(ncname,'y','long_name', 'north-south coordinate polar stereographic projection (EPSG:3413) ');

% Write the derived variables
% Note that I have to transpose all my sigma_1, d_sigma_1, etc.
% in order to have their x, y be correctly oriented for QGIS to read them
% in the right place
% Have not figured out how to attach EPSG:3413 projection to these,
% as of now I have to assign the projection manually in QGIS when I open
% this netcdf file to view there.
nccreate(ncname,'principalstress_1','Dimensions', {'x', length(vels.x), 'y', length(vels.y)})
   ncwrite(ncname,'principalstress_1',e1');
   ncwriteatt(ncname,'principalstress_1','description','More-compressional principal strain rate')    
   ncwriteatt(ncname,'principalstress_1','units','yr-1')

nccreate(ncname,'d_principalstress1','Dimensions', {'x', length(vels.x), 'y', length(vels.y)})
    ncwrite(ncname,'d_principalstress1',d_e1');    
    ncwriteatt(ncname,'d_principalstress1','description','Uncertainty in more-compressional principal strain rate')    
       ncwriteatt(ncname,'d_principalstress1','units','yr-1')


% Create the ncfile
ncname = 'GreenlandWinterPrincipalStrainRates_Extensional_0670.nc';
system(['rm ' ncname])  

% Write each variable in turn
% note that the dimensions are not 'row' and 'column' or 'r' and 'c'
% they are 'x' and 'y' and 't' (if I had t)
%
nccreate(ncname,'x', 'Dimensions',{'x', length(vels.x)});
    ncwrite(ncname,'x',vels.x);    
    ncwriteatt(ncname,'x','description','Polar Sterographic x')
    ncwriteatt(ncname,'x','xnits','meters')
    ncwriteatt(ncname,'x','long_name', 'east-west coordinate polar stereographic projection (EPSG:3413) ');

nccreate(ncname,'y', 'Dimensions',{'y', length(vels.y)})
    ncwrite(ncname,'y',vels.y');
    ncwriteatt(ncname,'y','description','Polar Sterographic y')
    ncwriteatt(ncname,'y','units','meters')
    ncwriteatt(ncname,'y','long_name', 'north-south coordinate polar stereographic projection (EPSG:3413) ');

nccreate(ncname,'principalstress_1','Dimensions', {'x', length(vels.x), 'y', length(vels.y)})
   ncwrite(ncname,'principalstress_1',e1');
   ncwriteatt(ncname,'principalstress_1','description','More-extensional principal strain rate')
   ncwriteatt(ncname,'principalstress_1','units','yr-1')

nccreate(ncname,'d_principalstress1','Dimensions', {'x', length(vels.x), 'y', length(vels.y)})
    ncwrite(ncname,'d_principalstress1',d_e1');    
    ncwriteatt(ncname,'d_principalstress1','description','Uncertainty in more-extensional principal strain rate')   
    ncwriteatt(ncname,'d_principalstress1','units','yr-1')
%
