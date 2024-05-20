function vels = readNSIDC_0670(xlim, ylim)
%
% 
% This function returns the average ice sheet velocity for any place in Greenland
% from the NASA MEaSUREs dataset at
% https://nsidc.org/data/nsidc-0670/versions/1
% 
% script by kristin poinar
%
%
%%%%%%%%%%%%%%%%%%%% INPUTS INFO %%%%%%%%%%%%%%%%%%%%
% xlim      desired locations in x, meters (north polar stereographic / EPSG:3413)
% ylim      desired locations in y, meters (north polar stereographic / EPSG:3413)
% These are each two-element vectors, [xmin xmax] and [ymin ymax]
% for example, for Jakobshavn, it could be
%   xlim = [-217000 -15000] 
%   ylim = [-2330000 -2220000]
%
%%%%%%%%%%%%%%%%%%%% FILES NEEDED %%%%%%%%%%%%%%%%%%%%
% greenland_vel_mosaic250_vx_v1.tif
% greenland_vel_mosaic250_vy_v1.tif
% greenland_vel_mosaic250_ex_v1.tif
% greenland_vel_mosaic250_ey_v1.tif
% from the url in the header
% this is about 1 GB of data

%%%%%%%%%%%%%%%%%%%% OUTPUTS INFO %%%%%%%%%%%%%%%%%%%%
% vx        velocity in x direction, m/yr
% vy        velocity in y direction, m/yr
% V         speed, m/yr
% ex        uncertainty (error) in x direction, m/yr
% ey        uncertainty (error) in y direction, m/yr
% theta     flow azimuth (radians)
% These variables are all fields in the structure vels

filepath = './nsidc-0670/';  % the directory that contains the tif files
filestem = 'greenland_vel_mosaic250_';        % unless you changed the file names, no need to alter this line
        dx = 250;
        dy = 250;
        vels.x = -645125 : dx : (859875-dx);
        vels.y = -3370125: dy :(-640125-dy);
        vels.vx = flipud(imread(strcat(filepath,filestem,'vx_v1.tif')));
        vels.vy = flipud(imread(strcat(filepath,filestem,'vy_v1.tif')));
        vels.ex = flipud(imread(strcat(filepath,filestem,'ex_v1.tif')));
            % nsidc.org/data/nsidc-0670 Error sources: Add 3% of vel mag
            vels.ex = sqrt(vels.ex.^2 + 0.03*vels.vx.^2); 
        vels.ey = flipud(imread(strcat(filepath,filestem,'ey_v1.tif')));
            % nsidc.org/data/nsidc-0670 Error sources: Add 3% of vel mag
            vels.ey = sqrt(vels.ey.^2 + 0.03*vels.vy.^2); 
            
        % Clip to the x, y window requested
        j = find(vels.x > xlim(1) & vels.x < xlim(end));
        i = find(vels.y > ylim(1) & vels.y < ylim(end));
        vels.vx = vels.vx(i,j);
        vels.vy = vels.vy(i,j);
        vels.ex = vels.ex(i,j);
        vels.ey = vels.ey(i,j);
        vels.x = vels.x(j);
        vels.y = vels.y(i);
        %
        % Remove flagged pixels of no data
        temp = find(vels.vx < -1e9);
        vels.vx(temp) = NaN; vels.vy(temp) = NaN; vels.ex(temp) = NaN; vels.ey(temp) = NaN;
        clear dx dy temp
        %
        %
        % Total vel
        vels.V = sqrt(vels.vx.^2 + vels.vy.^2);
        %
        % Flow azimuth
        vels.theta = atan(vels.vy./vels.vx);