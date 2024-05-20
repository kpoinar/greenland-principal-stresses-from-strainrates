function [ep,e1,e3,shear,de] = principalstrainrates_savitzkygolay(U,V,dx,dy,dU,dV,L)
%
% Principal strain rates (e1,e3) based on velocity field (U,V)
% Based on Harper & Humphrey 1998
% Kristin Poinar April 2016
%
% ADAPTED by kristin on 26 november 2023
% to use Savitzky-Golay filter as described by Brent Minchew
% in doi:10.1017/jog.2018.47
% Processes controlling the downstream evolution of ice rheology in glacier 
% shear margins: case study on Rutford Ice Stream, West Antarctica 
% BRENT M. MINCHEW,1* COLIN R. MEYER,2 ALEXANDER A. ROBEL,3,4 
% G. HILMAR GUDMUNDSSON,1† MARK SIMONS3
%
% "From the transects of horizontal speed, we calculated the lateral shear 
% strain rates by using a second-order Savitzky– Golay filter and a 2 km 
% (≈h) window. The Savitzky–Golay filter is an efficient method for fitting 
% piecewise polynomials of a given order within a window of given size."
%
% From Matlab help for sgolay,
% under Savitzky-Golay Differentiation:
% Estimate the first three derivatives of the sinusoid using the 
% Savitzky-Golay method. Use 25-sample frames and fifth order polynomials. 
% Divide the columns by powers of dt to scale the derivatives correctly.
% [b,g] = sgolay(5,25);
% dx = zeros(length(x),4);
% for p = 0:3
%   dx(:,p+1) = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
% end
%
% I want a first derivative, so p = 1.
%
%
% Get sizes for allocation
nx = size(U,2);
ny = size(U,1);
%
dudx = 0*U;
dvdy = 0*V;
dudy = 0*U;
dvdx = 0*V;
%
% Construct the filter:
lfiltx = round(L/dx);
    % filter width needs to be odd
    if ~mod(lfiltx,2), lfiltx = lfiltx+1; end
[~,gx] = sgolay(1,lfiltx);

lfilty = round(L/dy);
    % filter width needs to be odd
    if ~mod(lfilty,2), lfilty = lfilty+1; end
[~,gy] = sgolay(1,lfilty);
%
%
% Differentiate U and V to get dudx and dvdx and dudy and dvdy
% Apply the filter to each ROW of U, V to get x derivatives
p = 1;
for cc = 1:ny
  dudx(cc,:) = conv(U(cc,:), factorial(p)/(-dx)^p * gx(:,p+1), 'same');
  dvdx(cc,:) = conv(V(cc,:), factorial(p)/(-dx)^p * gx(:,p+1), 'same');
%   figure(10); clf; hold on; plot(U(cc,:))
%                    plot(conv(U(cc,:), factorial(0)/(-dx)^0 * gx(:,1), 'same'))
%                    plot(dudx(cc,:)*500)
%   title(sprintf('%d of %d',cc,nx))
%   pause(0.05)
end
%
%
%
% Differentiate U and V to get dudx and dvdx and dudy and dvdy
% Apply the filter to each COLUMN of U, V to get y derivatives
for cc = 1:nx
  dudy(:,cc) = conv(U(:,cc), factorial(p)/(-dy)^p * gy(:,p+1), 'same');
  dvdy(:,cc) = conv(V(:,cc), factorial(p)/(-dy)^p * gy(:,p+1), 'same');
end


%
%
% Calculate shear strain rate, angle with the vertical, and principal
% strain rates, based on easy formulas from Cuffey & Paterson (A.6-A.7)
shear = 0.5*(dudy + dvdx);
% Principal strain axes (Harper & Humphrey Eqns 3-4) Magnitudes (Eqn 3a,3b)
e1 = 0.5*(dudx + dvdy) - sqrt(0.25*(dudx-dvdy).^2 + shear.^2);
%e2 = 0;  % no strain in the vertical direction
e3 = 0.5*(dudx + dvdy) + sqrt(0.25*(dudx-dvdy).^2 + shear.^2);
%
% primary principal strain rate: the greater of e1 or e3
ep = e1;
ff = abs(e3)>abs(e1);
ep(ff) = e3(ff);
%
% How to get error?  I might assume it's the same as linear differentiation, 
% d_e = 1/2/dx * sqrt(dU.^2 + dV.^2)
de = 1/2/dx * sqrt(dU.^2 + dV.^2);
