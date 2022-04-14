function [E,F,G] = NLSW_EFGcalc(nx,ny,h,eta,P,Q,dx,dy)
% This function file will calculate the right hand sides of the continuity (E), x-momentum (F), and y-momentum (G) equations.
% You should not have to change this file for CP#2
% INPUT
%		nx - # of points in x direction
%		ny - # of points in y direction
%		h - constant water depth
%		eta - eta (free surface elevation) surface matrix size nx by ny
%		u - u (x-velocity) surface matrix size nx by ny
%		v - v (y-velocity) surface matrix size nx by ny
%		dx - grid size in x-direction
%		dy - grid size in y-direction
%
% OUTPUT
%		E,F,G - size nx by ny

E=zeros(nx,ny);
F=zeros(nx,ny);
G=zeros(nx,ny);

% group variables into temp matrices
H=h+eta;
Pflux=P.^2./H+9.81*H.^2/2;

for i=1:nx
   for j=1:ny
      % evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
      if i==1
         P_x=(P(i+1,j)-P(i,j))/dx;
         Pflux_x=(Pflux(i+1,j)-Pflux(i,j))/dx;
      elseif i==nx
         P_x=(P(i,j)-P(i-1,j))/dx;
         Pflux_x=(Pflux(i,j)-Pflux(i-1,j))/dx;
      else
         P_x=(P(i+1,j)-P(i-1,j))/dx/2;
         Pflux_x=(Pflux(i+1,j)-Pflux(i-1,j))/dx/2;
      end

      E(i,j)=-P_x; 
      F(i,j)=-Pflux_x;
      G(i,j)=0; %-9.81*eta_y-0.5*u_sq_y-0.5*v_sq_y;
      
   end
end

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

