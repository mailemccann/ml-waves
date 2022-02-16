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
Pflux_F=P.^2./H+9.81*H.^2/2;
Pflux_G=P.*Q./H;

Qflux_F=P.*Q./H;
Qflux_G=Q.^2./H+9.81*H.^2/2;


for i=1:nx
   for j=1:ny
      % evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
      eta_x=0;
      eta_y=0;
      if i==1
         P_x=(P(i+1,j)-P(i,j))/dx;
         Pflux_x=(Pflux_F(i+1,j)-Pflux_F(i,j))/dx;
         Qflux_x=(Qflux_F(i+1,j)-Qflux_F(i,j))/dx;
      elseif i==nx
         P_x=(P(i,j)-P(i-1,j))/dx;
         Pflux_x=(Pflux_F(i,j)-Pflux_F(i-1,j))/dx;
         Qflux_x=(Qflux_F(i,j)-Qflux_F(i-1,j))/dx;
      else
         eta_x=(eta(i+1,j)-eta(i-1,j))/dx/2;
         P_x=(P(i+1,j)-P(i-1,j))/dx/2;
         Pflux_x=(Pflux_F(i+1,j)-Pflux_F(i-1,j))/dx/2;
         Qflux_x=(Qflux_F(i+1,j)-Qflux_F(i-1,j))/dx/2;
      end
      
      if j==1
         Q_y=(Q(i,j+1)-Q(i,j))/dy;
         Pflux_y=(Pflux_G(i,j+1)-Pflux_G(i,j))/dy;
         Qflux_y=(Qflux_G(i,j+1)-Qflux_G(i,j))/dy;
      elseif j==ny
         Q_y=(Q(i,j)-Q(i,j-1))/dy;
         Pflux_y=(Pflux_G(i,j)-Pflux_G(i,j-1))/dy;
         Qflux_y=(Qflux_G(i,j)-Qflux_G(i,j-1))/dy;
      else
         eta_y=(eta(i,j+1)-eta(i,j-1))/dy/2;
         Q_y=(Q(i,j+1)-Q(i,j-1))/dy/2;
         Pflux_y=(Pflux_G(i,j+1)-Pflux_G(i,j-1))/dy/2;
         Qflux_y=(Qflux_G(i,j+1)-Qflux_G(i,j-1))/dy/2;
      end

      E(i,j)=-P_x-Q_y; 
      F(i,j)=-Pflux_x-Pflux_y;
      G(i,j)=-Qflux_x-Qflux_y;
      
   end
end

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

