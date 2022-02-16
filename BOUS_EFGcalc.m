function [E,F,G,Fstar,Gstar] = NLSW_EFGcalc(nx,ny,h,eta,P,Q,dx,dy)
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
%		E,F,G,F*,G* - size nx by ny

E=zeros(nx,ny);
F=zeros(nx,ny);
G=zeros(nx,ny);
Fstar=zeros(nx,ny);
Gstar=zeros(nx,ny);

B=1/15;

% group variables into temp matrices
H=h+eta;
Pflux_F=P.^2./H+9.81*H.^2/2;
Pflux_G=P.*Q./H;

Qflux_F=P.*Q./H;
Qflux_G=Q.^2./H+9.81*H.^2/2;

% precalc some derivatives
P_x=zeros(nx,ny);
Q_y=zeros(nx,ny);
eta_x=zeros(nx,ny);
eta_y=zeros(nx,ny);
eta_xx=zeros(nx,ny);
eta_yy=zeros(nx,ny);
for i=1:nx
   for j=1:ny
      % evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
      if i==1
         P_x(i,j)=(P(i+1,j)-P(i,j))/dx;
         eta_x(i,j)=(eta(i+1,j)-eta(i,j))/dx;
         eta_xx(i,j)=0;
      elseif i==nx
         P_x(i,j)=(P(i,j)-P(i-1,j))/dx;
         eta_x(i,j)=(eta(i,j)-eta(i-1,j))/dx;
         eta_xx(i,j)=0;
      else 
         P_x(i,j)=(P(i+1,j)-P(i-1,j))/dx/2;
         eta_x(i,j)=(eta(i+1,j)-eta(i-1,j))/dx/2;
         eta_xx(i,j)=(eta(i+1,j)-2.*eta(i,j)+eta(i-1,j))/dx^2;
      end
      
      if j==1
         Q_y(i,j)=(Q(i,j+1)-Q(i,j))/dy;
         eta_y(i,j)=(eta(i,j+1)-eta(i,j))/dy;
         eta_yy(i,j)=0;
      elseif j==ny
         Q_y(i,j)=(Q(i,j)-Q(i,j-1))/dy;
         eta_y(i,j)=(eta(i,j)-eta(i,j-1))/dy;
         eta_yy(i,j)=0;
      else 
         Q_y(i,j)=(Q(i,j+1)-Q(i,j-1))/dy/2;
         eta_y(i,j)=(eta(i,j+1)-eta(i,j-1))/dy/2;
         eta_yy(i,j)=(eta(i,j+1)-2.*eta(i,j)+eta(i,j-1))/dy^2;
      end
   end
end

for i=1:nx
   for j=1:ny
      % evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
      if i==1
         Pflux_x=(Pflux_F(i+1,j)-Pflux_F(i,j))/dx;
         Qflux_x=(Qflux_F(i+1,j)-Qflux_F(i,j))/dx;
         Q_xy=       (Q_y(i+1,j)    -Q_y(i,j))/dx;
         eta_xxx= (eta_xx(i+1,j) -eta_xx(i,j))/dx;
         eta_xyy= (eta_yy(i+1,j) -eta_yy(i,j))/dx;
      elseif i==nx
         Pflux_x=(Pflux_F(i,j)-Pflux_F(i-1,j))/dx;
         Qflux_x=(Qflux_F(i,j)-Qflux_F(i-1,j))/dx;
         Q_xy=       (Q_y(i,j)    -Q_y(i-1,j))/dx;
         eta_xxx= (eta_xx(i,j) -eta_xx(i-1,j))/dx;
         eta_xyy= (eta_yy(i,j) -eta_yy(i-1,j))/dx;
      else
         Pflux_x=(Pflux_F(i+1,j)-Pflux_F(i-1,j))/dx/2;
         Qflux_x=(Qflux_F(i+1,j)-Qflux_F(i-1,j))/dx/2;
         Q_xy=       (Q_y(i+1,j)    -Q_y(i-1,j))/dx/2.;
         eta_xxx= (eta_xx(i+1,j) -eta_xx(i-1,j))/dx/2.;
         eta_xyy= (eta_yy(i+1,j) -eta_yy(i-1,j))/dx/2.;
      end
      
      if j==1
         Pflux_y=(Pflux_G(i,j+1)-Pflux_G(i,j))/dy;
         Qflux_y=(Qflux_G(i,j+1)-Qflux_G(i,j))/dy;
         P_xy=       (P_x(i,j+1)    -P_x(i,j))/dy;
         eta_yyy= (eta_yy(i,j+1) -eta_yy(i,j))/dy;
         eta_xxy= (eta_xx(i,j+1) -eta_xx(i,j))/dy;
      elseif j==ny
         Pflux_y=(Pflux_G(i,j)-Pflux_G(i,j-1))/dy;
         Qflux_y=(Qflux_G(i,j)-Qflux_G(i,j-1))/dy;
         P_xy=       (P_x(i,j)    -P_x(i,j-1))/dy;
         eta_yyy= (eta_yy(i,j) -eta_yy(i,j-1))/dy;
         eta_xxy= (eta_xx(i,j) -eta_xx(i,j-1))/dy;
      else
         Pflux_y=(Pflux_G(i,j+1)-Pflux_G(i,j-1))/dy/2;
         Qflux_y=(Qflux_G(i,j+1)-Qflux_G(i,j-1))/dy/2;
         P_xy=       (P_x(i,j+1)    -P_x(i,j-1))/dy/2.;
         eta_yyy= (eta_yy(i,j+1) -eta_yy(i,j-1))/dy/2.;
         eta_xxy= (eta_xx(i,j+1) -eta_xx(i,j-1))/dy/2.;
      end

      E(i,j)=-P_x(i,j)-Q_y(i,j); 
      F(i,j)=-Pflux_x-Pflux_y+B*9.81*h(i,j)^3*(eta_xxx+eta_xyy);   % (note that these need to be updated for variable bathymetry)
      G(i,j)=-Qflux_x-Qflux_y+B*9.81*h(i,j)^3*(eta_yyy+eta_xxy);    % (note that these need to be updated for variable bathymetry)
      
      Fstar(i,j)=(B+1/3)*h(i,j)^2*Q_xy;  % (note that these need to be updated for variable bathymetry)
      Gstar(i,j)=(B+1/3)*h(i,j)^2*P_xy;  % (note that these need to be updated for variable bathymetry)    
      
   end
end

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

