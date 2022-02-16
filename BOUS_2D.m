%function CP2_Heun_iter(dx,dy,endx,endy,endt,h,courant)
% This function solves the nonlinear shallow water equations using Euler's method.  As set-up, it gives no output,
% although it creates an animation plot of the free surface.
% INPUT
%			dx - the grid size in the x-direction, in meters
%			dy - the grid size in the y-direction, in meters
%			endx - the total domain length in the x-direction, in meters
%			endy - the total domain length in the y-direction, in meters
%			endt - the total time to run the simulation, in seconds
%			h - the constant water depth, in meters
%			courant - the courant number, controls the time step
% As an example call, use:  >>CP2_Heun_iter(1,1,50,50,100,1,0.5)

clear all
ho=1.;
dx=0.1*ho;
dy=0.1*ho;
endx=20;
endy=20;
endt=20;
courant=0.25;

clf  % clear the current figure

x=[0:dx:endx]';  % create x column vector
y=[0:dy:endy]';  % create y column vector

dt=courant*dx/sqrt(9.81*ho);  % calculate time step
t=[0:dt:endt]';  % create t column vector

nx=length(x);  % total number of grid points in x
ny=length(y);  % total number of grid points in y
nt=length(t);  % total number of time steps

eta=zeros(nx,ny,2);  % initialize the eta, u, and v matrices (see CP#1 for explanation of this step)
P=zeros(nx,ny,2);
Q=zeros(nx,ny,2);
Us=zeros(nx,ny,2);
Vs=zeros(nx,ny,2);
h=zeros(nx,ny);
h(:,:)=ho;  % constant water depth
        
eta_NSW=zeros(nx,ny,2);  % initialize the eta, u, and v matrices (see CP#1 for explanation of this step)
P_NSW=zeros(nx,ny,2);
Q_NSW=zeros(nx,ny,2);

% Initial condition for eta, u and v are zero everywhere for the initial condition
amp=0.15;
for i=1:nx
    for j=1:ny
        loc=sqrt((x(i)-endx/2)^2+(y(j)-endy/2)^2);
        eta(i,j,1)=amp*exp(-1.0*loc^2);
    end
end

surf(x,y,eta(:,:,1)')  %surface plot of the initial free surface
caxis([-.1 .2]);  % set the color axis limits
axis([0 endx 0 endy -0.1 .1])  % axis limits of plot
view(-30,60)% rotate the surface so that we are looking from an angle, note that view(0,90) is the plan view
shading faceted  %  modify the surface shadings
pause(0.001)

eta_previous=eta(:,:,1);
eta_NSW=eta;

Epn0=zeros(nx,ny);
Epn1=zeros(nx,ny);
Epn2=zeros(nx,ny);

Fpn0=zeros(nx,ny);
Fpn1=zeros(nx,ny);
Fpn2=zeros(nx,ny);

Gpn0=zeros(nx,ny);
Gpn1=zeros(nx,ny);
Gpn2=zeros(nx,ny);

Fstar0=zeros(nx,ny);
Fstar1=zeros(nx,ny);
Fstar2=zeros(nx,ny);

Gstar0=zeros(nx,ny);
Gstar1=zeros(nx,ny);
Gstar2=zeros(nx,ny);

Epn0_NSW=zeros(nx,ny);
Epn1_NSW=zeros(nx,ny);
Epn2_NSW=zeros(nx,ny);

Fpn0_NSW=zeros(nx,ny);
Fpn1_NSW=zeros(nx,ny);
Fpn2_NSW=zeros(nx,ny);

Gpn0_NSW=zeros(nx,ny);
Gpn1_NSW=zeros(nx,ny);
Gpn2_NSW=zeros(nx,ny);

% Tri-diagonal coef's (note that these need to be updated for variable
% bathymetry)

B=1/15;
Ax= -1*(B+1/3)*h.^2/dx^2;
Bx=1+2*(B+1/3)*h.^2/dx^2;
Cx= -1*(B+1/3)*h.^2/dx^2;

% Correct for boundaries
Ax(1,:)=0;
Bx(1,:)=1;
Cx(1,:)=0;

Ax(nx,:)=0;
Bx(nx,:)=1;
Cx(nx,:)=0;

Ay= -1*(B+1/3)*h.^2/dy^2;
By=1+2*(B+1/3)*h.^2/dy^2;
Cy= -1*(B+1/3)*h.^2/dy^2;

% Correct for boundaries
Ay(:,1)=0;
By(:,1)=1;
Cy(:,1)=0;

Ay(:,ny)=0;
By(:,ny)=1;
Cy(:,ny)=0;


for n=1:nt-1  % start the time marching
    ['Time step ', num2str(n),' of ', num2str(nt)]  % plot some information to the command window
    
    % Predictor Step
    [Epn0,Fpn0,Gpn0,Fstar0,Gstar0] = BOUS_EFGcalc(nx,ny,h,eta(:,:,1),P(:,:,1),Q(:,:,1),dx,dy);
    
    [Epn0_NSW,Fpn0_NSW,Gpn0_NSW] = NLSW_EFGcalc(nx,ny,h,eta(:,:,1),P(:,:,1),Q(:,:,1),dx,dy);
    
    if n<3
        for i=1:nx
            for j=1:ny
                eta(i,j,2)=eta(i,j,1)+Epn0(i,j)*dt;
                Us(i,j,2)=Us(i,j,1)+Fpn0(i,j)*dt+Fstar0(i,j)-Fstar1(i,j);
                Vs(i,j,2)=Vs(i,j,1)+Gpn0(i,j)*dt+Gstar0(i,j)-Gstar1(i,j);
                
                eta_NSW(i,j,2)=eta(i,j,1)+Epn0_NSW(i,j)*dt;
                P_NSW(i,j,2)=P(i,j,1)+Fpn0_NSW(i,j)*dt;
                Q_NSW(i,j,2)=Q(i,j,1)+Gpn0_NSW(i,j)*dt;
                
            end
        end
    else
        for i=1:nx
            for j=1:ny
                eta(i,j,2)=eta(i,j,1)+dt/12*(23*Epn0(i,j)-16*Epn1(i,j)+5*Epn2(i,j));
                Us(i,j,2) =Us(i,j,1) +dt/12*(23*Fpn0(i,j)-16*Fpn1(i,j)+5*Fpn2(i,j))+Fstar0(i,j)-Fstar1(i,j);
                Vs(i,j,2) =Vs(i,j,1) +dt/12*(23*Gpn0(i,j)-16*Gpn1(i,j)+5*Gpn2(i,j))+Gstar0(i,j)-Gstar1(i,j);
                
                eta_NSW(i,j,2)=eta(i,j,1)+dt/12*(23*Epn0_NSW(i,j)-16*Epn1_NSW(i,j)+5*Epn2_NSW(i,j));
                P_NSW(i,j,2)=P(i,j,1)+dt/12*(23*Fpn0_NSW(i,j)-16*Fpn1_NSW(i,j)+5*Fpn2_NSW(i,j));
                Q_NSW(i,j,2)=Q(i,j,1)+dt/12*(23*Gpn0_NSW(i,j)-16*Gpn1_NSW(i,j)+5*Gpn2_NSW(i,j));
            end
        end
    end

    % Bous BC's
    % left wall
    eta(1,:,2)=eta(2,:,2);
    Us(1,:,2)=-Us(2,:,2);
    Vs(1,:,2)=Vs(2,:,2);
    
    %right wall
    eta(nx,:,2)=eta(nx-1,:,2);
    Us(nx,:,2)=-Us(nx-1,:,2);
    Vs(nx,:,2)=Vs(nx-1,:,2);
    
    %bottom wall
    eta(:,1,2)=eta(:,2,2);
    Us(:,1,2)=Us(:,2,2);
    Vs(:,1,2)=-Vs(:,2,2);
    
    %top wall
    eta(:,ny,2)=eta(:,ny-1,2);
    Us(:,ny,2)=Us(:,ny-1,2);
    Vs(:,ny,2)=-Vs(:,ny-1,2);

    % Tri-diag solver to get P from Us  
    for j=1:ny
        P(:,j,2)=tridiag(Bx(:,j),Ax(:,j),Cx(:,j), Us(:,j,2));
    end
    
    % Tri-diag solver to get Q from Vs
    for i=1:nx
        Q(i,:,2)=tridiag(By(i,:),Ay(i,:),Cy(i,:), Vs(i,:,2));
    end

    % NSW BC's
    % left wall
    eta_NSW(1,:,2)=eta_NSW(2,:,2);
    P_NSW(1,:,2)=-P_NSW(2,:,2);
    Q_NSW(1,:,2)=Q_NSW(2,:,2);
    
    %right wall
    eta_NSW(nx,:,2)=eta_NSW(nx-1,:,2);
    P_NSW(nx,:,2)=-P_NSW(nx-1,:,2);
    Q_NSW(nx,:,2)=Q_NSW(nx-1,:,2);
    
    %bottom wall
    eta_NSW(:,1,2)=eta_NSW(:,2,2);
    P_NSW(:,1,2)=P_NSW(:,2,2);
    Q_NSW(:,1,2)=-Q_NSW(:,2,2);
    
    %top wall
    eta_NSW(:,ny,2)=eta_NSW(:,ny-1,2);
    P_NSW(:,ny,2)=P_NSW(:,ny-1,2);
    Q_NSW(:,ny,2)=-Q_NSW(:,ny-1,2);
    

    % This code is for corrector - for a semi-implicit method
%     e_a=1;  % initialize error to some large value to start while loop
%     e_tol=1e-6;  % error tolerance
%     count=0;
%     Ec=Ep;
%     while e_a>e_tol
%         
%         count=count+1;
%         eta_old=eta(:,:,2);  % place values from last iteration into temp variable
%         u_old=u(:,:,2);      % we need to do this so we can calculate the error, ea
%         v_old=v(:,:,2);
%         
%         
%         % Corrector Step
%         eta_previous=eta(:,:,1);
%         [Ec,Fc,Gc] = EFGcalc(nx,ny,h,eta(:,:,2),u(:,:,2),v(:,:,2),dx,dy);
%         
%         
%         for i=1:nx
%             for j=1:ny
%                 eta(i,j,2)=eta(i,j,1)+0.5*(Ep(i,j)+Ec(i,j))*dt;
%                 u(i,j,2)=u(i,j,1)+0.5*(Fp(i,j)+Fc(i,j))*dt;
%                 
%             end
%         end
%         % enforce bc's
%         eta(1,:,2)=eta(2,:,2);
%         u(1,:,2)=-u(2,:,2);
%         eta(nx,:,2)=eta(nx-1,:,2);
%         u(nx,:,2)=-u(nx-1,:,2);
%         
%         e_a=sum(abs((sum(abs(eta_old))-sum(abs(eta(:,:,2))))/sum(abs(eta(:,:,2))))); % error calc
%         
%         if count>20
%             e_a=0;
%         end
%         
%     end
    
    % Plot the newly calculate surface
    clf
    subplot(4,1,1:3)
    surf(x,y,eta(:,:,2)')
    % caxis([-.1 .1]);
    xlabel('x (m)')
    ylabel('y (m)')
    %axis equal
    axis([0 endx 0 endy -amp/2 1.2*amp])
    view(-30,60)
    shading faceted
    
    subplot(4,1,4)
    plot(x(2:nx-1),P(2:nx-1,round(ny/2),2)'-P_NSW(2:nx-1,round(ny/2),2)')
    % caxis([-.1 .1]);
   % axis([0 endx -amp/2 1.2*amp])
    grid on
    pause(.001)
    
    % Shift all data one time level back.  This is a trick that allows use to use much smaller matrix
    % sizes for these variables.  We are calculating these values at all x,y and t points, so, the "total"
    % calculated domains are nx by ny by nt.  However, since we only need two time levels to perform
    % the Eulers method, we only keep two time levels in memory.  Memory savings here are great.  For example, if
    % we to calculate 5000 time steps on grid with nx=100 and ny=100, we save 1.2 GB (gigabytes) of memory space
    % by only keeping two levels in memory.
    eta(:,:,1)=eta(:,:,2);
    P(:,:,1)=P(:,:,2);
    Q(:,:,1)=Q(:,:,2);
    Us(:,:,1)=Us(:,:,2);
    Vs(:,:,1)=Vs(:,:,2);
    
    Epn2=Epn1;
    Epn1=Epn0;
    Fpn2=Fpn1;
    Fpn1=Fpn0;
    Gpn2=Gpn1;
    Gpn1=Gpn0;
    Fstar2=Fstar1;
    Fstar1=Fstar0;
    Gstar2=Gstar1;
    Gstar1=Gstar0;
    
    eta_NSW(:,:,1)=eta_NSW(:,:,2);
    P_NSW(:,:,1)=P_NSW(:,:,2);
    Q_NSW(:,:,1)=Q_NSW(:,:,2);
    Epn2_NSW=Epn1_NSW;
    Epn1_NSW=Epn0_NSW;
    Fpn2_NSW=Fpn1_NSW;
    Fpn1_NSW=Fpn0_NSW;
    Gpn2_NSW=Gpn1_NSW;
    Gpn1_NSW=Gpn0_NSW;
    
    % Check for instability
    max_eta=max(max(eta(:,:,2)));
    instability_threshold=0.15;
    if max_eta> instability_threshold
        error(['ERROR: SIMULATION GONE WILD @ t=', num2str(t(n))])
    end
    
end







            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

