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
h=1.;
dx=0.1*h;
dy=10;
endx=30;
endy=10;
endt=20;
courant=0.25;

clf  % clear the current figure

x=[0:dx:endx]';  % create x column vector
y=[0:dy:endy]';  % create y column vector

dt=courant*dx/sqrt(9.81*h);  % calculate time step
t=[0:dt:endt]';  % create t column vector

nx=length(x);  % total number of grid points in x
ny=length(y);  % total number of grid points in y
nt=length(t);  % total number of time steps

eta=zeros(nx,ny,2);  % initialize the eta, u, and v matrices (see CP#1 for explanation of this step)
P=zeros(nx,ny,2);
Q=zeros(nx,ny,2);


% Initial condition for eta, u and v are zero everywhere for the initial condition
amp=0.15;
for i=1:nx
    for j=1:ny
        loc=sqrt((x(i)-endx/2)^2); %+(y(j)-endy/2)^2);
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
Epn0=zeros(nx,ny);
Epn1=zeros(nx,ny);
Epn2=zeros(nx,ny);
Fpn0=zeros(nx,ny);
Fpn1=zeros(nx,ny);
Fpn2=zeros(nx,ny);


for n=1:nt-1  % start the time marching
    ['Time step ', num2str(n),' of ', num2str(nt)]  % plot some information to the command window
    
    % Predictor Step
    [Epn0,Fpn0,Gpn0] = NLSW_EFGcalc(nx,ny,h,eta(:,:,1),P(:,:,1),Q(:,:,1),dx,dy);
    
    if n<3
        for i=1:nx
            for j=1:ny
                eta(i,j,2)=eta(i,j,1)+Epn0(i,j)*dt;
                P(i,j,2)=P(i,j,1)+Fpn0(i,j)*dt;
            end
        end
    else
        for i=1:nx
            for j=1:ny
                eta(i,j,2)=eta(i,j,1)+dt/12*(23*Epn0(i,j)-16*Epn1(i,j)+5*Epn2(i,j));
                P(i,j,2)=P(i,j,1)+dt/12*(23*Fpn0(i,j)-16*Fpn1(i,j)+5*Fpn2(i,j));
            end
        end
    end
        
    
    % enforce bc's
    eta(1,:,2)=eta(2,:,2);
    P(1,:,2)=-P(2,:,2);
    eta(nx,:,2)=eta(nx-1,:,2);
    P(nx,:,2)=-P(nx-1,:,2);

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
    plot(x,eta(:,:,2)')
    % caxis([-.1 .1]);
    axis([0 endx -amp/2 1.2*amp])
    % view(-30,60)
    % shading faceted
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
    Epn2=Epn1;
    Epn1=Epn0;
    Fpn2=Fpn1;
    Fpn1=Fpn0;
    
    % Check for instability
    max_eta=max(max(eta(:,:,2)));
    instability_threshold=0.15;
    if max_eta> instability_threshold
        error(['ERROR: SIMULATION GONE WILD @ t=', num2str(t(n))])
    end
    
end







            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

