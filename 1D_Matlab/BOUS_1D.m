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
dy=10;
endx=30;
endy=10;
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

eta_NSW=zeros(nx,ny,2);  % initialize the eta, u, and v matrices (see CP#1 for explanation of this step)
P_NSW=zeros(nx,ny,2);
Q_NSW=zeros(nx,ny,2);

% Initial condition for eta, u and v are zero everywhere for the initial condition
amp=0.15;
for i=1:nx
    for j=1:ny
        if i<(nx/2)
            loc1=sqrt((x(i)-endx/4)^2); %+(y(j)-endy/2)^2);
            eta(i,j,1)=.2*exp(-1.0*loc1^2); % sin IC: sin(x(i)*pi/2); % amp*sin(x(i)*pi/2); %
        else
            loc2=sqrt((x(i)-3*endx/4)^2); %+(y(j)-endy/2)^2);
            eta(i,j,1)=.07*exp(-1.0*loc2^2); % sin IC: sin(x(i)*pi/2); % amp*sin(x(i)*pi/2); %
        end
        h(i,j) = 1+.01*sin(x(i));
%         if x(i) <= endx/2
%             h(i,j)=.2+x(i)*0.025;  % constant water depth
%         elseif x(i) >= endx/2
%             h(i,j)=.2+(endx-x(i))*0.025;
%         end
    end
end

clf
plot(x,-h)
hold on
title('ICs / BCs')
plot(x,eta(:,:,1))
legend('eta', 'depth')
%%

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

Epn0_NSW=zeros(nx,ny);
Epn1_NSW=zeros(nx,ny);
Epn2_NSW=zeros(nx,ny);
Fpn0_NSW=zeros(nx,ny);
Fpn1_NSW=zeros(nx,ny);
Fpn2_NSW=zeros(nx,ny);


% Tri-diagonal coef's
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


for n=1:nt-1  % start the time marching
    ['Time step ', num2str(n),' of ', num2str(nt)]  % plot some information to the command window
    
    % Predictor Step
    [Epn0,Fpn0,Gpn0] = BOUS_EFGcalc(nx,ny,h,eta(:,:,1),P(:,:,1),Q(:,:,1),dx,dy);
    
    [Epn0_NSW,Fpn0_NSW,Gpn0_NSW] = NLSW_EFGcalc(nx,ny,h,eta(:,:,1),P(:,:,1),Q(:,:,1),dx,dy);
    
    plot(x,Us(:,1,2))
    
    if n<3
        for i=1:nx
            for j=1:ny
                eta(i,j,2)=eta(i,j,1)+Epn0(i,j)*dt;
                Us(i,j,2)=Us(i,j,1)+Fpn0(i,j)*dt;
                
                eta_NSW(i,j,2)=eta(i,j,1)+Epn0_NSW(i,j)*dt;
                P_NSW(i,j,2)=P(i,j,1)+Fpn0_NSW(i,j)*dt;
                
            end
        end
    else
        for i=1:nx
            for j=1:ny
                eta(i,j,2)=eta(i,j,1)+dt/12*(23*Epn0(i,j)-16*Epn1(i,j)+5*Epn2(i,j));
                Us(i,j,2)=Us(i,j,1)+dt/12*(23*Fpn0(i,j)-16*Fpn1(i,j)+5*Fpn2(i,j));
                
                eta_NSW(i,j,2)=eta(i,j,1)+dt/12*(23*Epn0_NSW(i,j)-16*Epn1_NSW(i,j)+5*Epn2_NSW(i,j));
                P_NSW(i,j,2)=P(i,j,1)+dt/12*(23*Fpn0_NSW(i,j)-16*Fpn1_NSW(i,j)+5*Fpn2_NSW(i,j));
            end
        end
    end
    
    % Tri-diag solver to get P from Us
    % Correct for boundaries
    Us(1,:,2)=-Us(2,:,2);
    Us(nx,:,2)=-Us(nx-1,:,2);
    
   for j=1:ny
       P(:,j,2)=tridiag(Bx(:,j),Ax(:,j),Cx(:,j), Us(:,j,2));
   end
    
    % enforce bc's
    eta(1,:,2)=eta(2,:,2);
    P(1,:,2)=-P(2,:,2);
    eta(nx,:,2)=eta(nx-1,:,2);
    P(nx,:,2)=-P(nx-1,:,2);
    
    eta_NSW(1,:,2)=eta_NSW(2,:,2);
    P_NSW(1,:,2)=-P_NSW(2,:,2);
    eta_NSW(nx,:,2)=eta_NSW(nx-1,:,2);
    P_NSW(nx,:,2)=-P_NSW(nx-1,:,2);
   

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
    subplot(2,1,1)
    plot(x,eta(:,1,2)',x,eta_NSW(:,1,2)')
    axis([0 endx -amp/2 1.2*amp])
    
    subplot(2,1,2)
    plot(x,P(:,1,2)')
    hold on
    plot(x,10*(P(:,1,2)'-P_NSW(:,1,2)'),'--')
    legend('Bous P','10*(Bous P - NSW P)')
    axis([0 endx -amp/ho*sqrt(9.81*ho) amp/ho*sqrt(9.81*ho)])
    
    pause(.001)
    
    %%%%%% SAVE TO CSV %%%%%%%%%%%%%%
    
    % Boussinesq
    % Difference between P NSW and Bouss, difference between eta NSW and
    % bouss
    if n==1
        xmat = x';
        etamat = eta_NSW(:,1,2)';
        Pmat = P_NSW(:,1,2)';
        Qmat = Q_NSW(:,1,2)';
        hmat = h(:,1)';
        umat = Us(:,1,2)';
        vmat = Vs(:,1,2)';

        Pdiff = P(:,1,2)'-P_NSW(:,1,2)';
        etadiff = eta(:,1,2)'-eta_NSW(:,1,2)';
        Qdiff = Q(:,1,2)'-Q_NSW(:,1,2)';
    else
        xmat = vertcat(xmat,x');
        etamat = vertcat(etamat,eta_NSW(:,1,2)');
        Pmat = vertcat(Pmat,P_NSW(:,1,2)');
        Qmat = vertcat(Qmat,Q_NSW(:,1,2)');
        hmat = vertcat(hmat,h(:,1)');
        umat = vertcat(umat,Us(:,1,2)');
        vmat = vertcat(vmat,Vs(:,1,2)');

        Pdiff = vertcat(Pdiff,P(:,1,2)'-P_NSW(:,1,2)');
        etadiff = vertcat(etadiff,eta(:,1,2)'-eta_NSW(:,1,2)');
        Qdiff = vertcat(Qdiff,Q(:,1,2)'-Q_NSW(:,1,2)');
    end
        
    
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
%     clf
%     subplot(2,1,1)
%     
%     surf(x,y,eta(:,:,1)')  %surface plot of the initial free surface
%     caxis([-.1 .2]);  % set the color axis limits
%     axis([0 endx 0 endy -0.1 .1])  % axis limits of plot
%     view(-30,60)% rotate the surface so that we are looking from an angle, note that view(0,90) is the plan view
%     shading faceted  %  modify the surface shadings
    
    
%     plot(x,eta(:,1,2)')
%     hold on
%     plot(x,eta_NSW(:,1,2)','--')
%     axis([0 endx -amp/2 1.2*amp])
%     
%     subplot(2,1,2)
%     plot(x,P(:,1,2)')
%     hold on
%     plot(x,10*(P(:,1,2)'-P_NSW(:,1,2)'),'--')
%     legend('Bous P','10*(Bous P - NSW P)')
%     axis([0 endx -amp/ho*sqrt(9.81*ho) amp/ho*sqrt(9.81*ho)])
%     
%     pause(.001)
    
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
    
    eta_NSW(:,:,1)=eta_NSW(:,:,2);
    P_NSW(:,:,1)=P_NSW(:,:,2);
    Epn2_NSW=Epn1_NSW;
    Epn1_NSW=Epn0_NSW;
    Fpn2_NSW=Fpn1_NSW;
    Fpn1_NSW=Fpn0_NSW;
    
    % Check for instability
    max_eta=max(max(eta(:,:,2)));
    instability_threshold=0.2;
    if max_eta> instability_threshold
        % Save to csvs
        writematrix(xmat, "xmat11.csv")
        writematrix(hmat, "hmat1.csv")
        writematrix(etamat,"etamat11.csv")
        writematrix(Qmat,"Qmat11.csv")
        writematrix(Pmat,"Pmat11.csv")
        writematrix(umat,"umat11.csv")
        writematrix(vmat,"vmat11.csv")
        writematrix(Pdiff,"Pdiff11.csv")
        writematrix(etadiff,"etadiff11.csv")
        error(['ERROR: SIMULATION GONE WILD @ t=', num2str(t(n))])
    end
    
end

% Save to csvs
writematrix(xmat, "xmat11.csv")
writematrix(hmat, "hmat11.csv")
writematrix(etamat,"etamat11.csv")
writematrix(Qmat,"Qmat11.csv")
writematrix(Pmat,"Pmat11.csv")
writematrix(umat,"umat11.csv")
writematrix(vmat,"vmat11.csv")
writematrix(Pdiff,"Pdiff11.csv")
writematrix(etadiff,"etadiff11.csv")