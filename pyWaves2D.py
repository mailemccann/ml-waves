import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

'''
This script contains functions that solve the nonlinear 
shallow water equations using Euler's method.  

As set-up, 
it gives no output, although it creates an animation plot
of the free surface.

INPUT
    dx - the grid size in the x-direction, in meters
    dy - the grid size in the y-direction, in meters
    endx - the total domain length in the x-direction, in meters
    endy - the total domain length in the y-direction, in meters
    endt - the total time to run the simulation, in seconds
    h - the constant water depth, in meters
    courant - the courant number, controls the time step
As an example call, use:  >>CP2_Heun_iter(1,1,50,50,100,1,0.5)
'''
def solveNLSW_BOUS(disp_plot=True): #dx,dy,endx,endy,endt,h,courant
    # Set plot = False if you don't want to display the surface plot

    ho = 1. 
    dx = 0.1*ho 
    dy = 0.1*ho 
    endx = 20 
    endy = 20 
    endt = 20 
    courant = 0.25 

    x = np.arange(0, endx+dx, dx).astype(np.double)[:,np.newaxis]   # create x column vector
    y = np.arange(0, endy+dy, dy).astype(np.double)[:,np.newaxis]  # create y column vector
    dt = courant*dx/np.sqrt(9.81*ho)   # calculate time step
    t  = np.arange(0, endt+dt, dt).astype(np.double)[:,np.newaxis]  # create t column vector

    nx = len(x)   # total number of grid points in x
    ny = len(y)   # total number of grid points in y
    nt = len(t)   # total number of time steps

    eta = np.zeros((nx,ny,2))
    P = np.zeros((nx,ny,2))
    Q = np.zeros((nx,ny,2))
    Us = np.zeros((nx,ny,2))
    Vs = np.zeros((nx,ny,2))   # initialize the eta, u, and v matrices (see CP#1 for explanation of this step)
    h = np.ones((nx,ny))*ho

    eta_NSW = np.zeros((nx,ny,2))
    P_NSW = np.zeros((nx,ny,2))
    Q_NSW = np.zeros((nx,ny,2))   # initialize the eta, u, and v matrices (see CP#1 for explanation of this step)

    # Initial condition for eta, u and v are zero everywhere for the initial condition
    amp = 0.15
    for ii in range(nx):
        for jj in range(ny):
            loc = np.sqrt((x[ii]-endx/2)**2+(y[jj]-endy/2)**2)
            eta[ii,jj,0] = amp*np.exp(-1.0*loc**2)
            eta_NSW[ii,jj,0] = amp*np.exp(-1.0*loc**2)

    # Init Bouss variables
    Epn0=np.zeros((nx,ny)) 
    Epn1=np.zeros((nx,ny)) 
    Epn2=np.zeros((nx,ny)) 

    Fpn0=np.zeros((nx,ny)) 
    Fpn1=np.zeros((nx,ny)) 
    Fpn2=np.zeros((nx,ny)) 

    Gpn0=np.zeros((nx,ny)) 
    Gpn1=np.zeros((nx,ny)) 
    Gpn2=np.zeros((nx,ny)) 

    Fstar0=np.zeros((nx,ny)) 
    Fstar1=np.zeros((nx,ny)) 
    Fstar2=np.zeros((nx,ny)) 

    Gstar0=np.zeros((nx,ny)) 
    Gstar1=np.zeros((nx,ny)) 
    Gstar2=np.zeros((nx,ny)) 

    # Init NSW variables
    Epn0_NSW=np.zeros((nx,ny)) 
    Epn1_NSW=np.zeros((nx,ny)) 
    Epn2_NSW=np.zeros((nx,ny)) 

    Fpn0_NSW=np.zeros((nx,ny)) 
    Fpn1_NSW=np.zeros((nx,ny)) 
    Fpn2_NSW=np.zeros((nx,ny)) 

    Gpn0_NSW=np.zeros((nx,ny)) 
    Gpn1_NSW=np.zeros((nx,ny)) 
    Gpn2_NSW=np.zeros((nx,ny)) 

    # Tri-diagonal coef's, (note that these need to be updated for variable bathymetry)
    B=1/15 
    Ax= -1*(B+1/3)*np.square(h)/dx**2 
    Bx=1+2*(B+1/3)*np.square(h)/dx**2 
    Cx= -1*(B+1/3)*np.square(h)/dx**2 

    # Correct for boundaries
    Ax[0,:]=0 
    Bx[0,:]=1 
    Cx[0,:]=0 

    Ax[-1,:]=0 
    Bx[-1,:]=1 
    Cx[-1,:]=0 

    Ay= -1*(B+1/3)*np.square(h)/dy**2 
    By=1+2*(B+1/3)*np.square(h)/dy**2 
    Cy= -1*(B+1/3)*np.square(h)/dy**2 

    # Correct for boundaries
    Ay[:,0]=0 
    By[:,0]=1 
    Cy[:,0]=0 

    Ay[:,-1]=0 
    By[:,-1]=1 
    Cy[:,-1]=0

    if disp_plot:
        fig, ax = plt.subplots(1,3,subplot_kw={"projection": "3d"})

    for n in range(nt):  # start the time marching
        print('Time step ' + str(n) + ' of ' + str(nt))  # plot some information to the command window

        # Predictor Step
        [Epn0,Fpn0,Gpn0,Fstar0,Gstar0] = BOUS_EFGcalc(nx,ny,h,eta[:,:,0],P[:,:,0],Q[:,:,0],dx,dy) 
        
        [Epn0_NSW,Fpn0_NSW,Gpn0_NSW] = NLSW_EFGcalc(nx,ny,h,eta[:,:,0],P[:,:,0],Q[:,:,0],dx,dy) 

        if n<2:
            for i in range(nx):
                for j in range(ny):
                    eta[i,j,1]=eta[i,j,0]+Epn0[i,j]*dt 
                    Us[i,j,1]=Us[i,j,0]+Fpn0[i,j]*dt+Fstar0[i,j]-Fstar1[i,j] 
                    Vs[i,j,1]=Vs[i,j,0]+Gpn0[i,j]*dt+Gstar0[i,j]-Gstar1[i,j] 
                    
                    eta_NSW[i,j,1]=eta[i,j,0]+Epn0_NSW[i,j]*dt 
                    P_NSW[i,j,1]=P[i,j,0]+Fpn0_NSW[i,j]*dt 
                    Q_NSW[i,j,1]=Q[i,j,0]+Gpn0_NSW[i,j]*dt 
        else:
            for i in range(nx):
                for j in range(ny):
                    eta[i,j,1]=eta[i,j,0]+dt/12*(23*Epn0[i,j]-16*Epn1[i,j]+5*Epn2[i,j]) 
                    Us[i,j,1] =Us[i,j,0] +dt/12*(23*Fpn0[i,j]-16*Fpn1[i,j]+5*Fpn2[i,j])+Fstar0[i,j]-Fstar1[i,j] 
                    Vs[i,j,1] =Vs[i,j,0] +dt/12*(23*Gpn0[i,j]-16*Gpn1[i,j]+5*Gpn2[i,j])+Gstar0[i,j]-Gstar1[i,j] 
                    
                    eta_NSW[i,j,1]=eta[i,j,0]+dt/12*(23*Epn0_NSW[i,j]-16*Epn1_NSW[i,j]+5*Epn2_NSW[i,j]) 
                    P_NSW[i,j,1]=P[i,j,0]+dt/12*(23*Fpn0_NSW[i,j]-16*Fpn1_NSW[i,j]+5*Fpn2_NSW[i,j]) 
                    Q_NSW[i,j,1]=Q[i,j,0]+dt/12*(23*Gpn0_NSW[i,j]-16*Gpn1_NSW[i,j]+5*Gpn2_NSW[i,j]) 

        # Bous BC's
        # left wall
        eta[0,:,1] = eta[1,:,1]
        Us[0,:,1] = -Us[1,:,1]
        Vs[0,:,1] = Vs[1,:,1]
        
        #right wall
        eta[-1,:,1] = eta[-2,:,1] 
        Us[-1,:,1] = -Us[-2,:,1]  
        Vs[-1,:,1] = Vs[-2,:,1] 
        
        #bottom wall
        eta[:,0,1] = eta[:,1,1] 
        Us[:,0,1]  = Us[:,1,1] 
        Vs[:,0,1] = -Vs[:,1,1] 
        
        #top wall
        eta[:,-1,1] = eta[:,-2,1]
        Us[:,-1,1] = Us[:,-2,1] 
        Vs[:,-1,1] = -Vs[:,-2,1]

        # Tri-diag solver to get P from Us 
        for jj in range(ny):
            P[:,jj,1]=tridiag(Bx[:,jj],Ax[:,jj],Cx[:,jj], Us[:,jj,1]) 
         
        
        # Tri-diag solver to get Q from Vs
        for ii in range(nx):
            Q[ii,:,1]=tridiag(By[ii,:],Ay[ii,:],Cy[ii,:], Vs[ii,:,1]) 

        # NSW BC's
        # left wall
        eta_NSW[0,:,1] = eta_NSW[1,:,1]
        P_NSW[0,:,1] = -P_NSW[1,:,1]
        Q_NSW[0,:,1] = Q_NSW[1,:,1]
        
        #right wall
        eta_NSW[-1,:,1] = eta_NSW[-2,:,1]
        P_NSW[-1,:,1] = -P_NSW[-2,:,1]
        Q_NSW[-1,:,1] = Q_NSW[-2,:,1]
        
        #bottom wall
        eta_NSW[:,0,1] = eta_NSW[:,1,1]
        P_NSW[:,0,1] = P_NSW[:,1,1]
        Q_NSW[:,0,1] = -Q_NSW[:,1,1]
        
        #top wall
        eta_NSW[:,-1,1] = eta_NSW[:,-2,1]
        P_NSW[:,-1,1] = P_NSW[:,-2,1]
        Q_NSW[:,-1,1] = -Q_NSW[:,-2,1]

        if disp_plot:
            #if (n*dt)%1 == 0:
            # Plot 
            ax[0].cla()
            #ax.plot(x[1:-2],np.asarray(100*(Q[0,1:-2,1]-Q_NSW[0,1:-2,1])),label = '100*(Bous Q - NSW Q)')
            X, Y = np.meshgrid(x, y)
            Z = eta_NSW[:,:,0] #P[1:(nx-2),1:(ny-2),1] #-P_NSW[1:(nx-2),1:(ny-2),1])*100
            surf0 = ax[0].plot_surface(X, Y, Z, linewidth=0, antialiased=False, cmap = 'Blues',vmin=-.1, vmax = .15)
            #ax.plot(x[1:-2],np.asarray(100*(P[1:(nx-2),int(ny/2),1]-P_NSW[1:(nx-2),int(ny/2),1])),label = '100*(Bous P - NSW P)')
            #ax.legend()
            #fig.colorbar(surf, shrink=0.5, aspect=5)
            ax[0].set_xlabel('x (m)')
            ax[0].set_ylabel('y (m)')
            ax[0].set_zlim([-.5,.2])
            ax[0].set_title(r'$\eta NLSW')

            ax[1].cla()
            X, Y = np.meshgrid(x, y)
            Z = eta[:,:,1] #P[1:(nx-2),1:(ny-2),1] #-P_NSW[1:(nx-2),1:(ny-2),1])*100
            surf1 = ax[1].plot_surface(X, Y, Z, linewidth=0, antialiased=False, cmap = 'Blues', vmin=-.1, vmax = .15)
            #fig.colorbar(surf, shrink=0.5, aspect=5)
            ax[1].set_xlabel('x (m)')
            ax[1].set_ylabel('y (m)')
            ax[1].set_xlabel('x (m)')
            ax[1].set_ylabel('y (m)')
            ax[1].set_zlim([-.5,.2])
            ax[1].set_title(r'$\eta Boussinesq')

            ax[2].cla()
            #X, Y = np.meshgrid(x, y)
            Z = (P[:,:,1] -P_NSW[:,:,1])*100 #eta[:,:,1] #P[1:(nx-2),1:(ny-2),1] #-P_NSW[1:(nx-2),1:(ny-2),1])*100
            surf2 = ax[2].plot_surface(X,Y,Z, linewidth=0, antialiased=False, cmap = 'RdBu')
            #fig.colorbar(surf, shrink=0.5, aspect=5)
            ax[2].set_xlabel('x (m)')
            ax[2].set_ylabel('y (m)')
            ax[2].set_title(r'$\delta P_Bouss - P_NLSW *100')

            fig.suptitle('time = ' + str(np.round(n*dt,3)) + ' s')
            plt.pause(0.05)

        # Shift all data one time level back.  This is a trick that allows use to use much smaller matrix
        # sizes for these variables.  We are calculating these values at all x,y and t points, so, the "total"
        # calculated domains are nx by ny by nt.  However, since we only need two time levels to perform
        # the Eulers method, we only keep two time levels in memory.  Memory savings here are great.  For example, if
        # we to calculate 5000 time steps on grid with nx=100 and ny=100, we save 1.2 GB (gigabytes) of memory space
        # by only keeping two levels in memory.

        eta[:,:,0]=eta[:,:,1]
        P[:,:,0]=P[:,:,1] 
        Q[:,:,0]=Q[:,:,1]
        Us[:,:,0]=Us[:,:,1]
        Vs[:,:,0]=Vs[:,:,1]
        
        Epn2=deepcopy(Epn1)
        Epn1=deepcopy(Epn0)
        Fpn2=deepcopy(Fpn1)
        Fpn1=deepcopy(Fpn0) 
        Gpn2=deepcopy(Gpn1) 
        Gpn1=deepcopy(Gpn0) 
        Fstar2=deepcopy(Fstar1)
        Fstar1=deepcopy(Fstar0) 
        Gstar2=deepcopy(Gstar1) 
        Gstar1=deepcopy(Gstar0) 
        
        eta_NSW[:,:,0]=eta_NSW[:,:,1]
        P_NSW[:,:,0]=P_NSW[:,:,1]
        Q_NSW[:,:,0]=Q_NSW[:,:,1]
        Epn2_NSW=deepcopy(Epn1_NSW) 
        Epn1_NSW=deepcopy(Epn0_NSW)
        Fpn2_NSW=deepcopy(Fpn1_NSW)
        Fpn1_NSW=deepcopy(Fpn0_NSW)
        Gpn2_NSW=deepcopy(Gpn1_NSW) 
        Gpn1_NSW=deepcopy(Gpn0_NSW) 

        # Check for instability
        '''
        max_eta=max(max(eta[:,:,2])) 
        instability_threshold=0.15 
        if max_eta> instability_threshold:
            error(['ERROR: SIMULATION GONE WILD @ t=', num2str(t(n))])
        '''

    cv2.destroyAllWindows()
    video.release()


def BOUS_EFGcalc(nx,ny,h,eta,P,Q,dx,dy):
    # This function file will calculate the right hand sides of the continuity (E), x-momentum (F), and y-momentum (G) equations.
    # You should not have to change this file for CP#2
    # INPUT
    #		nx - # of points in x direction
    #		ny - # of points in y direction
    #		h - water depth nx by ny
    #		eta - eta (free surface elevation) surface matrix size nx by ny
    #		u - u (x-velocity) surface matrix size nx by ny
    #		v - v (y-velocity) surface matrix size nx by ny
    #		dx - grid size in x-direction
    #		dy - grid size in y-direction
    #
    # OUTPUT
    #		E,F,G,F*,G* - size nx by ny

    E=np.zeros((nx,ny)) 
    F=np.zeros((nx,ny)) 
    G=np.zeros((nx,ny))
    Fstar=np.zeros((nx,ny)) 
    Gstar=np.zeros((nx,ny)) 

    B=1/15 

    # group variables into temp matrices
    H=h+eta 
    Pflux_F = np.square(P)/H+9.81*np.square(H)/2 
    Qflux_G = np.square(Q)/H+9.81*np.square(H)/2

    Pflux_G=P*Q/H
    Qflux_F=P*Q/H

    # precalc some derivatives
    P_x=np.zeros((nx,ny)) 
    Q_y=np.zeros((nx,ny)) 
    eta_x=np.zeros((nx,ny)) 
    eta_y=np.zeros((nx,ny)) 
    eta_xx=np.zeros((nx,ny)) 
    eta_yy=np.zeros((nx,ny)) 

    for i in range(nx):
        for j in range(ny):
        # evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
            if i==0:
                P_x[i,j]=(P[i+1,j]-P[i,j])/dx 
                eta_x[i,j]=(eta[i+1,j]-eta[i,j])/dx 
                eta_xx[i,j]=0 
            elif i==nx-1:
                P_x[i,j]=(P[i,j]-P[i-1,j])/dx 
                eta_x[i,j]=(eta[i,j]-eta[i-1,j])/dx 
                eta_xx[i,j]=0 
            else: 
                P_x[i,j]=(P[i+1,j]-P[i-1,j])/dx/2 
                eta_x[i,j]=(eta[i+1,j]-eta[i-1,j])/dx/2 
                eta_xx[i,j]=(eta[i+1,j]-2*eta[i,j]+eta[i-1,j])/(dx**2)
            
            if j==0:
                Q_y[i,j]=(Q[i,j+1]-Q[i,j])/dy 
                eta_y[i,j]=(eta[i,j+1]-eta[i,j])/dy 
                eta_yy[i,j]=0 
            elif j==ny-1:
                Q_y[i,j]=(Q[i,j]-Q[i,j-1])/dy 
                eta_y[i,j]=(eta[i,j]-eta[i,j-1])/dy 
                eta_yy[i,j]=0 
            else:
                Q_y[i,j]=(Q[i,j+1]-Q[i,j-1])/dy/2 
                eta_y[i,j]=(eta[i,j+1]-eta[i,j-1])/dy/2 
                eta_yy[i,j]=(eta[i,j+1]-2*eta[i,j]+eta[i,j-1])/(dy**2)

    for i in range(nx):
        for j in range(ny):
            # evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
            if i==0:
                Pflux_x=(Pflux_F[i+1,j]-Pflux_F[i,j])/dx 
                Qflux_x=(Qflux_F[i+1,j]-Qflux_F[i,j])/dx 
                Q_xy=       (Q_y[i+1,j]    -Q_y[i,j])/dx 
                eta_xxx= (eta_xx[i+1,j] -eta_xx[i,j])/dx 
                eta_xyy= (eta_yy[i+1,j] -eta_yy[i,j])/dx 
            elif i==nx-1:
                Pflux_x=(Pflux_F[i,j]-Pflux_F[i-1,j])/dx 
                Qflux_x=(Qflux_F[i,j]-Qflux_F[i-1,j])/dx 
                Q_xy=       (Q_y[i,j]    -Q_y[i-1,j])/dx 
                eta_xxx= (eta_xx[i,j] -eta_xx[i-1,j])/dx 
                eta_xyy= (eta_yy[i,j] -eta_yy[i-1,j])/dx 
            else:
                Pflux_x=(Pflux_F[i+1,j]-Pflux_F[i-1,j])/dx/2 
                Qflux_x=(Qflux_F[i+1,j]-Qflux_F[i-1,j])/dx/2 
                Q_xy=       (Q_y[i+1,j]    -Q_y[i-1,j])/dx/2
                eta_xxx= (eta_xx[i+1,j] -eta_xx[i-1,j])/dx/2
                eta_xyy= (eta_yy[i+1,j] -eta_yy[i-1,j])/dx/2

            
            if j==0:
                Pflux_y=(Pflux_G[i,j+1]-Pflux_G[i,j])/dy 
                Qflux_y=(Qflux_G[i,j+1]-Qflux_G[i,j])/dy 
                P_xy=       (P_x[i,j+1]    -P_x[i,j])/dy 
                eta_yyy= (eta_yy[i,j+1] -eta_yy[i,j])/dy 
                eta_xxy= (eta_xx[i,j+1] -eta_xx[i,j])/dy 
            elif j==ny-1:
                Pflux_y=(Pflux_G[i,j]-Pflux_G[i,j-1])/dy 
                Qflux_y=(Qflux_G[i,j]-Qflux_G[i,j-1])/dy 
                P_xy=       (P_x[i,j]    -P_x[i,j-1])/dy 
                eta_yyy= (eta_yy[i,j] -eta_yy[i,j-1])/dy 
                eta_xxy= (eta_xx[i,j] -eta_xx[i,j-1])/dy 
            else:
                Pflux_y=(Pflux_G[i,j+1]-Pflux_G[i,j-1])/dy/2 
                Qflux_y=(Qflux_G[i,j+1]-Qflux_G[i,j-1])/dy/2 
                P_xy=       (P_x[i,j+1]    -P_x[i,j-1])/dy/2
                eta_yyy= (eta_yy[i,j+1] -eta_yy[i,j-1])/dy/2
                eta_xxy= (eta_xx[i,j+1] -eta_xx[i,j-1])/dy/2


            E[i,j]=-P_x[i,j]-Q_y[i,j]  
            F[i,j]=-Pflux_x - Pflux_y + B*9.81*(h[i,j]**3)*(eta_xxx+eta_xyy)    # (note that these need to be updated for variable bathymetry)
            G[i,j]=-Qflux_x - Qflux_y + B*9.81*(h[i,j]**3)*(eta_yyy+eta_xxy)     # (note that these need to be updated for variable bathymetry)
            
            Fstar[i,j]=(B+1/3)*(h[i,j]**2)*Q_xy   #(note that these need to be updated for variable bathymetry)
            Gstar[i,j]=(B+1/3)*(h[i,j]**2)*P_xy   # (note that these need to be updated for variable bathymetry)    

    return [E,F,G,Fstar,Gstar]
            
def NLSW_EFGcalc(nx,ny,h,eta,P,Q,dx,dy):
    # This function file will calculate the right hand sides of the continuity (E), x-momentum (F), and y-momentum (G) equations.
    # INPUT
    #		nx - # of points in x direction
    #		ny - # of points in y direction
    #		h - constant water depth
    #		eta - eta (free surface elevation) surface matrix size nx by ny
    #		u - u (x-velocity) surface matrix size nx by ny
    #		v - v (y-velocity) surface matrix size nx by ny
    #		dx - grid size in x-direction
    #		dy - grid size in y-direction
    #
    # OUTPUT
    #		E,F,G - size nx by ny

    E=np.zeros((nx,ny)) 
    F=np.zeros((nx,ny)) 
    G=np.zeros((nx,ny)) 

    # group variables into temp matrices
    H=h+eta 
    Pflux_F = np.square(P)/H + 9.81*np.square(H)/2 
    Pflux_G = P*Q/H

    Qflux_F = P*Q/H
    Qflux_G = np.square(Q)/H + 9.81*np.square(H)/2 

    for i in range(nx):
        for j in range(ny):
            # evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
            eta_x=0
            eta_y=0
            if i==0:
                P_x=(P[i+1,j]-P[i,j])/dx 
                Pflux_x=(Pflux_F[i+1,j]-Pflux_F[i,j])/dx 
                Qflux_x=(Qflux_F[i+1,j]-Qflux_F[i,j])/dx 
            elif i==nx-1:
                P_x=(P[i,j]-P[i-1,j])/dx 
                Pflux_x=(Pflux_F[i,j]-Pflux_F[i-1,j])/dx 
                Qflux_x=(Qflux_F[i,j]-Qflux_F[i-1,j])/dx 
            else:
                eta_x = (eta[i+1,j]-eta[i-1,j])/dx/2
                P_x=(P[i+1,j]-P[i-1,j])/dx/2 
                Pflux_x=(Pflux_F[i+1,j]-Pflux_F[i-1,j])/dx/2 
                Qflux_x=(Qflux_F[i+1,j]-Qflux_F[i-1,j])/dx/2 

            if j==0:
                Q_y=(Q[i,j+1]-Q[i,j])/dx 
                Pflux_y=(Pflux_G[i,j+1]-Pflux_G[i,j])/dx 
                Qflux_y=(Qflux_G[i,j+1]-Qflux_G[i,j])/dx 
            elif j==ny-1:
                Q_y=(Q[i,j]-Q[i,j-1])/dx 
                Pflux_y=(Pflux_G[i,j]-Pflux_G[i,j-1])/dx 
                Qflux_y=(Qflux_G[i,j]-Qflux_G[i,j-1])/dx 
            else:
                eta_y = (eta[i,j+1]-eta[i,j-1])/dx/2
                Q_y = (Q[i,j+1]-Q[i,j-1])/dx/2 
                Pflux_y=(Pflux_G[i,j+1]-Pflux_G[i,j-1])/dx/2 
                Qflux_y=(Qflux_G[i,j+1]-Qflux_G[i,j-1])/dx/2 

        E[i,j] = -P_x - Q_y  
        F[i,j] = -Pflux_x - Pflux_y
        G[i,j] = -Qflux_x - Qflux_y  

    return [E,F,G]


def tridiag( a, b, c, f ):
    #  Solve the  n x n  tridiagonal system for y:
    #
    #  [ a(1)  c(1)                                  ] [  y(1)  ]   [  f(1)  ]
    #  [ b(2)  a(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
    #  [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
    #  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
    #  [                    ...    ...    ...        ] [        ]   [        ]
    #  [                        b(n-1) a(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
    #  [                                 b(n)  a(n)  ] [  y(n)  ]   [  f(n)  ]
    #
    #  f must be a vector (row or column) of length n
    #  a, b, c must be vectors of length n (note that b(1) and c(n) are not used)
    # some additional information is at the end of the file

    n = len(f) 
    v = np.zeros((n))    
    y = np.zeros((n)) 
    w = a[0] 
    y[0] = f[0]/w 

    for i in range(1,n):
        v[i-1] = c[i-1]/w 
        w = a[i] - b[i]*v[i-1] 
        y[i] = ( f[i] - b[i]*y[i-1])/w 

    for j in range(n-2,0,-1):
        y[j] = y[j] - v[j]*y[j+1] 

    return y


if __name__ == '__main__':
    solveNLSW_BOUS(disp_plot=True)
