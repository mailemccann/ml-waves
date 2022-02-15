import numpy as np
import matplotlib.pyplot as plt

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
def solveNLSW_BOUS(): #dx,dy,endx,endy,endt,h,courant
    ho = 1. 
    dx = 0.1*ho 
    dy = 10 
    endx = 30 
    endy = 10 
    endt = 20 
    courant = 0.25 

    x = np.arange(0, endx+dx, dx).astype(np.double)[:,np.newaxis]   # create x column vector
    y = np.arange(0, endy+dy, dy).astype(np.double)[:,np.newaxis]  # create y column vector
    dt = courant*dx/np.sqrt(9.81*ho)   # calculate time step
    t  = np.arange(0, endt+dt, dt).astype(np.double)[:,np.newaxis]  # create t column vector

    nx = len(x)   # total number of grid points in x
    ny = len(y)   # total number of grid points in y
    nt = len(t)   # total number of time steps

    eta = P = Q = Us = Vs = np.zeros((nx,ny,2))   # initialize the eta, u, and v matrices (see CP#1 for explanation of this step)
    h = np.zeros((nx,ny)) 

    eta_NSW = P_NSW = Q_NSW = np.zeros((nx,ny,2))   # initialize the eta, u, and v matrices (see CP#1 for explanation of this step)

    # Initial condition for eta, u and v are zero everywhere for the initial condition
    amp = 0.15
    for i in range(nx):
        for j in range(ny):
            loc = np.sqrt((x[i]-endx/2)**2) #+(y(j)-endy/2)**2)
            eta[i,j,0] = amp*np.exp(-1.0*loc**2)
            h[i,j] = ho  # constant water depth


    eta_previous=eta[:,:,0]
    eta_NSW=eta

    Epn0=np.zeros((nx,ny)) 
    Epn1=np.zeros((nx,ny)) 
    Epn2=np.zeros((nx,ny)) 
    Fpn0=np.zeros((nx,ny)) 
    Fpn1=np.zeros((nx,ny)) 
    Fpn2=np.zeros((nx,ny)) 

    Epn0_NSW=np.zeros((nx,ny)) 
    Epn1_NSW=np.zeros((nx,ny)) 
    Epn2_NSW=np.zeros((nx,ny)) 
    Fpn0_NSW=np.zeros((nx,ny)) 
    Fpn1_NSW=np.zeros((nx,ny)) 
    Fpn2_NSW=np.zeros((nx,ny)) 

    # Tri-diagonal coef's
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

    for n in range(nt):  # start the time marching
        print('Time step ' + str(n) + ' of ' + str(nt))  # plot some information to the command window
        
        # Predictor Step
        [Epn0,Fpn0,Gpn0] = BOUS_EFGcalc(nx,ny,h,eta[:,:,0],P[:,:,0],Q[:,:,0],dx,dy) 
        
        [Epn0_NSW,Fpn0_NSW,Gpn0_NSW] = NLSW_EFGcalc(nx,ny,h,eta[:,:,0],P[:,:,0],Q[:,:,0],dx,dy) 
        
        if n<3:
            for i in range(nx):
                for j in range(ny):
                    eta[i,j,1]=eta[i,j,0]+Epn0[i,j]*dt 
                    Us[i,j,1]=Us[i,j,0]+Fpn0[i,j]*dt 
                    
                    eta_NSW[i,j,1]=eta[i,j,0]+Epn0_NSW[i,j]*dt 
                    P_NSW[i,j,1]=P[i,j,0]+Fpn0_NSW[i,j]*dt 
        else:
            for i in range(nx):
                for j in range(ny):
                    eta[i,j,1]=eta[i,j,0]+dt/12*(23*Epn0[i,j]-16*Epn1[i,j]+5*Epn2[i,j]) 
                    Us[i,j,1]=Us[i,j,0]+dt/12*(23*Fpn0[i,j]-16*Fpn1[i,j]+5*Fpn2[i,j]) 
                    
                    eta_NSW[i,j,1]=eta[i,j,0]+dt/12*(23*Epn0_NSW[i,j]-16*Epn1_NSW[i,j]+5*Epn2_NSW[i,j]) 
                    P_NSW[i,j,1]=P[i,j,0]+dt/12*(23*Fpn0_NSW[i,j]-16*Fpn1_NSW[i,j]+5*Fpn2_NSW[i,j]) 

        
        # Tri-diag solver to get P from Us
        # Correct for boundaries
        Us[0,:,1]=-Us[1,:,1]
        Us[-1,:,1]=-Us[-2,:,1] 
        
        for j in range(ny):
            P[:,j,1]=tridiag(Bx[:,j],Ax[:,j],Cx[:,j], Us[:,j,1]) 
        
        # enforce bc's
        eta[0,:,1]=eta[1,:,1]
        P[0,:,1]=-P[1,:,1]
        eta[-1,:,1]=eta[-2,:,1] 
        P[-1,:,1]=-P[-2,:,1]
        
        eta_NSW[0,:,1]=eta_NSW[1,:,1]
        P_NSW[0,:,1]=-P_NSW[1,:,1]
        eta_NSW[-1,:,1]=eta_NSW[-2,:,1]
        P_NSW[-1,:,1]=-P_NSW[-2,:,1]

        # Plot 
        plt.plot(x,np.asarray(P[:,0,1]), label = 'Bous P')
        plt.plot(x,np.asarray(10*(P[:,0,1]-P_NSW[:,0,1])),'--', label = '10*(Bous P - NSW P)')
        plt.pause(0.05)

        plt.show()


        
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

        Epn2=Epn1 
        Epn1=Epn0 
        Fpn2=Fpn1 
        Fpn1=Fpn0 
        
        eta_NSW[:,:,0]=eta_NSW[:,:,1]
        P_NSW[:,:,0]=P_NSW[:,:,1]
        Epn2_NSW=Epn1_NSW 
        Epn1_NSW=Epn0_NSW 
        Fpn2_NSW=Fpn1_NSW 
        Fpn1_NSW=Fpn0_NSW 
        
        # Check for instability
        '''
        max_eta=max(max(eta[:,:,2])) 
        instability_threshold=0.15 
        if max_eta> instability_threshold:
            error(['ERROR: SIMULATION GONE WILD @ t=', num2str(t(n))])
        '''


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
    #		E,F,G - size nx by ny

    E=np.zeros((nx,ny)) 
    F=np.zeros((nx,ny)) 
    G=np.zeros((nx,ny)) 

    B=1/15 

    # group variables into temp matrices
    H=h+eta 
    Pflux=np.square(P)/H+9.81*np.square(H)/2 

    for i in range(nx):
        for j in range(ny):
            # evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
            if i==0:
                P_x=(P[i+1,j]-P[i,j])/dx 
                Pflux_x=(Pflux[i+1,j]-Pflux[i,j])/dx 
            elif i==nx-1:
                P_x=(P[i,j]-P[i-1,j])/dx 
                Pflux_x=(Pflux[i,j]-Pflux[i-1,j])/dx 
            else:
                P_x=(P[i+1,j]-P[i-1,j])/dx/2 
                Pflux_x=(Pflux[i+1,j]-Pflux[i-1,j])/dx/2 
            
            if i<2:
                eta_xxx=(eta[i+3,j]-3*eta[i+2,j]+3*eta[i+1,j]-eta[i,j])/(dx**3)   # forward
            elif i>nx-3:
                eta_xxx=(eta[i,j]-3*eta[i-1,j]+3*eta[i-2,j]-eta[i-3,j])/(dx**3)    # backward
            else:
                eta_xxx=(eta[i+2,j]-2*eta[i+1,j]+2*eta[i-1,j]-eta[i-2,j])/(2*dx**3)    # centered

            E[i,j]=-P_x  
            F[i,j]=-Pflux_x+B*9.81*h[i,j]**3*eta_xxx 
            G[i,j]=0  #-9.81*eta_y-0.5*u_sq_y-0.5*v_sq_y 

    return [E,F,G]
            
def NLSW_EFGcalc(nx,ny,h,eta,P,Q,dx,dy):
    # This function file will calculate the right hand sides of the continuity (E), x-momentum (F), and y-momentum (G) equations.
    # You should not have to change this file for CP#2
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
    Pflux=np.square(P)/H + 9.81*np.square(H)/2 

    for i in range(nx):
        for j in range(ny):
            # evaluate derivatives, using forward or backward derivatives at the boundaries (vertical walls)
            if i==0:
                P_x=(P[i+1,j]-P[i,j])/dx 
                Pflux_x=(Pflux[i+1,j]-Pflux[i,j])/dx 
            elif i==nx-1:
                P_x=(P[i,j]-P[i-1,j])/dx 
                Pflux_x=(Pflux[i,j]-Pflux[i-1,j])/dx 
            else:
                P_x=(P[i+1,j]-P[i-1,j])/dx/2 
                Pflux_x=(Pflux[i+1,j]-Pflux[i-1,j])/dx/2 

        E[i,j]=-P_x  
        F[i,j]=-Pflux_x 
        G[i,j]=0  #-9.81*eta_y-0.5*u_sq_y-0.5*v_sq_y 

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
    y = v 
    w = a[0] 
    y[0] = f[0]/w 

    for i in range(2,n):
        v[i-1] = c[i-1]/w 
        w = a[i] - b[i]*v[i-1] 
        y[i] = ( f[i] - b[i]*y[i-1])/w 

    for j in range(n-2,1,-1):
        y[j] = y[j] - v[j]*y[j+1] 

    return y


if __name__ == '__main__':
    solveNLSW_BOUS()