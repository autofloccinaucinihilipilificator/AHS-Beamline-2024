import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

N=int(input('N?'))
resolution=float(input('resolution?'))
timestep=float(input('timestep?'))
kBoltz=1
q=1
epislon0=1
mu0=1
mass=1
v0=np.array([1,0,0])
p0=np.array([0,0,0]) #if you change this you have to also change the initial conditions function
rho=q/(resolution**3)

J = np.array([[[[0,0,0] for k in range(-N,N+1)] for j in range(-N,N+1)] for i in range(-N,N+1)])

#sets initial conditions based on charge density
def intialCond(i,j,k):
    if i==0 and j==0 and k==0:
        return([0,0,0])
    else:
        return([kBoltz*q*i/((i**2+j**2+k**2)**(3/2)),kBoltz*q*j/((i**2+j**2+k**2)**(3/2)),kBoltz*q*k/((i**2+j**2+k**2)**(3/2))])

#average value of the eletric field at a point
def avgValEfield(p):
    return(IndCheck(efield,2*N,N+(0.5+(p[0]/resolution))//1,N+(0.5+(p[1]/resolution))//1,N+(0.5+(p[2]/resolution))//1))
    #average over 8 points: return((IndCheck(efield,2*N,1+N+((p[1]/resolution))//1,1+N+((p[2]/resolution))//1,1+N+((p[3]/resolution))//1)[i]+IndCheck(efield,2*N,1+N+((p[1]/resolution))//1,N+((p[2]/resolution))//1,1+N+((p[3]/resolution))//1)[i]+IndCheck(efield,2*N,N+((p[1]/resolution))//1,1+N+((p[2]/resolution))//1,1+N+((p[3]/resolution))//1)[i]+IndCheck(efield,2*N,1+N+((p[1]/resolution))//1,1+N+((p[2]/resolution))//1,N+((p[3]/resolution))//1)[i]+IndCheck(efield,2*N,N+((p[1]/resolution))//1,N+((p[2]/resolution))//1,1+N+((p[3]/resolution))//1)[i]+IndCheck(efield,2*N,N+((p[1]/resolution))//1,1+N+((p[2]/resolution))//1,N+((p[3]/resolution))//1)[i]+IndCheck(efield,2*N,1+N+((p[1]/resolution))//1,N+((p[2]/resolution))//1,N+((p[3]/resolution))//1)[i]+IndCheck(efield,2*N,N+((p[1]/resolution))//1,N+((p[2]/resolution))//1,N+((p[3]/resolution))//1)[i])/8)

#average value of the magnetic field at a point
def avgValMfield(p):
    return(IndCheck(mfield,2*N-1,N+(p[0]/resolution)//1,N+(p[1]/resolution)//1,N+(p[2]/resolution)//1))

#completed index check
def IndCheck(L,n,i,j,k):
    if 0<=i<=n and 0<=j<=n and 0<=k<=n:
        return(np.array(L[int(i)][int(j)][int(k)]))
    else:
        return(np.array([0,0,0]))

#current density check
def Jcheck(i,j,k):
    if i == N+(0.5+(p0[0]/resolution))//1 and j == N+(0.5+(p0[1]/resolution))//1 and k == N+(0.5+(p0[2]/resolution))//1:
        return(rho * v0)
    else:
        return(np.array([0,0,0]))

#directional derivatives
def ddx(L,n,i,j,k,xyz):
    return(np.array((IndCheck(L,n,i,j,k)[xyz] + IndCheck(L,n,i,j+1,k)[xyz] + IndCheck(L,n,i,j,k+1)[xyz] + IndCheck(L,n,i,j+1,k+1)[xyz])/(4*resolution)))

def ddy(L,n,i,j,k,xyz):
    return(np.array((IndCheck(L,n,i,j,k)[xyz] + IndCheck(L,n,i+1,j,k)[xyz] + IndCheck(L,n,i,j,k+1)[xyz] + IndCheck(L,n,i+1,j,k+1)[xyz])/(4*resolution)))

def ddz(L,n,i,j,k,xyz):
    return(np.array((IndCheck(L,n,i,j,k)[xyz] + IndCheck(L,n,i+1,j,k)[xyz] + IndCheck(L,n,i,j+1,k)[xyz] + IndCheck(L,n,i+1,j+1,k)[xyz])/(4*resolution)))


#curl
def curl(F,s,e,n):
    return(np.array([[[[ddy(F,n,i,j,k,2)-ddz(F,n,i,j,k,1),ddz(F,n,i,j,k,0)-ddx(F,n,i,j,k,2),ddx(F,n,i,j,k,1)-ddy(F,n,i,j,k,0)] for k in range(s,e)] for j in range(s,e)] for i in range(s,e)]))

#make list of points for electric field
efield = np.array([[[intialCond(i*resolution,j*resolution,k*resolution) for k in range(-N,N+1)] for j in range(-N,N+1)] for i in range(-N,N+1)])

#make list of points for magnetic field
mfield = np.array([[[[0,0,0] for k in range(-N,N)] for j in range(-N,N)] for i in range(-N,N)])

plt.savefig('plot.png')

iterations=int(input('iterations?'))

while iterations>0:
    print(iterations)
    print(p0)
    print(v0)
    #time derivative of magnetic field
    mfield = mfield - timestep * curl(efield,0,2*N,2*N)
    #efield = [[[[(efield[i][j][k][xyz]-timestep*curl(efield,0,2*N,2*N)[i][j][k][xyz] for xyz in range(2))] for k in range(2*N)] for j in range(2*N)] for i in range(2*N)]
    #time derivative of electric field
    efield = efield + (timestep/epislon0) * ((curl(mfield,-1,2*N,2*N-1)/mu0) - J)
    #mfield = [[[[(mfield[i][j][k][xyz]+timestep*(1/epislon0)*(((curl(mfield,-1,2*N,2*N-1)[i][j][k][xyz])/mu0)-J[i][j][k][xyz])) for xyz in range(2)] for k in range(2*N)] for j in range(2*N)] for i in range(2*N)]
    #update particle velocity
    v0 = v0 + (timestep*q/mass) * (avgValEfield(p0)+np.cross(v0,avgValMfield(p0)))
    #update particle position
    p0 = p0 + timestep * v0
    #update current density
    J = [[[Jcheck(i,j,k) for k in range(-N,N+1)] for j in range(-N,N+1)] for i in range(-N,N+1)]
    iterations -= 1
print(iterations)
print(p0)
print(v0)

ax = plt.figure().add_subplot(projection='3d')

#Make the grid
x, y, z = np.meshgrid(np.arange(-N, N, 1),
                      np.arange(-N, N, 1),
                      np.arange(-N, N, 1))
#Make the direction data for the arrows
u = np.array([[[mfield[int(x[i][j][k]+N)][int(y[i][j][k]+N)][int(z[i][j][k]+N)][0] for k in range(-N,N)] for j in range(-N,N)] for i in range(-N,N)])
v = np.array([[[mfield[int(x[i][j][k]+N)][int(y[i][j][k]+N)][int(z[i][j][k]+N)][1] for k in range(-N,N)] for j in range(-N,N)] for i in range(-N,N)])
w = np.array([[[mfield[int(x[i][j][k]+N)][int(y[i][j][k]+N)][int(z[i][j][k]+N)][2] for k in range(-N,N)] for j in range(-N,N)] for i in range(-N,N)])

#plots the current magnetic field (super hard to see lol, can't even tell if its right)
ax.quiver(x, y, z, u, v, w, length=1/resolution, normalize=True)
plt.savefig('plot.png')
