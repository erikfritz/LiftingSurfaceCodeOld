# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 11:47:06 2020

@author: fritzek
"""
#%% packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os
os.chdir("C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python/LiftingSurfaceCode")
from vortexRing import vortexRingVec
import geomParser
os.chdir("C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python/LiftingSurfaceCode")
import time

#%% input
# geometry
dspan = 10.
d1wingStart = np.array([0.,0.,-dspan/2])
d1wingEnd = np.array([0.,0.,dspan/2])
iElemI = 10
iElemJ = 20

# simulation parameters
dTimeBegin = 0.
dTimeEnd = 5.
dTimeStep = 0.05
d1Time = np.arange(dTimeBegin,dTimeEnd,dTimeStep)
iTimeSteps = len(d1Time)

dQ = 20.
dalpha = 5.
d2Q = np.ones([iTimeSteps,1])*np.array([np.cos(dalpha*np.pi/180)*dQ,np.sin(dalpha*np.pi/180)*dQ,0.])
drho = 1.225
dmu = 1.802e-5

# plotting
bPlotDev = False

#%% geometrical calculations
# wing
d1zLE = d1zTE = np.linspace(d1wingStart[2],d1wingEnd[2],100)
#d1elChord = np.ones_like(d1zLE) #straight chord
d1elChord = np.sqrt(1-np.square(d1zLE)/(np.square(dspan/2))) #elliptical chord
d1xLE = -0.25*d1elChord
d1yLE = np.zeros_like(d1xLE)
d1xTE = 0.75*d1elChord
d1yTE = np.zeros_like(d1xTE)

dRe = drho*dQ*(d1xTE-d1xLE)/dmu

# airfoil
d1data = np.fromfile('C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python/Airfoils/NACA0012.dat',dtype=float,sep=' ')
d1airfoilX = d1data[::2]
d1airfoilY = d1data[1::2]
iLEidx = np.argmin(d1airfoilX)
dLEval = d1airfoilX[iLEidx]
d1dataSSx = np.flip(d1airfoilX[0:iLEidx+1],0)
d1dataSSy = np.flip(d1airfoilY[0:iLEidx+1],0)
d1dataPSx = d1airfoilX[iLEidx:-1]
d1dataPSy = d1airfoilY[iLEidx:-1]

d2airfoilSS = np.vstack((np.linspace(0,1,101),np.interp(np.linspace(0,1,101),d1dataSSx,d1dataSSy)))
d2airfoilPS = np.vstack((np.linspace(0,1,101),np.interp(np.linspace(0,1,101),d1dataPSx,d1dataPSy)))
d2cambLine = np.vstack((np.linspace(0,1,101),d2airfoilPS[1,:]+0.5*(d2airfoilSS[1,:]-d2airfoilPS[1,:])))
d1cambLineGrid = np.interp(np.linspace(0,1,iElemI+1),d2cambLine[0,:],d2cambLine[1,:])

# panelling
d1spanPos = np.cos(np.linspace(np.pi,0,iElemJ+1))*dspan/2
d2gridX = np.zeros([iElemJ+1,iElemI+1])
d2gridY = np.zeros([iElemJ+1,iElemI+1])
d2gridZ = np.zeros([iElemJ+1,iElemI+1])
d2locLE = np.zeros([iElemJ+1,3])
d2locTE = np.zeros([iElemJ+1,3])
d1locChord = np.zeros(iElemJ+1)

for j in range(0,iElemJ+1):
    d2locLE[j,:] = np.array((np.interp(d1spanPos[j],d1zLE,d1xLE), np.interp(d1spanPos[j],d1zLE,d1yLE), d1spanPos[j]))
    d2locTE[j,:] = np.array((np.interp(d1spanPos[j],d1zTE,d1xTE), np.interp(d1spanPos[j],d1zTE,d1yTE), d1spanPos[j]))
    d1locChord[j] = np.sqrt(np.sum(np.square(d2locLE[j,:]-d2locTE[j,:])))
    d2locChordLine = np.linspace(d2locLE[j,:],d2locTE[j,:],iElemI+1)
    
#    d2gridX[j,:] = d2locChordLine[:,0]+d1cambLineGrid*np.sin(dalpha*np.pi/180)*d1locChord[j]
#    d2gridY[j,:] = d2locChordLine[:,1]+d1cambLineGrid*np.cos(dalpha*np.pi/180)*d1locChord[j]
    d2gridX[j,:] = d2locChordLine[:,0]+d1cambLineGrid*d1locChord[j]
    d2gridY[j,:] = d2locChordLine[:,1]+d1cambLineGrid*d1locChord[j]
    d2gridZ[j,:] = np.ones(iElemI+1)*d1spanPos[j]
    
d2gridY = np.nan_to_num(d2gridY)

#%%create meshgrid for VRs and CPs
#vortex rings
d2VRgridX,d2VRgridY,d2VRgridZ = geomParser.vrGrid(d2gridX,d2gridY,d2gridZ)
d2nX,d2nY,d2nZ = geomParser.nVec(d2VRgridX,d2VRgridY,d2VRgridZ)
d1nX = np.reshape(d2nX,(1,np.size(d2nX)),order='F').T
d1nY = np.reshape(d2nY,(1,np.size(d2nY)),order='F').T
d1nZ = np.reshape(d2nZ,(1,np.size(d2nZ)),order='F').T
#collocation points
d2CPgridX,d2CPgridY,d2CPgridZ=geomParser.cpGrid(d2gridX,d2gridY,d2gridZ)

#%%simulation
tic = time.time()
toc = np.zeros([iTimeSteps])
print('Progress   Simulated time   Wall time   No. of vortex elements')

st1bound = {}
st1wake = {}

for t in range(0,iTimeSteps):
    st1bound[d1Time[t]] = {}
    st1wake[d1Time[t]] = {}

    #influence coefficients bound vortices on themselves
    d2VelIndUB2B,d2VelIndVB2B,d2VelIndWB2B,_,_,_ = vortexRingVec(d2CPgridX,d2CPgridY,d2CPgridZ,
                                                                 d2VRgridX,d2VRgridY,d2VRgridZ,
                                                                 np.ones(np.size(d2CPgridX)))
    d2aIJ = np.tile(d1nX,(1,iElemI*iElemJ))*d2VelIndUB2B + np.tile(d1nY,(1,iElemI*iElemJ))*d2VelIndVB2B + np.tile(d1nZ,(1,iElemI*iElemJ))*d2VelIndWB2B
            
    d1TEX = d2VRgridX[:,-1]
    d1TEY = d2VRgridY[:,-1]
    d1TEZ = d2VRgridZ[:,-1]

    if t==0:
        st1wake[d1Time[t]]['X'] = d2VRgridX[:,-1]
        st1wake[d1Time[t]]['Y'] = d2VRgridY[:,-1]
        st1wake[d1Time[t]]['Z'] = d2VRgridZ[:,-1]
        st1wake[d1Time[t]]['Gamma'] = []
        
        d1RHS = -(d2Q[t,0]*d1nX + d2Q[t,1]*d1nY + d2Q[t,2]*d1nZ)
        st1bound[d1Time[t]]['Gamma'] = np.linalg.solve(d2aIJ, d1RHS)

    else: 
        #W2W
        if np.size(st1wake[d1Time[t-1]]['X']) > iElemJ+1: #wake consists of multiple rows of vortex rings
            d2velIndUW2W,d2velIndVW2W,d2velIndWW2W,_,_,_= vortexRingVec(st1wake[d1Time[t-1]]['X'],st1wake[d1Time[t-1]]['Y'],st1wake[d1Time[t-1]]['Z'],
                                                                        st1wake[d1Time[t-1]]['X'],st1wake[d1Time[t-1]]['Y'],st1wake[d1Time[t-1]]['Z'],
                                                                        st1wake[d1Time[t-1]]['Gamma'])
            d1velIndUW2W = np.sum(d2velIndUW2W,axis=1)
            d1velIndVW2W = np.sum(d2velIndVW2W,axis=1)
            d1velIndWW2W = np.sum(d2velIndWW2W,axis=1)
            
            st1wake[d1Time[t]]['Gamma'] = np.vstack([st1bound[d1Time[t-1]]['Gamma'][-iElemJ:],st1wake[d1Time[t-1]]['Gamma']])
            
        else: #wake consists of only one row of vortex rings
            d1velIndUW2W = np.zeros_like(d1TEX)
            d1velIndVW2W = np.zeros_like(d1TEY)
            d1velIndWW2W = np.zeros_like(d1TEZ)
            
            st1wake[d1Time[t]]['Gamma'] = st1bound[d1Time[t-1]]['Gamma'][-iElemJ:]

        #B2W
        d2velIndUB2W,d2velIndVB2W,d2velIndWB2W,_,_,_= vortexRingVec(st1wake[d1Time[t-1]]['X'],st1wake[d1Time[t-1]]['Y'],st1wake[d1Time[t-1]]['Z'],
                                                                    d2VRgridX,d2VRgridY,d2VRgridZ,
                                                                    st1bound[d1Time[t-1]]['Gamma'])
        d1velIndUB2W = np.sum(d2velIndUB2W,axis=1)
        d1velIndVB2W = np.sum(d2velIndVB2W,axis=1)
        d1velIndWB2W = np.sum(d2velIndWB2W,axis=1)
        
        #Convect wake
        d2wakeConvX = np.reshape((d2Q[t,0]*np.ones_like(d1velIndUB2W)+d1velIndUW2W+d1velIndUB2W)*dTimeStep,[iElemJ+1,t],order='F')
        d2wakeConvY = np.reshape((d2Q[t,1]*np.ones_like(d1velIndVB2W)+d1velIndVW2W+d1velIndVB2W)*dTimeStep,[iElemJ+1,t],order='F')
        d2wakeConvZ = np.reshape((d2Q[t,2]*np.ones_like(d1velIndWB2W)+d1velIndWW2W+d1velIndWB2W)*dTimeStep,[iElemJ+1,t],order='F')
       
        st1wake[d1Time[t]]['X'] = np.append([d1TEX],st1wake[d1Time[t-1]]['X'].T+d2wakeConvX.T,axis=0).T
        st1wake[d1Time[t]]['Y'] = np.append([d1TEY],st1wake[d1Time[t-1]]['Y'].T+d2wakeConvY.T,axis=0).T
        st1wake[d1Time[t]]['Z'] = np.append([d1TEZ],st1wake[d1Time[t-1]]['Z'].T+d2wakeConvZ.T,axis=0).T
        
        #W2B
        d2velIndUW2B,d2velIndVW2B,d2velIndWW2B,_,_,_= vortexRingVec(d2CPgridX,d2CPgridY,d2CPgridZ,
                                                            st1wake[d1Time[t]]['X'],st1wake[d1Time[t]]['Y'],st1wake[d1Time[t]]['Z'],
                                                            st1wake[d1Time[t]]['Gamma'])
        d1velIndUW2B = np.sum(d2velIndUW2B,axis=1)
        d1velIndVW2B = np.sum(d2velIndVW2B,axis=1)
        d1velIndWW2B = np.sum(d2velIndWW2B,axis=1)
        
        #circulation
        d1RHS = -(d2Q[t,0]*d1nX + d2Q[t,1]*d1nY + d2Q[t,2]*d1nZ + 
                 np.reshape(d1velIndUW2B,(np.size(d1velIndUW2B),1))*d1nX + 
                 np.reshape(d1velIndVW2B,(np.size(d1velIndVW2B),1))*d1nY + 
                 np.reshape(d1velIndWW2B,(np.size(d1velIndWW2B),1))*d1nZ)
        st1bound[d1Time[t]]['Gamma'] = np.linalg.solve(d2aIJ, d1RHS)
        
        d2Gamma = st1bound[d1Time[t]]['Gamma'].reshape([iElemJ,iElemI],order='F')
        d2GammaCor = np.hstack([d2Gamma[:,0].reshape([iElemJ,1]), d2Gamma[:,1:] - d2Gamma[:,0:-1]])
        st1bound[d1Time[t]]['GammaCor'] = d2GammaCor.reshape([1,np.size(d2Gamma)],order='F')
        
        #lift
        d2dX = 0.75*(d2VRgridX[1:,0:-1]-d2VRgridX[0:-1,0:-1])+0.25*(d2VRgridX[1:,1:]-d2VRgridX[0:-1,1:])
        d2dY = 0.75*(d2VRgridY[1:,0:-1]-d2VRgridY[0:-1,0:-1])+0.25*(d2VRgridY[1:,1:]-d2VRgridY[0:-1,1:])
        d2dZ = 0.75*(d2VRgridZ[1:,0:-1]-d2VRgridZ[0:-1,0:-1])+0.25*(d2VRgridZ[1:,1:]-d2VRgridZ[0:-1,1:])
        d2eNorm = np.sqrt(np.square(d2dX)+np.square(d2dY)+np.square(d2dZ))
        d2e1 = d2dX/d2eNorm
        d2e2 = d2dY/d2eNorm
        d2e3 = d2dZ/d2eNorm
        
        d2VtotU = d2Q[t,0] + np.sum(d2VelIndUB2B,axis=1).reshape([iElemJ,iElemI],order='F') + d1velIndUW2B.reshape([iElemJ,iElemI],order='F')
        d2VtotV = d2Q[t,1] + np.sum(d2VelIndVB2B,axis=1).reshape([iElemJ,iElemI],order='F') + d1velIndVW2B.reshape([iElemJ,iElemI],order='F')
        d2VtotW = d2Q[t,2] + np.sum(d2VelIndWB2B,axis=1).reshape([iElemJ,iElemI],order='F') + d1velIndWW2B.reshape([iElemJ,iElemI],order='F')
              
#        d2Fx = drho*(d2VtotW*d2e2 - d2VtotV*d2e3)*d2GammaCor
#        d2Fy = drho*(d2VtotU*d2e3 - d2VtotW*d2e1)*d2GammaCor
#        d2Fz = drho*(d2VtotV*d2e1 - d2VtotU*d2e2)*d2GammaCor
        d2Fx = drho*(d2VtotV*d2e3 - d2VtotW*d2e2)*d2GammaCor
        d2Fy = drho*(d2VtotW*d2e1 - d2VtotU*d2e3)*d2GammaCor
        d2Fz = drho*(d2VtotU*d2e2 - d2VtotV*d2e1)*d2GammaCor
       
        st1bound[d1Time[t]]['L'] = -d2Fy/np.cos(np.arctan(d2VtotV/d2VtotU))
        st1bound[d1Time[t]]['Ltot'] = np.sum(np.sum(st1bound[d1Time[t]]['L'],axis=1)*np.diff(d1spanPos))
        st1bound[d1Time[t]]['Fax'] = np.sum(np.sum(d2Fx,axis=1)*np.diff(d1spanPos))
        #lift coefficient
        d1dc = np.interp(d1spanPos[0:-1]+0.5*np.diff(d1spanPos),d1zLE,d1elChord)
        
#        st1bound[d1Time[t]]['Cl'] = np.sum(st1bound[d1Time[t]]['L'],axis=1)/(0.5*drho*np.sum(np.square(d2Q[t]))*d1dc)
        st1bound[d1Time[t]]['Cl'] = np.sum(st1bound[d1Time[t]]['L']/(0.5*drho*(np.square(d2VtotU)+np.square(d2VtotV))),axis=1)/d1dc

    toc[t] = time.time() - tic
    print(str("{:.2f}".format(round(100*(t+1)/iTimeSteps,2)))+'%\t\t'+
          str("{:.2f}".format(round(d1Time[t],3)))+'s\t\t'+
          str("{:.2f}".format(round(toc[t],3)))+'s\t\t'+
          str(iElemI*iElemJ+(t-1)*iElemJ))
    
#%% plotting
#geometry
import matlab.engine
eng = matlab.engine.start_matlab()

d1CPgridX = np.reshape(d2CPgridX,(1,np.size(d2CPgridX)),order='F')
d1CPgridY = np.reshape(d2CPgridY,(1,np.size(d2CPgridY)),order='F')
d1CPgridZ = np.reshape(d2CPgridZ,(1,np.size(d2CPgridZ)),order='F')

if bPlotDev:
    for tPlot in range(1,iTimeSteps):
        R = eng.plotGeometry(matlab.double(d1xLE.tolist()), matlab.double(d1yLE.tolist()), matlab.double(d1zLE.tolist()),
                         matlab.double(d1xTE.tolist()), matlab.double(d1yTE.tolist()), matlab.double(d1zTE.tolist()),
                         matlab.double(d1CPgridX.tolist()), matlab.double(d1CPgridY.tolist()), matlab.double(d1CPgridZ.tolist()),
                         matlab.double(d2gridX.tolist()), matlab.double(d2gridY.tolist()), matlab.double(d2gridZ.tolist()),
                         matlab.double(d2VRgridX.tolist()), matlab.double(d2VRgridY.tolist()), matlab.double(d2VRgridZ.tolist()),
                         matlab.double(st1wake[d1Time[tPlot]]['X'].tolist()),
                         matlab.double(st1wake[d1Time[tPlot]]['Y'].tolist()),
                         matlab.double(st1wake[d1Time[tPlot]]['Z'].tolist()),
                         matlab.double(st1wake[d1Time[tPlot]]['Gamma'].tolist()),
                         matlab.double([d2Q[-1,0]]),matlab.double([dTimeStep]),iTimeSteps)
else:
    tPlot = t
    R = eng.plotGeometry(matlab.double(d1xLE.tolist()), matlab.double(d1yLE.tolist()), matlab.double(d1zLE.tolist()),
                         matlab.double(d1xTE.tolist()), matlab.double(d1yTE.tolist()), matlab.double(d1zTE.tolist()),
                         matlab.double(d1CPgridX.tolist()), matlab.double(d1CPgridY.tolist()), matlab.double(d1CPgridZ.tolist()),
                         matlab.double(d2gridX.tolist()), matlab.double(d2gridY.tolist()), matlab.double(d2gridZ.tolist()),
                         matlab.double(d2VRgridX.tolist()), matlab.double(d2VRgridY.tolist()), matlab.double(d2VRgridZ.tolist()),
                         matlab.double(st1wake[d1Time[tPlot]]['X'].tolist()),
                         matlab.double(st1wake[d1Time[tPlot]]['Y'].tolist()),
                         matlab.double(st1wake[d1Time[tPlot]]['Z'].tolist()),
                         matlab.double(st1wake[d1Time[tPlot]]['Gamma'].tolist()),
                         matlab.double([d2Q[-1,0]]),matlab.double([dTimeStep]),iTimeSteps)

#circulation and lift
fig2, (ax2,ax3) = plt.subplots(1,2)
ax2.clear()
ax3.clear()
for i in range(0,iElemI):
    ax2.plot(d1spanPos[0:-1]+0.5*np.diff(d1spanPos), st1bound[d1Time[t]]['GammaCor'].reshape([iElemJ,iElemI],order='F')[:,i], label='i = %s' % i)
    ax3.plot(d1spanPos[0:-1]+0.5*np.diff(d1spanPos), st1bound[d1Time[t]]['L'][:,i])
    
ax2.plot(d1spanPos[0:-1]+0.5*np.diff(d1spanPos), np.sum(st1bound[d1Time[t]]['GammaCor'].reshape([iElemJ,iElemI],order='F'),axis=1),'k', label='total')
ax3.plot(d1spanPos[0:-1]+0.5*np.diff(d1spanPos), np.sum(st1bound[d1Time[t]]['L'],axis=1),'k')
ax2.set_xlabel('Span [m]')
ax2.set_ylabel('Circulation $\Gamma$ [$\mathrm{m}^2/\mathrm{s}$]')
ax3.set_xlabel('Span [m]')
ax3.set_ylabel('Lift $L$ [N/m]')
ax2.legend()

#lift coefficient
fig3 = plt.figure('Lift coefficient')
ax4 = fig3.add_subplot()
ax4.clear()
ax4.plot(d1spanPos[0:-1]+0.5*np.diff(d1spanPos), st1bound[d1Time[t]]['Cl'],'k')
ax4.set_xlabel('Span [m]')
ax4.set_ylabel('$c_l$ [-]')

#forces
fig4 = plt.figure('Forces')
ax5 = fig4.add_subplot()
ax5.clear()
d1Ltot = np.empty(0)
d1Fax = np.empty(0)
for i in d1Time[1:]:
    d1Ltot = np.append(d1Ltot,st1bound[i]['Ltot'])
    d1Fax = np.append(d1Fax,st1bound[i]['Fax'])
ax5.plot(d1Time[1:],d1Ltot, label='Total lift')
ax5.plot(d1Time[1:],d1Fax, label='Axial force')
ax5.set_xlabel('Time [s]')
ax5.set_ylabel('Forces [N]')
ax5.legend()

d1Rz = d1spanPos[0:-1]+0.5*np.diff(d1spanPos)
d1Rg = np.sum(st1bound[d1Time[t]]['GammaCor'].reshape([iElemJ,iElemI],order='F'),axis=1)
d1Rc = st1bound[d1Time[t]]['Cl']
d2out = np.vstack([d1Rz,d1Rg,d1Rc]).T
d2out1 = np.vstack([d1Rz,d1Rg,d1Rc])
np.savetxt('test.txt',d2out)
for i in range(0,len(d1Rz)):
    print('('+str(d1Rz[i])+','+str(d1Rg[i])+')')
    
for i in range(0,len(d1Rz)):
    print('('+str(d1Rz[i])+','+str(d1Rc[i])+')')