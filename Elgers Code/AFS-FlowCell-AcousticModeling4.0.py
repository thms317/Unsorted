# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 09:12:20 2015

@author: Elger Vlieg
"""

'''

The program has four modes

1)FreeDraw
Draw an acoustic wave or standing wave of frequency (f) at time (t) in the 
flow cell

2)FrequencyPerfomance
Calculate the effective gradient at the sample over a linear sampled frequency 
range

3)LayerThicknessPerformance
Calculate the maximum and minimum gradient at the sample over a linear sampled 
frequency range and a linear sampled layer thickness range

4)SaveCoefficients
Creates two LaTeX compatible tables in text files with the transmission and
reflection coefficients of the flow cell

Note that all processes can be aborted by inserting crtl+c in the therminal!
'''


########################    Flow Cell parameters     ##########################
'''
Flowcellname = "Your FC name here"

parameters go like [layer1, layer2, layer3...]
                    [Air, ITO, Piezo, ITO, pdms...]

thickness = [1, 0.6, 1000, 0.6]         #um some arbitrary number for the outer ones
speed = [343.14, 4100, 6935, 4100]      #speed of sound m/s
density = [0.0012754, 7.14, 4.64, 7.14] #g/cm3

Liq = 4                 #number of liquid layer
LeftRight = "Right"     #String. The side of the liquid layer the sample is at (Left or Right)
Piezo = 2               #number of piezoelectric layer
'''
'''
Flowcellname = "New design"
thickness = [0,1, 1000,1,  60, 170, 0] 
speed = [343.14,4100, 6935,4100,  1482.28, 5.61*10**3, 343.14] 
density = [0.0012754,7.14, 4.64,7.14,  0.9982071, 2.51, 0.0012754] 
Liq = 5 
LeftRight = "Right" 
Piezo = 3 
'''

Flowcellname = "Prototype Flow Cell with oil immersion objective"
thickness = [0,1, 1000, 1,  500, 300, 170, 0] 
speed = [343.14,4100, 6935,4100,  1076.5, 1482.28, 5.61*10**3, 1315.7] 
density = [0.0012754,7.14, 4.64,7.14,  0.97, 0.9982071, 2.51, 1.093] 
Liq = 6 
LeftRight = "Right" 
Piezo = 3 

'''
Flowcellname = "Prototype Flow Cell with air objective"
thickness = [0,1, 1000,1,  500, 300, 170, 0] 
speed = [343.14,4100, 6935,4100,  1076.5, 1482.28, 5.61*10**3, 343.14] 
density = [0.0012754,7.14, 4.64,7.14,  0.97, 0.9982071, 2.51, 0.0012754] 
Liq = 6 
LeftRight = "Right" 
Piezo = 3 
'''


####################       Please choose one mode       #######################

#True is active
FreeDraw =                  0   #Boolean
FrequencyPerformance =      0   #Boolean    
LayerThicknessPerformance = 1   #Boolean
SaveCoefficients =          0   #Boolean


##################    Coefficient table specifications     ####################


MainWorkspace = "C:\Users\Elger\Dropbox\Bachelor Onderzoek\FlowCellCoefficients\\"
                        #Directory where the coefficients are saved
dec = 6                 #Amount of decimals in the coefficients
deletezeros = 0         #Boolean. All ending zeros are removed if true


#######################    Diagram configuration     ##########################


samplerate = 10000      #Drawn waves sample points(determines smoothness of wave lines)

LegendX = 0.65          #Legend X-position (0 against left side, 1 over the right side)
LegendY = 0.90          #Legend Y-position (0 under the bottom, 1 against top) 

cut = 1                 #Boolean, cuts of spikes for the generated pressure indication
                        #in the FrequencyPerformance mode

#####################    Caculation specifications     ########################


#------------------------------Free Draw--------------------------------------#

f = 6                   #Frequency to be drawn MHz #doesnt always work for f<<1
df = 1./(f*10**6)       #will show up as multiple lines when fitCheck is enabled
t = [0,df/6.]           #Array of floats. All the times to draw for (only for non standing waves). For the right
                        #outgoing wave we have P = cos( 2pi f t).
standing = 1            #Boolean. Draw standing waves or real acoustic pressure profile at times in t? 
                        #Max amplitude is drawn.
fitCheck = 0            #Draw standing profile for multiple (pre-defined) times to check for errors

#-----------------Frequency and thickness performance-------------------------#

bf=0                    #First frequency MHz
ef=20                   #Last frequency MHz
sf=2001                #Amount of linear samples to take. Please sample no lower than 1 kHz/sample
                        #+- 0.003sec/sample at default settings(Intel Core i5 3210M @ 2.50GHz)  

#--------------------------Thickness performance------------------------------#

bt = 0                  #First Thickness um
et = 1200               #Last Thickness  um
st = 121                #Amount of linear samples to take. Warning st*sf*0.003sec for calculation time at default settings!  
Layer = Piezo           #What layer to differ


#####################    Performance configuration     ########################


ErrorLimit = 10**-2     #Relative to pressure/velocity amplitude fitting error that is considered as big
     
NonFitTres = 10**-7     #Performance threshold for the instant phase parameter calculation
                        #(lower the tighter the selection). More than 10**-5 not recommended.
size = 1000             #The phase is fitted when the instant calculation is not accurate enough. 
                        #The residue is linear calculated over 'size' samples   
accuracy =  4           #The linear residue minimization is 'accuracy-1' times repeated over the interval 
                        #of the previous minimum with its neighbours
 
'''
Default
ErrorLimit = 10**-2    
NonFitTres = 10**-7    The higher the higher the calculation speed. Higher than 10**-5 not recommended
                       at default settings holds a performance increase of +- 40%
size = 1000            The lower the higher "" ""
accuracy = 4           The lower the higher "" ""
'''



###############################################################################
###############################################################################

  
import numpy as np
from matplotlib import pyplot as plt
import cmath as cm
import datetime as dt
  
def e(f,d,c,p):
    return np.e**(1j*2*np.pi*f*d/c)
    
def WF(x,A,B,e0,e1,f,c):
    return A*np.cos(2*np.pi * f * x/c + e0) + B*np.cos(-2*np.pi * f * x/c + e1)

def VEL   (x,A,B,e0,e1,f,c):
    return A*np.cos(2*np.pi * f * x/c + e0) - B*np.cos(-2*np.pi * f * x/c + e1)
        
def SWF(x,A,e,c,f):   
    return np.abs((A*np.cos(np.pi*2*f*x/c + e)))
    
def PV(x,A,B,e0,e1,f,c,p):
    return 1/(p*c) * (A*np.cos(2*np.pi * f * x/c + e0) - B*np.cos(-2*np.pi * f * x/c + e1))
    
 
def EffectiveGradient(x,A,B,e0,e1,f,c):
    
    k = np.pi * 2 * f / c
    e = (e0-e1)/2.
    
    return -2 *k*A*B*np.sin(k*x + e) * np.cos(k*x+e)
    
def BuildTransferMatrix(f,d,c,p):
        
    c11 = 1/2. * (e(f,d,c,p) + e(f,-d,c,p))  
    c12 = (p*c/2.) * (e(f,-d,c,p)-e(f,d,c,p))
    c21 = 1/(2.*p*c) * (e(f,-d,c,p)-e(f,d,c,p))
    c22 = 1/2. * (e(f,-d,c,p)+e(f,d,c,p))
          
    return np.array( ((c11,c12), (c21, c22)) )
    
def BuildInverseTransferMatrix(f,d,c,p):
         
    c11 = 1/2. * (e(f,d,c,p) + e(f,-d,c,p))  
    c12 = -(p*c/2.) * (e(f,-d,c,p)-e(f,d,c,p))
    c21 = -1/(2.*p*c) * (e(f,-d,c,p)-e(f,d,c,p))
    c22 = 1/2. * (e(f,-d,c,p)+e(f,d,c,p))
           
    return np.array( ((c11,c12), (c21, c22)) )
    
    
def GetMatrices(f,d,c,p):
    
    M = np.zeros((np.size(d)-2)*4).reshape((np.size(d)-2),2,2) * 1j
    
    for i in range(np.size(d)-2):
        M[i] = BuildTransferMatrix(f,d[i+1],c[i+1],p[i+1])
    
    return M

def GetInverseMatrices(f,d,c,p):
    
    Mi = np.zeros((np.size(d)-2)*4).reshape((np.size(d)-2),2,2) * 1j
    
    for i in range(np.size(d)-2):
        Mi[i] = BuildInverseTransferMatrix(f,d[i+1],c[i+1],p[i+1])
    
    return Mi
    
    
def CalculateMatchingConstants(Piezo,p,c,d,f):
    
    Mi = GetInverseMatrices(f,d,c,p)
    
    M0 = [1,0],[0,1]
    M0 = np.array(M0)
    for i in range(Piezo-2):

        M0 = np.dot(Mi[i],M0)

    return [M0[0,0] - M0[0,1]*(1/(p[0]*c[0])), M0[1,0] - M0[1,1]*(1/(p[0]*c[0]))]
    
    
def GeneratedTransferVector(Piezo,p,c,d,f):
    
    cc = CalculateMatchingConstants(Piezo,p,c,d,f)
    
    r = cc[0]/cc[1]
    
    k = 2*np.pi * f / c[Piezo - 1]
    d1 = d[Piezo - 1]
    p1 = p[Piezo - 1]
    c1 = c[Piezo - 1]
    
    df =1/(1+r/(p1*c1)+np.e**(-1j*k*d1)*(r/(p1*c1)-1))
    Pt = r/(1j*p1*c1) * np.sin(k*d1) - np.cos(k*d1)
    Pb = 1j*p1*c1*np.sin(k*d1) + r*np.cos(k*d1)
    
    return df*np.array([Pt,Pb])
    
 
def GetGeneratedWave(Piezo,p,c,d,f,PU):  
    
    V = GeneratedTransferVector(Piezo,p,c,d,f) 
    
    return np.dot(V,PU[Piezo-1])
    
    
def GetBorderSpecs(f,d,c,p,t,Piezo):
    PU = np.zeros((np.size(d)-1)*2).reshape((np.size(d)-1),2)   * 1j
    
    M = GetMatrices(f,d,c,p)
    
    PiezoVector = GeneratedTransferVector(Piezo,p,c,d,f)
    
    tempi = np.size(d) - 1
    PU[np.size(d)-2] = [np.e**(1j*2*np.pi*f*t), (np.e**(1j*2*np.pi*f*t))/(p[tempi]*c[tempi])]
    
    
    for i in range(np.size(d)-2):
        index = np.size(d) - 3 - i

        
        M00 = M[index,0,0]
        M01 = M[index,0,1]
        M10 = M[index,1,0]
        M11 = M[index,1,1]
        
        PU0 = PU[index+1,0]
        PU1 = PU[index+1,1]
        
        
        PU[index] = [PU0*M00 + PU1*M01,PU0*M10+PU1*M11]
        
        if(index == Piezo-2):
            PiezoFactor = PiezoVector[0] * PU0 + PiezoVector[1] * PU1 
            k = 2*np.pi*f/c[index+1]
            
            PU[index][0] = PU[index][0] + PiezoFactor* (1-np.e**(-1j*k*d[index+1]))
            PU[index][1] = PU[index][1] + PiezoFactor*(-1-np.e**(-1j*k*d[index+1]))/(p[index+1]*c[index+1])
        

            
    return PU


def MakeRealXDistance(d):
    
    temp = np.zeros(np.size(d)-1)
    
    for i in range(np.size(temp)):
        temp[i] = np.sum(d[0:i+1])
    
    temp = temp - d[0]
    
    return temp


def GetMinIndex(l):
    temp= l[0]
    index = 0
    for i in range(np.size(l)):
        if(l[i] < temp):
            temp = l[i]
            index = i
    
    
    return index,temp

def GetMaxIndex(l):
    temp= l[0]
    index = 0
    for i in range(np.size(l)):
        if(l[i] > temp):
            temp = l[i]
            index = i
    
    
    return index,temp
    
def FitPhaseShift0(f,c,p,P1,P2,U1,U2,a,b):
    
    const = (P1+U1*c*p)/(P2+U2*c*p)

    
    x = np.linspace(0,np.pi*2,size+1)
    y = np.abs(const - (np.cos(2*np.pi*f*a/(c*1.) + x)/np.cos(2*np.pi*f*b/(c*1.) + x)))
    index,mini = GetMinIndex(y)
    
    for i in range(accuracy-1):
        j = i+1
        dp = 2*np.pi/(size**j)
        x2 = np.linspace(x[index]-dp,x[index]+dp,size+1)
        y2 = np.abs(const - (np.cos(2*np.pi*f*a/(c*1.) + x2)/np.cos(2*np.pi*f*b/(c*1.) + x2)))
        index,mini = GetMinIndex(y2)
        x = x2
    
    phase = x[index] 
    mini = mini
    
    return phase,mini/const

def FitPhaseShift1(f,c,p,P1,P2,U1,U2,a,b):
    
    const = (P1-U1*c*p)/(P2-U2*c*p)
    

    x = np.linspace(0,np.pi*2,size+1)
    y = np.abs(const - (np.cos(-2*np.pi*f*a/(c*1.) + x)/np.cos(-2*np.pi*f*b/(c*1.) + x)))
    index,mini = GetMinIndex(y)

    
    for i in range(accuracy-1):
        j = i+1
        dp = 2*np.pi/(size**j)
        x2 = np.linspace(x[index]-dp,x[index]+dp,size+1)
        y2 = np.abs(const - (np.cos(-2*np.pi*f*a/(c*1.) + x2)/np.cos(-2*np.pi*f*b/(c*1.) + x2)))
        index,mini = GetMinIndex(y2)
        x = x2
        
    phase = x[index] 
    mini = mini
    
    return phase,mini



def ArccosCof(f,c,p,a,b,C):
    k = np.pi * 2 * f/c   
    return -(np.sin(k*a) - C*np.sin(k*b))/(np.sqrt(-2*C*(np.sin(k*a)*np.sin(k*b) + np.cos(k*a)*np.cos(k*b)) + 1 + C**2))

def ArccosCof2(f,c,p,a,b,C):
    k = np.pi * 2 * f/c   
    return (np.sin(k*a) - C*np.sin(k*b))/(np.sqrt(-2*C*(np.sin(k*a)*np.sin(k*b) + np.cos(k*a)*np.cos(k*b)) + 1 + C**2))

def GetPhaseShift(f,c,p,a,b,C):
    k = np.pi * 2 * f/c    
    return np.arccos(-(np.sin(k*a) - C*np.sin(k*b))/(np.sqrt(-2*C*(np.sin(k*a)*np.sin(k*b) + np.cos(k*a)*np.cos(k*b)) + 1 + C**2)))
    
    
def GetWaves(PU,d,c,f,p,Piezo,pg,TO,Liq):
    
    PU = np.real(PU)
    AmpPhase = np.zeros((np.size(PU[:,1]))*2).reshape((np.size(PU[:,1])),2)
    ABee = np.zeros((np.size(PU[:,1])-1)*4).reshape((np.size(PU[:,1])-1),4)
    
    x = MakeRealXDistance(d)
    error = False
    
    a0 = 1000
    b0 = 1000
    a1 = 1000
    b1 = 1000
    e0 = 1000
    e1 = 1000
    
    for i in range(np.size(AmpPhase[:,1])-1):
        
        if(TO):
            i = Liq - 2
            
        P1 = PU[i,0]
        U1 = PU[i,1]
        P2 = PU[i+1,0]
        U2 = PU[i+1,1]
        
        if(i == Piezo-2):
            P1 = P1 - np.real(pg)
            U1 = U1 + np.real(pg)/(p[Piezo-1]*c[Piezo-1])
            P2 = P2 - np.real(pg)
            U2 = U2 - np.real(pg)/(p[Piezo-1]*c[Piezo-1])
        
        a=x[i]
        b=x[i+1]


        c0 = (P1+c[i+1]*p[i+1]*U1)/(P2+c[i+1]*p[i+1]*U2)
        cof = ArccosCof2(f,c[i+1],p[i+1],a,b,c0)
        
        
        if (np.abs(cof) < 1):
            e0 = np.arccos(cof)
            b0 = np.cos(2*np.pi*f*b/c[i+1] + e0)
            a0 = np.cos(2*np.pi*a*f/c[i+1] + e0)                        
            
            
        if (np.abs(cof) > 1 or (np.abs((a0/b0 - c0)/c0)) > NonFitTres ):           
            e0,error0 = FitPhaseShift0(f,c[i+1],p[i+1],P1,P2,U1,U2,a,b)
            b0 = np.cos(2*np.pi*f*b/c[i+1] + e0)
            a0 = np.cos(2*np.pi*a*f/c[i+1] + e0)            
                       

        c1 = (P1-c[i+1]*p[i+1]*U1)/(P2-c[i+1]*p[i+1]*U2)
        cof1 = ArccosCof(f,c[i+1],p[i+1],-a,-b,c1)            
        
        if (np.abs(cof1) < 1):      
            e1 = np.arccos(cof1)
            b1 = np.cos(-2*np.pi*f*b/c[i+1] + e1) 
            a1 = np.cos(-2*np.pi*f*a/c[i+1] + e1) 
            
        if (np.abs(cof1) > 1 or (np.abs((a1/b1 - c1)/c1)) > NonFitTres ):
            e1,error1 = FitPhaseShift1(f,c[i+1],p[i+1],P1,P2,U1,U2,a,b)       
            b1 = np.cos(-2*np.pi*f*b/c[i+1] + e1) 
            a1 = np.cos(-2*np.pi*f*a/c[i+1] + e1) 
            
        
        A = (P1 + U1*p[i+1]*c[i+1])/ (2*a0)
        
        A1 = (P2+U2*p[i+1]*c[i+1])/(2*b0)
        
        B = (P1 - U1*p[i+1]*c[i+1])/ (2*a1)
        
        B1 = (P2 - U2*p[i+1]*c[i+1])/(2*b1)
        
        
        A = (A+A1)/2.
        B = (B+B1)/2.
        
        tempP1 = np.abs((WF(a,A,B,e0,e1,f,c[i+1]) - P1)/P1)
        tempP2 = np.abs((WF(b,A,B,e0,e1,f,c[i+1]) - P2)/P2)
        if( tempP1 > ErrorLimit or  tempP2 > ErrorLimit):
            error = True
            
            if(not TO):
                print ''
                print 'LARGE FITTING ERROR'
                print '##Layer ',i+1,'##'
                print ''
                print 'Relative acoustic pressure error'
                print 'DP1: ',(WF(a,A,B,e0,e1,f,c[i+1]) - P1)/P1 
                print 'DP2: ',(WF(b,A,B,e0,e1,f,c[i+1]) - P2)/P2
        
        tempU1 = np.abs((PV(a,A,B,e0,e1,f,c[i+1],p[i+1]) - U1)/U1)
        tempU2 = np.abs((PV(b,A,B,e0,e1,f,c[i+1],p[i+1]) - U2)/U2)
        if(( tempU1 > ErrorLimit or  tempU2 > ErrorLimit)):
            error = True
            
            
            if(not TO):
                print ''
                print 'LARGE FITTING ERROR'
                print '##Layer ',i+1,'##'
                print ''  
                print 'Relative normal velocity error'
                print 'DU1: ',(PV(a,A,B,e0,e1,f,c[i+1],p[i+1]) - U1)/U1
                print 'DU2: ',(PV(b,A,B,e0,e1,f,c[i+1],p[i+1]) - U2)/U2
            
        
        
        if(A < 0):
            e0 = e0+np.pi
            A = A*-1
    
        if(B < 0):
            e1 = e1+np.pi
            B = B*-1    
        
        ABee[i] = [A,B,e0,e1]
        

        if(i < Piezo-2):
            AmpPhase[i] = [A+B, (-e1+e0)/2.]

            
        if(i > Piezo-2):
            AmpPhase[i+1] = [A+B, (-e1+e0)/2.]
       
       
        if(i == Piezo-2):
             
            b = B*np.e**(1j*(e1))         
            
            
            cosAe = cm.polar(pg+b)
            
            C = cosAe[0]
            ee1 = cosAe[1] 
            
            AmpPhase[i] = [A+C, (-ee1+e0)/2.]
   
   
        if(i == Piezo-2):
            
            dp = d[Piezo-1]
            k = np.pi * 2 * f/c[Piezo-1]
            
            a = A*np.e**(1j*(e0))         
            g = pg*np.e**(-1j*k*dp)                
            
            cosAe = cm.polar(g+a)
            
            D = cosAe[0]
            ee0 = cosAe[1]           
            AmpPhase[i+1] = [B+D, (-e1+ee0)/2.]
            
        if(TO):
            return AmpPhase,ABee,error
        
    
    return AmpPhase,ABee,error


def MakeImageStanding(d,f,AmpPhase,c):   
    x = MakeRealXDistance(d)
    
    
    maxy = 0

    for i in range(np.size(d) - 2):
        
        if (i < Piezo -2):
            tempx = np.linspace(x[i],x[i+1],samplerate)
    
            y = SWF(tempx, AmpPhase[i,0],AmpPhase[i,1],c[i+1],f)
        
            
                
            plt.plot(tempx,y)
            if(np.max(y) > maxy):
                maxy = np.max(y)
                
        if (i > Piezo-2):
            tempx = np.linspace(x[i],x[i+1],samplerate)
            
            y = SWF(tempx, AmpPhase[i+1,0],AmpPhase[i+1,1],c[i+1],f)
        
            
                
            plt.plot(tempx*10**6,y)
            
            if(np.max(y) > maxy):
                maxy = np.max(y)     

        if (i == Piezo-2):
            
            tempx = np.linspace(x[i],x[i]+(x[i+1]-x[i])/2.,samplerate)
    
            y = SWF(tempx, AmpPhase[i,0],AmpPhase[i,1],c[i+1],f)
            
            plt.plot(tempx*10**6,y)
            if(np.max(y) > maxy):
                maxy = np.max(y)     

            tempx = np.linspace(x[i]+(x[i+1]-x[i])/2.,x[i+1],samplerate)
    
            y = SWF(tempx, AmpPhase[i+1,0],AmpPhase[i+1,1],c[i+1],f)
             
            plt.plot(tempx*10**6,y, color = 'green')
            
            if(np.max(y) > maxy):
                maxy = np.max(y)   
                
                   
    for i in range(np.size(x)):
        plt.plot([x[i]*10**6,x[i]*10**6],[0,maxy])
    
    tempx = [x[Piezo-2]*10**6 + (x[Piezo-1]-x[Piezo-2])*10**6/2.,x[Piezo-2]*10**6 + (x[Piezo-1]-x[Piezo-2])*10**6/2.]
    plt.plot( tempx ,[0,maxy], color = 'black', ls = '--')

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlim(-np.max(x)*10**6/7. , np.max(x)*10**6 + np.max(x)*10**6/7.)
    plt.ylabel('Normalized acoustic pressure')
    plt.xlabel('Distance $(\mu m)$')
    plt.title('TM modelled standing pressure profile: f = ' + str(f/10.**6)+'MHz')



def gradient(ABee,Liq,LeftRight,d,f):
    
        
    x = MakeRealXDistance(d)        
    
    tempLiq = Liq
    
    if(Piezo > Liq):
        tempLiq = Liq-1
    
    Ae = ABee[tempLiq-2,0]
    Be = ABee[tempLiq-2,1]
    e0e = ABee[tempLiq-2,2]
    e1e = ABee[tempLiq-2,3]
    

    
    if(LeftRight == 'Left'):
        X = x[Liq-2]

        
    if(LeftRight == 'Right'):
        X = x[Liq-1]
    
    
    return EffectiveGradient(X,Ae,Be,e0e,e1e,f,c[Liq-1])

def cutC(A,C):
    
    val = np.max(np.abs(A))
    
    for i in range(np.size(C) - 1):
        if (C[i] > val):
            C[i] = C[i-1]
    
    return C


def NeglectNanAverage(A):
    
    temp = 0
    for i in A: 
        if(i>0 or i<0):
            temp = temp + i
    
    return temp
    
def FreqEfficiency(Liq,LeftRight,d,c,p,bf,ef,sf,Piezo,count,TO):
    
    F = np.linspace(bf,ef,sf)
    A = np.zeros(np.size(F))
    B = np.zeros(np.size(F))
    C = np.zeros(np.size(F))
    

    errorcountf  = 0
    for i in range(np.size(F)):
        
        PU = GetBorderSpecs(F[i],d,c,p,0,Piezo)
        pg = GetGeneratedWave(Piezo,p,c,d,F[i],PU)
        
        AmpPhase, ABee, error = GetWaves(PU,d,c,F[i],p,Piezo,pg,TO,Liq)
        
        amp = cm.polar(pg)[0]
        
        
        A[i] = gradient(ABee,Liq,LeftRight,d,F[i])/amp
        B[i] = gradient(ABee,Liq,LeftRight,d,F[i])/F[i]
        C[i] = amp*ef/20000.
        
        
        if(count):
            print i+1 , '/', sf
        
        if(error and not count):
            print 'At ',F[i],'Hz'
            errorcountf = errorcountf + 1
            A[i] = 0
    
    if(cut):
        C = cutC(A,C)
    
    if(count):
        B = B/(2*(np.abs(NeglectNanAverage(B)/NeglectNanAverage(C))) )     
        plt.plot(F*10**-6,A, label = "Driving efficiency")
        plt.plot(F*10**-6,B, label = "Rescaled gradient $(Hz^{-1})$", color = 'g')
        plt.plot(F*10**-6,C, label = "Indication of pressure generated in transducer $(m \cdot Hz)$", color = 'r')
        

        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.title("Driving efficiency as function of frequency")
        plt.xlabel("Frequency $(MHz)$")
        plt.ylabel("Generated pressure to effective gradient scaling parameter $(m^{-1})$")
        plt.legend(bbox_to_anchor=(LegendX, LegendY), loc=2, borderaxespad=0.)
    
    
    
    maxindex,maxa = GetMaxIndex(A)
    minindex,mina = GetMinIndex(A)
    

    
    return F[maxindex],F[minindex],maxa,mina,errorcountf
    

def ImageWave(d,f,ABee,c,Piezo,pg,Liq):
    
    x = MakeRealXDistance(d)
    
    
    maxy = 0
    
            
    for i in range(np.size(d) - 2):
        tempx = np.linspace(x[i],x[i+1],samplerate)
           
            
        y = WF(tempx,ABee[i,0],ABee[i,1],ABee[i,2],ABee[i,3],f,c[i+1])

        if (i == Piezo - 2):
            
            for j in range(np.size(tempx)/2):
                
                y[j] = y[j] + np.real(pg * e(f,-tempx[j]+x[i],c[i+1],p[i+1]))
                y[np.size(y)-1-j] = y[np.size(y)-1-j] + np.real(pg * e(f,-tempx[j]+x[i],c[i+1],p[i+1]))
        
      
        plt.plot(tempx*10**6,y)
                
        if(np.max(np.abs(y)) > maxy):
            maxy = np.max(np.abs(y))
    
    for i in range(np.size(x)):
        
        
        plt.plot([x[i]*10**6,x[i]*10**6],[-maxy,maxy], lw = 2)

    xp = (x[Piezo-2] + x[Piezo-1])/2.
    
    plt.plot([xp*10**6,xp*10**6],[-maxy,maxy], color = 'black', linestyle = "--")
    
    plt.plot([-np.max(x)*10**6/11., np.max(x)*10**6 + np.max(x)*10**6/11.],[0,0])
     
    
     
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))    
    plt.xlim(-np.max(x)*10**6/7. , np.max(x)*10**6 + np.max(x)*10**6/7.)
    plt.ylabel('Normalized acoustic pressure')
    plt.xlabel('Distance $(\mu m)$')
    plt.title('TM modelled acoustic pressure profile: f = ' + str(f/10.**6)+'MHz')
    
    
def ThicknessDep (bf,ef,bt,et,sf,st,d,p,c,Piezo,Liq,Layer,LeftRight):
    
    sizeL = d[Layer - 1] * 10**6
    errorcount = 0
    T = np.linspace(bt,et,st)
    MAXA = np.zeros(st)
    MINA = np.zeros(st)

    print ''
    print 'Frequencies with large errors are being displayed'
    print ''
    
    index = 0
    Bh = int(dt.datetime.now().strftime("%H"))
    Bm = int(dt.datetime.now().strftime("%M"))
    Bs = int(dt.datetime.now().strftime("%S"))    
    Bd = int(dt.datetime.now().strftime("%d"))  
    
    for t in T:
        d[Layer - 1] = t
        fmax, fmin, maxa, mina,errorcountf = FreqEfficiency(Liq,LeftRight,d,c,p,bf,ef,sf,Piezo,False,True)
        
        if(maxa > 0):
            MAXA[index] = maxa
        else:
            MAXA[index] = 0
        
        if(mina < 0):
            MINA[index] = mina
        else:
            MINA[index] = 0

        Nh = int(dt.datetime.now().strftime("%H"))
        Nm = int(dt.datetime.now().strftime("%M"))
        Ns = int(dt.datetime.now().strftime("%S"))
        Nd = int(dt.datetime.now().strftime("%d"))    
        
        s = Ns- Bs
        
        if(s<0):
            s = 60+s
            Nm = Nm - 1
            
        m = Nm- Bm
        
        if(m<0):
            m = 60+m
            Nh = Nh - 1

        h = Nh - Bh
        
        if(h<0):
            h = 24+h
            Nd = Nd - 1
        
        day = Nd - Bd 
            
        index = index + 1 
        print index, '/',st, "    // Time running:", day, "d", h, "h",m, "m", s,'s'
        errorcount= errorcount+errorcountf
        
    print ''
    print '$$$Error Report$$$'
    print 'Amount of major fitting errors =', errorcount, '/',sf*st
    print 'In ', st, ' thickness samples with', sf, 'frequency samples each'
    print ''
    print 'Gradients found at large errors are overwritten with zero!'
    
    
    inMax, valMax = GetMaxIndex(MAXA)
    inMin, valMin = GetMinIndex(MINA)    
    
    print ''
    print 'Maximal positive gradient driving efficiency found: ', valMax
    print 'At', T[inMax]*10**6, 'um'
    print ''
    print 'Maximal negative gradient driving efficiency found: ', valMin
    print 'At', T[inMin]*10**6, 'um'
    
    
    plt.plot(T*10**6,MAXA, label = 'Positive gradient driving efficiency')
    plt.plot(T*10**6,MINA, label = 'Negative gradient driving efficiency')
    
    plt.plot([bt*10**6,et*10**6],[0,0], color = 'black')
    
    plt.plot([sizeL,sizeL],[np.min(MINA),np.max(MAXA)], color = 'black', ls = '--', label = 'Specified original thickness')
    
    temptitle = 'Response between '+str(bf/10**6)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title(temptitle+' MHz and '+str(ef/10**6)+' MHz in ' +str(sf) + ' linear samples')
    plt.xlabel("Layer thickness $(\mu m)$")
    plt.ylabel("Generated pressure to effective gradient scaling parameter $(m^{-1})$")
    plt.legend(bbox_to_anchor=(LegendX, LegendY), loc=2, borderaxespad=0.)
    


def DrawFreqEff(Liq,LeftRight,d,c,p,bf,ef,sf,Piezo):

    
    fmax,fmin,maxa,mina,error = FreqEfficiency(Liq,LeftRight,d,c,p,bf,ef,sf,Piezo,True,False)
    
         

    print ''
    print 'Best frequency for negative gradient:'
    if(mina < 0):
        print fmin/10**6, 'MHz'
        
    else:
        print 'Not found!'   
        
    print ''
    print 'Best frequency for positive gradient:'
    if(maxa > 0):
        print fmax/10**6, 'MHz'
        
    else:
        print 'Not found!'   

def DrawWave(Liq,d,c,p,Piezo,f,t,standing):

    if(standing and fitCheck):
        t = [0,df/1.,df/2.,df/3.,df/4.,df/5.,df/6.,df/7.,df/8.,df/9.,df/10.]
        
    elif(standing):
        t = [df/6.]
    
    for tn in t:
        
        PU = GetBorderSpecs(f,d,c,p,tn,Piezo)
        
        pg = GetGeneratedWave(Piezo,p,c,d,f,PU)
        
        AmpPhase,ABee,error =  GetWaves(PU,d,c,f,p,Piezo,pg,False,Liq)
        

        
        if(standing == True):
            MakeImageStanding(d,f,AmpPhase,c)
            
        else:
            ImageWave(d,f,ABee,c,Piezo,pg,Liq)
    
    
def grabnumber(num,sdec,dec):
    

    num = num*(10**(sdec))
    
    if(sdec == dec):
        num= np.round(num,0)
        
    num = int(num)
    return num%10
    
def numberarray(num, dec):
        store = ["" for x in range(dec + 1)]

      
        
        for m in range(dec+1):
            store[m]=str(grabnumber(num,m,dec))
        
        return store

def getlastnumarray(store):
    temp = 0
    for l in range(np.size(store)):
        if(store[np.size(store) - l-1] == '0'):
            temp = temp + 1
        
        else:
            return (np.size(store) - temp)
    
    return np.size(store)
    
           
def writetofile(direc,name,T,dec):
    
    filename = direc + name    
    
    FW = open(filename,'w')
    
    temp = np.round(T,dec)
    
    
    for i in range(np.size(temp[0])):
        for j in range(np.size(temp[0])):
            
            if(j != i):
                if(temp[i][j] < 0):
                    FW.write('-')
                    temp[i][j] = temp[i][j]*-1
                    
                array = numberarray(temp[i][j],dec)
                
                
                if(i == 0 or j == 0 or deletezeros):
                    index = getlastnumarray(array)
                else:
                    index = np.size(array)
                
                for k in range(index):
                    
                    FW.write(array[k])
                    if(index>1 and k==0):
                        FW.write('.')
            
            else:
                if(i == 0):
                    FW.write('To From')
                    
                else:
                    for o in range(int((dec+1)*1.8)):
                        FW.write('-')
            
            
            if(j < np.size(temp[0])-1):
                    FW.write('&')
        FW.write('\\\\ \n')

def LayerToLayer(T, RT):
    for i in range(np.size(T[0])-1):
        for j in range(np.size(T[0])-1):
            if(i != j and j>i):
                temp1 = j - i
                temp2 = 1.
                for k in range(temp1):
                    temp2 = temp2*RT[i+k,0]
                T[i+1,j+1] = temp2
    return T
    
def MakeRTArray(RT):
    temp = np.size(RT[:][:])/2.
    T = np.zeros((temp+2)**2).reshape([(temp+2),(temp+2)])

    for i in range(np.size(T[0])):
        T[0,i] = i
        T[i,0] = i
        T[i,i] = 1
    return T

def MakeBackToFrontRT_T(temp):
    for i in range((np.size(temp[:][:]))/2):
            temp[i][0] = 2 - temp[i][0]
    return temp

def MakeSuperT(RT):
    T = MakeRTArray(RT)
    BF_RT_T = MakeBackToFrontRT_T(RT)   
    T = LayerToLayer(T,BF_RT_T)    
    T = np.transpose(T)    
    RT = MakeBackToFrontRT_T(BF_RT_T)    
    T = LayerToLayer(T,RT)
    return T
           
def makeRT(z_0):
    RT = np.zeros((np.size(z_0)-1)*2).reshape([(np.size(z_0)-1),2])


    for i in range(np.size(z_0)-1):
        z_00 = z_0[i]
        z_01 = z_0[i+1]     
        RT[i,0] = 2*z_01/(z_00+z_01)
        RT[i,1] = -(1-RT[i,0])
    return RT

def AddReflection(T,RT):
    
    NoI = np.size(T[0])-1
    
    for i in range(NoI):
        
        for j in range(NoI-1):
            if(j<(NoI-i-1)):
                
                
                T[i+1,NoI-j] = T[i+1,NoI-j-1] * RT[NoI-j-2,1]
    
    for i in range(NoI-1):
        
        for j in range(NoI-1):
            if(j<(NoI-i-1)):
                
                #print i
                T[NoI-i,j+1] = T[NoI-i,j+2] * -RT[j,1]
    
    return T
    

def SaveCofTable(c,p):
    
    direc = MainWorkspace + Flowcellname+"\\"

    if not os.path.exists(direc):
        os.makedirs(direc)


    nameT = "T coefficients of "+ Flowcellname +'.txt'      
    nameR = "R coefficients of "+ Flowcellname +'.txt'  
    
    z_0 = c*p
    RT = makeRT(z_0) 
    
    
    T = MakeSuperT(RT)     
    
    
    writetofile(direc,nameT,T,dec)

            
    R = AddReflection(T,RT)
            
    writetofile(direc,nameR,R,dec)    
    
def main():
    
    plt.close()
    
    bool1 = (FreeDraw and FrequencyPerformance)
    bool2 = (FreeDraw and LayerThicknessPerformance)
    bool3 = (LayerThicknessPerformance and FrequencyPerformance)
    bool4 = (FreeDraw and SaveCoefficients)
    bool5 = (FrequencyPerformance and SaveCoefficients)
    bool6 = (LayerThicknessPerformance and SaveCoefficients)

    
    if(bool1 or bool2 or bool3 or bool4 or bool5 or bool6):
        print ''
        print 'Please choose only one mode, FreeDraw, FrequencyPerfomance, LayerThicknessPerformance or SaveCoefficents'
        return 0
        
    if(LeftRight != "Right" and LeftRight != "Left"):
        print "LeftRight must be either 'Left' or 'Right'"
        return 0
    
    if(np.size(d) != np.size(c) or np.size(d) != np.size(p)):
        print "The Flow Cell layer parameter arrays need to be of equal size" 
        return 0         
      
      
    if(LayerThicknessPerformance):
        ThicknessDep (bf,ef,bt,et,sf,st,d,p,c,Piezo,Liq,Layer,LeftRight)

    if(FrequencyPerformance):    
        DrawFreqEff(Liq,LeftRight,d,c,p,bf,ef,sf,Piezo)

    if(FreeDraw):    
        DrawWave(Liq,d,c,p,Piezo,f,t,standing)
    
    if(SaveCoefficients):   
        SaveCofTable(c,p)
    
    print ''
    print 'Done!'
    print ''
    
    
    
d = np.array(thickness) * 10**(-6)
c = np.array(speed)
p = np.array(density) * 10**3
    
bf=bf *10**6
ef=ef *10**6    

bt = bt*10**-6
et = et*10**-6 

f = f * 10**6  
main()    



