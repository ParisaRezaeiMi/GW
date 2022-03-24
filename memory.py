# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sxs 
import scipy
import matplotlib.pyplot as plt
import numpy as np
from sympy.physics.wigner import wigner_3j
from scipy.integrate import cumtrapz

C=3.0*10**8
SM=1.989*10**30
Mpc=3.086*10**22
R=100*Mpc
M=100*SM
GG=6.67*10**(-11)

## the following function was not useful! check michael's code on utilsss
def ylm(l,m,theta,phi):
    ylmm = scipy.special.sph_harm(m, l, theta, phi) 
    
    return (ylmm) 


##defining G, which is an angular integral 
def G(l1,l2,l3,m1,m2,m3):
    g = (-1.)**(m1+m2)*np.sqrt(((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.))/(4.*np.pi)) \
                *wigner_3j(l1,l2,l3,0,2,-2)*wigner_3j(l1,l2,l3,-m1,m2,-m3)
                
    return (float(g)) 
##w is h^lm in the paper, loaded from the data

w = sxs.load("SXS:BBH:0123/Lev/rhOverM", extrapolation_order=2)

#plt.plot(w.t, w.data.view(float))

def Wgen(l,m):
    hlm=w[:, w.index(l,m)]
    return hlm
    
def multLM(l2,l3,m2,m3):
    multiplication=np.gradient(Wgen(l2,m2),w.t)*np.conjugate(np.gradient(Wgen(l3,m3),w.t))
    return multiplication
def integralLM(l2,l3,m2,m3):
    integ=cumtrapz(multLM(l2,l3,m2,m3), w.t)
    return integ
def h20memoryLM(l2,l3,m2,m3):
    ans=1/24**0.5 *G(2,l2,l3,0,m2,m3)* integralLM(l2,l3,m2,m3)
    return ans

#used for sYlm
def fac(n):
   result = 1
   for i in range(2, n+1):
      result = result *i
   return result
#used for sYlm
def dlms(l, m, s, Theta):
    sq = np.sqrt(fac(l+m)*fac(l-m)*fac(l+s)*fac(l-s))
    d = 0.
    for k in range(max(0,m-s),min(l+m,l-s)+1):
        d = d + (-1.)**k*np.sin(Theta/2.)**(2.*k+s-m)*np.cos(Theta/2.)**(2.*l+m-s-2.*k)/(fac(k)*fac(l+m-k)*fac(l-s-k)*fac(s-m+k))
    return sq*d

def sYlm(s,l,m,Theta,Phi):
    
    res = (-1.)**(-s)*np.sqrt((2.*l+1)/(4*np.pi))*dlms(l,m,-s,Theta)
    if res==0:
        return 0.
    else:
        return complex(res*np.cos(m*Phi), res*np.sin(m*Phi))

def hlmmemory(l1,m1):
    hlmmemory=0
    for l2 in range(2, 4):
        for l3 in range(2, 4):
            m2=-l2
            while m2 <= l2:
                m2 += 1
                m3=-l3
                while m3 <= l3:
                    hlmmemory=hlmmemory+ 1/24**0.5 *G(l1,l2,l3,m1,m2,m3)*integralLM(l2, l3, m2, m3)
                    m3 += 1
    return hlmmemory


                
hlm=0
totalhlm=0
hlmmemoryonly=0
hlmmemoryWoutspins=0
for l1 in range(2, 4):
        # print('l3',l3)
        m1=-l1
        while m1 <= l1:
            m1 += 1
            totalhlm=totalhlm+(Wgen(l1,m1)[:-1]+hlmmemory(l1,m1))*sYlm(-2,l1,m1,np.pi/2,0)
            hlm=hlm+Wgen(l1,m1)[:-1]*sYlm(-2,l1,m1,np.pi/2,0)
            hlmmemoryonly=hlmmemoryonly+hlmmemory(l1,m1)*sYlm(-2,l1,m1,np.pi/2,0)
            hlmmemoryWoutspins=hlmmemoryWoutspins+hlmmemory(l1,m1)


temp=GG*M/R/C**2.0
print('scale factor',temp)     



scaledhlmtotal= (GG*M/R/C**2.0)*totalhlm
scaledhlmmemory= (GG*M/R/C**2.0)*hlmmemoryonly


print('***********')          

print('hlm scaled',scaledhlmtotal)     
print('***********')          
     
print('hlm only memory',hlmmemoryonly) 
print('***********')          
         



#plt.plot(w.t[:-1],totalhlm-hlm)

plt.plot(GG*M/C**3.0*w.t[:-1],scaledhlmtotal)

#plt.plot(w.t[:-1],scaledhlmtotal)  papeeeer  1906.06263
#plt.plot(GG*M/C**3.0*w.t[:-1],scaledhlmmemory)

plt.show()



# =============================================================================
# COMMENTS:
# 

     
# =============================================================================
# h20memoryef=0
# for l2 in range(2, 6):
#     for l3 in range(2, 6):
#         # print('l3',l3)
#         m2=-l2
#         while m2 <= l2:
#             m2 += 1
#             m3=-l3
#             while m3 <= l3:
#                 h20memoryef=h20memoryef+ R/C/24**0.5 *G(2,l2,l3,0,m2,m3)*integralLM(l2, l3, m2, m3)
#                 m3 += 1
#            
# =============================================================================
                
                
# 
# w22= w[:, w.index(2,2)]
# w22dot=np.gradient(w22,w.t)
# w22dotc=np.conjugate(w22dot)
# mult22=w22dotc*w22dot
# 
# 
# ##the integral of the gradient of modes over time 
# #integral=np.cumsum(mult)*w.t
# integral = cumtrapz(mult22, w.t)
# #print('checkk:',integralLM(2,2)-integral)
# 
# #computing the memory effect of h22 which is the dominent mode
# h20memoryof22= R/C/24**0.5 *G(2,2,2,0,2,2)* integral
# h20memoryof22prime= h20memoryLM(2,2,2,2)
# #print('checkkk:',h20memoryof22prime-h20memoryof22)
# #?????? why w.t don't work in the following plot? what integration changes in the dimensions?
# #plt.plot(w.t[:-1],h20memoryof22)  
# #plt.show()
# 
# 
# #trying out different l & m's, to see which combinations contribute with the memory
# #first, we keep l' & l'' as 2. We change m' & m'', we get zero for G if they're not similar
# #print(' l=2')
# 
# # =============================================================================
# # print('G2,2=',G(2,2,2,0,2,2)) 
# # print('G2,-2=',G(2,2,2,0,2,-2))
# # print('G-2,-2=',G(2,2,2,0,-2,-2))
# # print('G-2,2=',G(2,2,2,0,-2,2))
# # print('G1,1=',G(2,2,2,0,1,1))
# # print('G1,-1=',G(2,2,2,0,1,-1))
# # print('G-1,1=',G(2,2,2,0,-1,1))
# # print('G-1,-1=',G(2,2,2,0,-1,-1))
# # print('G0,1=',G(2,2,2,0,0,1))
# # print('G-1,0=',G(2,2,2,0,-1,0))
# # print('G0,0=',G(2,2,2,0,0,0))
# # =============================================================================
# 
# 
# # =============================================================================
# # G2,2= 0.18022375157286857
# # G2,-2= 0.0
# # G-2,-2= 0.18022375157286857
# # G-2,2= 0.0
# # G1,1= -0.09011187578643429
# # G1,-1= 0.0
# # G-1,1= 0.0
# # G-1,-1= -0.09011187578643429
# # G0,1= 0.0
# # G-1,0= 0.0
# # G0,0= -0.18022375157286857
# # =============================================================================
# 
# #print(' l=1')
# 
# #first, we keep l' & l'' as 1. We change m' & m'', we get zero for all G.
# # =============================================================================
# # print('G1,-1=',G(2,1,1,0,1,-1))
# # print('G-1,-1=',G(2,1,1,0,-1,-1))
# # print('G-1,1=',G(2,1,1,0,-1,1))
# # print('G1,1=',G(2,1,1,0,1,1))
# # print('G0,0=',G(2,1,1,0,0,0))
# # print('G0,1=',G(2,1,1,0,0,1))
# # print('G-1,0=',G(2,1,1,0,-1,0))
# # =============================================================================
# 
# 
# # =============================================================================
# # G1,-1= -0.0
# # G-1,-1= 0.0
# # G-1,1= -0.0
# # G1,1= 0.0
# # G0,0= 0.0
# # G0,1= 0.0
# # G-1,0= -0.0
# # G0,0= 0.0
# # =============================================================================
# 
# 
# 
# 
# w21= w[:, w.index(2,1)]
# #print('w21=',w21)
# w21dot=np.gradient(w21,w.t)
# w21dotc=np.conjugate(w21dot)
# mult21=w21dotc*w21dot
#  
# w2n1= w[:, w.index(2,-1)]
# #print('w2n1=',w2n1)
# w2n1dot=np.gradient(w2n1,w.t)
# w2n1dotc=np.conjugate(w2n1dot)
# mult2n1=w2n1dotc*w2n1dot
# 
# 
# w20= w[:, w.index(2,0)]
# #print('w20=',w20)
# w20dot=np.gradient(w20,w.t)
# w20dotc=np.conjugate(w20dot)
# mult20=w20dotc*w20dot
# 
# 
# # =============================================================================
# # print('mult21=',mult21)
# # print('mult2n1=',mult2n1)
# # print('mult20=',mult20)
# # 
# # =============================================================================
# 
# 
# # #the integral of the gradient of modes over time 
# # #integral=np.cumsum(mult)*w.t
# integral21 = cumtrapz(mult21, w.t)
# integral2n1 = cumtrapz(mult2n1, w.t)
# integral20 = cumtrapz(mult20, w.t)
# 
# #computing the memory effect of h22 which is the dominent mode
# #h20memoryTot= R/C/24**0.5 * ( G(2,2,2,0,1,1)* integral21+ G(2,2,2,0,-1,-1)* integral2n1+G(2,2,2,0,0,0)* integral20
#  
# 
# h20memory21=R/C/24**0.5 * G(2,2,2,0,1,1)* integral21
# h20memory2n1=R/C/24**0.5 * G(2,2,2,0,-1,-1)* integral2n1
# h20memory20=R/C/24**0.5 * G(2,2,2,0,0,0)* integral20
# 
# h20memorytot=h20memoryof22+h20memory21+h20memory2n1+h20memory20
# #plt.plot(w.t[:-1],h20memorytot)
# #plt.plot(w.t[:-1],h20memoryof22)
# #print('h20', h20memoryof22)
# 
# #plt.show()
# 
# #ArithmeticErrorplt.plot(w.t[:-1],h20memoryef)
# 
# #####plt.plot(w.t[:-1],totalhlm)
# 
# ####plt.show()
# 
# 
# 
# 
# #print('ylm=',ylm(2,0,0,np.pi/3))
# 
# #eq.1 of micheal's paper, only for lm=20
# #h=(w20+h20memorytot)*ylm(2,0,0,np.pi/3)
# 
# #plt.plot(w.t,h)
# #plt.show()
# 
# =============================================================================
