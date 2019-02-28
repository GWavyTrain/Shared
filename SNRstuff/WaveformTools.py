import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import interpolate

import PhenomA as pa
import LISA as li

""" Constants """
C       = 299792458.         # m/s
YEAR    = 3.15581497632e7    # sec
TSUN    = 4.92549232189886339689643862e-6 # mass of sun in seconds (G=C=1)
MPC     = 3.08568025e22/C    # mega-Parsec in seconds

TOBS_MAX = 4*YEAR # Maximum observation period (LISA's nominal mission lifetime)

""" Cosmological values """
H0      = 69.6      # Hubble parameter today
Omega_m = 0.286     # density parameter of matter

def get_Dl(z, Omega_m, H0):
    """ calculate luminosity distance in geometric units """
    # see http://arxiv.org/pdf/1111.6396v1.pdf
    x0 = (1. - Omega_m)/Omega_m
    xZ = x0/(1. + z)**3

    Phi0  = (1. + 1.320*x0 + 0.4415*x0**2  + 0.02656*x0**3)
    Phi0 /= (1. + 1.392*x0 + 0.5121*x0**2  + 0.03944*x0**3)
    PhiZ  = (1. + 1.320*xZ + 0.4415*xZ**2  + 0.02656*xZ**3)
    PhiZ /= (1. + 1.392*xZ + 0.5121*xZ**2  + 0.03944*xZ**3)
    
    return 2.*C/H0*(1.0e-3*MPC)*(1. + z)/np.sqrt(Omega_m)*(Phi0 - PhiZ/np.sqrt(1. + z))


def get_z(z, Dl, Omega_m, H0):
    """ calculate redishift uisng root finder """
    
    return get_Dl(z, Omega_m, H0) - Dl
       
def SetFreqBounds(self, lisa):
    """ """
    
    Mc = self.M*self.eta**(3./5) # chirp mass

    # Determine start frequency of binary
    if (self.f_start == None): # T_merge was specified
        self.f_start = (5.*Mc/self.T_merge)**(3./8.)/(8.*np.pi*Mc)
    else:
        self.T_merge = 5.*Mc/(8.*np.pi*self.f_start*Mc)**(8./3.)

    # Determine the end frequency
    if (self.T_merge > lisa.Tobs):
        self.f_end = (5.*Mc/(np.abs(lisa.Tobs-self.T_merge)))**(3./8.)/(8.*np.pi*Mc)
    else:
        self.f_end = pa.get_freq(self.M, self.eta, "cut") # PhenomA cut-off frequency i.e. frequency upper bound    

    return
    
def calc_k(theta, phi):
    """ Calculate the unit-direction vector pointing towards the source """
    
    sth = np.sin(theta)

    k = -np.array([sth*np.cos(phi), sth*np.sin(phi), np.cos(theta)])

    return k
    
def calc_k_dot_r(k, rij):
    """ Dot product between unit-direction vector and the S/C unit-separation vectors """

    k_dot_r = k[0]*rij[0,:,:,:] + k[1]*rij[1,:,:,:] + k[2]*rij[2,:,:,:]

    return k_dot_r
    
def CalcStrain(self, lisa, theta=None, phi=None, iota=None): 
    """ """
    
#def CalcStrain(self, lisa, Delta_logf, theta=None, phi=None, iota=None):
    """ """
    Delta_logf = np.log(self.f_end) - np.log(self.f_start) 
    
    if (Delta_logf > 0.00005): # Generate a track   
        N = 500 # number of points
        f = np.logspace(np.log10(self.f_start), np.log10(self.f_end), N)
        
        if (theta == None and phi == None): # generate sky averaged response
            self.Figure_Type = 'track'
            X_char = np.sqrt(16./5*f)*pa.Aeff(f, self.M, self.eta, self.Dl)
            
        else: # Generate X Michelson channel
            self.Figure_Type = 'track_sky_dependent'
            Aeff = pa.Aeff(f, self.M, self.eta, self.Dl)
            
            # calculate LISA transfer function
            tStar = pa.dPsieff_df(f, self.M, self.eta, 0.0)/(2*np.pi) # stationary time of SPA
            
            k   = calc_k(theta, phi)
            x   = lisa.SC_Orbits(tStar)
            rij = lisa.SC_Seps(tStar, x)
            
            rij_OUTER_rij = rij.reshape((3,1, 3,3, N))*rij.reshape((1,3, 3,3, N))
            
            u = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)])
            v = np.array([np.sin(phi), -np.cos(phi), 0.0])
            ep = np.outer(u,u) - np.outer(v,v)
            ec = np.outer(u,v) + np.outer(v,u)
            
            dp12 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,0,1,:], ep)
            dc12 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,0,1,:], ec)
            dp21 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,1,0,:], ep)
            dc21 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,1,0,:], ec)
            dp13 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,0,2,:], ep)
            dc13 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,0,2,:], ec)
            dp31 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,2,0,:], ep)
            dc31 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,2,0,:], ec)
            
            kDOTrij = calc_k_dot_r(k, rij)
            
            kDOTr12 = kDOTrij[0,1,:] 
            kDOTr21 = kDOTrij[1,0,:]
            kDOTr13 = kDOTrij[0,2,:]
            kDOTr31 = kDOTrij[2,0,:]
            
            TransArg12 = f/(2*lisa.fstar)*(1. - kDOTr12)
            TransArg21 = f/(2*lisa.fstar)*(1. - kDOTr21)
            TransArg13 = f/(2*lisa.fstar)*(1. - kDOTr13)
            TransArg31 = f/(2*lisa.fstar)*(1. - kDOTr31)
            
            Trans12 = 0.5*np.sinc(TransArg12/np.pi)*np.exp(1j*TransArg12)*np.exp(1j*2*np.pi*f*np.dot(k,x[:,0,:])/C)
            Trans21 = 0.5*np.sinc(TransArg21/np.pi)*np.exp(1j*TransArg21)*np.exp(1j*2*np.pi*f*np.dot(k,x[:,1,:])/C)
            Trans13 = 0.5*np.sinc(TransArg13/np.pi)*np.exp(1j*TransArg13)*np.exp(1j*2*np.pi*f*np.dot(k,x[:,0,:])/C)
            Trans31 = 0.5*np.sinc(TransArg31/np.pi)*np.exp(1j*TransArg31)*np.exp(1j*2*np.pi*f*np.dot(k,x[:,2,:])/C)
            
            #########
#             Aeff_INTERP = interpolate.splrep(f, Aeff, s=0)
#             dPsieff_df_INTERP = interpolate.splrep(f, pa.dPsieff_df(f, self.M, self.eta, 0.0), s=0)
#             
#             epsilon = li.YEAR/10**6
#             t1_for = lisa.SC_Orbits(tStar + epsilon)[:,0,:]
#             t1_cen = lisa.SC_Orbits(tStar)[:,0,:]
#             t1_bac = lisa.SC_Orbits(tStar - epsilon)[:,0,:]
#             kDOTxdot1  = ( np.dot(k,t1_for) - np.dot(k,t1_bac) )/(2*epsilon)
#             kDOTxddot1 = ( np.dot(k,t1_for) - 2*np.dot(k,t1_cen) + np.dot(k,t1_bac) )/(2*epsilon)**2
#             
#             t2_for = lisa.SC_Orbits(tStar + epsilon)[:,1,:]
#             t2_cen = lisa.SC_Orbits(tStar)[:,1,:]
#             t2_bac = lisa.SC_Orbits(tStar - epsilon)[:,1,:]
#             kDOTxdot2  = ( np.dot(k,t2_for) - np.dot(k,t2_bac) )/(2*epsilon)
#             kDOTxddot2 = ( np.dot(k,t2_for) - 2*np.dot(k,t2_cen) + np.dot(k,t2_bac) )/(2*epsilon)**2
#             
#             t3_for = lisa.SC_Orbits(tStar + epsilon)[:,2,:]
#             t3_cen = lisa.SC_Orbits(tStar)[:,2,:]
#             t3_bac = lisa.SC_Orbits(tStar - epsilon)[:,2,:]
#             kDOTxdot3  = ( np.dot(k,t3_for) - np.dot(k,t3_bac) )/(2*epsilon)
#             kDOTxddot3 = ( np.dot(k,t3_for) - 2*np.dot(k,t3_cen) + np.dot(k,t3_bac) )/(2*epsilon)**2
#             
#             TransArg12 *= (1 + kDOTxdot1/C)
#             TransArg21 *= (1 + kDOTxdot2/C)
#             TransArg13 *= (1 + kDOTxdot1/C)
#             TransArg31 *= (1 + kDOTxdot3/C)
#             
#             Aeff_Pert_factor1 = (Aeff + f*interpolate.splev(f, Aeff_INTERP, der=1))/C
#             Aeff_Pert_factor2 = f*Aeff/(4*np.pi)*interpolate.splev(f, dPsieff_df_INTERP, der=1)/C
            ##########
            
            # iota = pi/2, psi = 0
#             y12_a = Trans12*(Aeff + Aeff_Pert_factor1*kDOTxdot1 + Aeff_Pert_factor1*kDOTxddot1)*0.5*dp12
#             y21_a = Trans21*(Aeff + Aeff_Pert_factor1*kDOTxdot2 + Aeff_Pert_factor1*kDOTxddot2)*0.5*dp21
#             y13_a = Trans13*(Aeff + Aeff_Pert_factor1*kDOTxdot1 + Aeff_Pert_factor1*kDOTxddot1)*0.5*dp13
#             y31_a = Trans31*(Aeff + Aeff_Pert_factor1*kDOTxdot3 + Aeff_Pert_factor1*kDOTxddot3)*0.5*dp31

#            original lines...
#            y12_a = 0.5*Trans12*Aeff*0.5*dp12
#            y21_a = 0.5*Trans21*Aeff*0.5*dp21
#            y13_a = 0.5*Trans13*Aeff*0.5*dp13
#            y31_a = 0.5*Trans31*Aeff*0.5*dp31

            # Hack for Kelly
            y12_a = 0.5*Trans12*Aeff*(0.5*(1 + np.cos(iota)**2)*dp12 + 1j*np.cos(iota)*dc12)
            y21_a = 0.5*Trans21*Aeff*(0.5*(1 + np.cos(iota)**2)*dp21 + 1j*np.cos(iota)*dc21)
            y13_a = 0.5*Trans13*Aeff*(0.5*(1 + np.cos(iota)**2)*dp13 + 1j*np.cos(iota)*dc13)
            y31_a = 0.5*Trans31*Aeff*(0.5*(1 + np.cos(iota)**2)*dp31 + 1j*np.cos(iota)*dc31)

            X_TDI = (y12_a - y13_a)*np.exp(-1j*f/lisa.fstar) + (y21_a - y31_a)
            
            XX_TDI = 0.5*np.abs(X_TDI)**2 # 8/5 -> 1/2 for only polarization average, for Kelly
            
            # iota = pi/2, psi = pi/4
#             y12_a = Trans12*(Aeff + Aeff_Pert_factor1*kDOTxdot1 + Aeff_Pert_factor1*kDOTxddot1)*0.5*dc12
#             y21_a = Trans21*(Aeff + Aeff_Pert_factor1*kDOTxdot2 + Aeff_Pert_factor1*kDOTxddot2)*0.5*dc21
#             y13_a = Trans13*(Aeff + Aeff_Pert_factor1*kDOTxdot1 + Aeff_Pert_factor1*kDOTxddot1)*0.5*dc13
#             y31_a = Trans31*(Aeff + Aeff_Pert_factor1*kDOTxdot3 + Aeff_Pert_factor1*kDOTxddot3)*0.5*dc31

#            y12_a = 0.5*Trans12*Aeff*0.5*dc12
#            y21_a = 0.5*Trans21*Aeff*0.5*dc21
#            y13_a = 0.5*Trans13*Aeff*0.5*dc13
#            y31_a = 0.5*Trans31*Aeff*0.5*dc31

            y12_a = 0.5*Trans12*Aeff*(-0.5*(1 + np.cos(iota)**2)*dc12 + 1j*np.cos(iota)*dp12)
            y21_a = 0.5*Trans21*Aeff*(-0.5*(1 + np.cos(iota)**2)*dc21 + 1j*np.cos(iota)*dp21)
            y13_a = 0.5*Trans13*Aeff*(-0.5*(1 + np.cos(iota)**2)*dc13 + 1j*np.cos(iota)*dp13)
            y31_a = 0.5*Trans31*Aeff*(-0.5*(1 + np.cos(iota)**2)*dc31 + 1j*np.cos(iota)*dp31)

            X_TDI = (y12_a - y13_a)*np.exp(-1j*f/lisa.fstar) + (y21_a - y31_a)
            
            XX_TDI += 0.5*np.abs(X_TDI)**2
            
            X_char = np.sqrt(4*f*XX_TDI)
            
    else:
    
        if (theta == None and phi == None): # generate sky averaged response
            N = 1
            f = np.array([self.f_start])
            
            self.Figure_Type = 'point'
            
            X_char = np.sqrt(16./5*pa.Aeff(f, self.M, self.eta, self.Dl)**2*np.sqrt(f)*(self.f_end - f))
            
        else:   
            N = 1
            f = np.array([self.f_start])

            self.Figure_Type = 'point_sky_dependent'
            Aeff = pa.Aeff(f, self.M, self.eta, self.Dl)
        
            # calculate LISA transfer function
            tStar = pa.dPsieff_df(f, self.M, self.eta, 0.0)/(2*np.pi) # stationary time of SPA
        
            k   = calc_k(theta, phi)
            x   = lisa.SC_Orbits(tStar)
            rij = lisa.SC_Seps(tStar, x)
        
            rij_OUTER_rij = rij.reshape((3,1, 3,3, N))*rij.reshape((1,3, 3,3, N))
        
            u = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)])
            v = np.array([np.sin(phi), -np.cos(phi), 0.0])
            ep = np.outer(u,u) - np.outer(v,v)
            ec = np.outer(u,v) + np.outer(v,u)
        
            dp12 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,0,1,:], ep)
            dc12 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,0,1,:], ec)
            dp21 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,1,0,:], ep)
            dc21 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,1,0,:], ec)
            dp13 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,0,2,:], ep)
            dc13 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,0,2,:], ec)
            dp31 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,2,0,:], ep)
            dc31 = np.einsum('nmk,nm->k', rij_OUTER_rij[:,:,2,0,:], ec)
        
            kDOTrij = calc_k_dot_r(k, rij)
        
            kDOTr12 = kDOTrij[0,1,:] 
            kDOTr21 = kDOTrij[1,0,:]
            kDOTr13 = kDOTrij[0,2,:]
            kDOTr31 = kDOTrij[2,0,:]
        
            TransArg12 = f/(2*lisa.fstar)*(1. - kDOTr12)
            TransArg21 = f/(2*lisa.fstar)*(1. - kDOTr21)
            TransArg13 = f/(2*lisa.fstar)*(1. - kDOTr13)
            TransArg31 = f/(2*lisa.fstar)*(1. - kDOTr31)
        
            Trans12 = 0.5*np.sinc(TransArg12/np.pi)*np.exp(1j*TransArg12)*np.exp(1j*2*np.pi*f*np.dot(k,x[:,0,:])/C)
            Trans21 = 0.5*np.sinc(TransArg21/np.pi)*np.exp(1j*TransArg21)*np.exp(1j*2*np.pi*f*np.dot(k,x[:,1,:])/C)
            Trans13 = 0.5*np.sinc(TransArg13/np.pi)*np.exp(1j*TransArg13)*np.exp(1j*2*np.pi*f*np.dot(k,x[:,0,:])/C)
            Trans31 = 0.5*np.sinc(TransArg31/np.pi)*np.exp(1j*TransArg31)*np.exp(1j*2*np.pi*f*np.dot(k,x[:,2,:])/C)
#            Trans12 = 0.5*np.sinc(TransArg12/np.pi)*np.exp(1j*TransArg12)
#            Trans21 = 0.5*np.sinc(TransArg21/np.pi)*np.exp(1j*TransArg21)
#            Trans13 = 0.5*np.sinc(TransArg13/np.pi)*np.exp(1j*TransArg13)
#            Trans31 = 0.5*np.sinc(TransArg31/np.pi)*np.exp(1j*TransArg31)

#            y12_a = Trans12*Aeff*0.5*dp12/2
#            y21_a = Trans21*Aeff*0.5*dp21/2
#            y13_a = Trans13*Aeff*0.5*dp13/2
#            y31_a = Trans31*Aeff*0.5*dp31/2

            y12_a = 0.5*Trans12*Aeff*(0.5*(1 + np.cos(iota)**2)*dp12 + 1j*np.cos(iota)*dc12)
            y21_a = 0.5*Trans21*Aeff*(0.5*(1 + np.cos(iota)**2)*dp21 + 1j*np.cos(iota)*dc21)
            y13_a = 0.5*Trans13*Aeff*(0.5*(1 + np.cos(iota)**2)*dp13 + 1j*np.cos(iota)*dc13)
            y31_a = 0.5*Trans31*Aeff*(0.5*(1 + np.cos(iota)**2)*dp31 + 1j*np.cos(iota)*dc31)
            
            X_TDI = (y12_a - y13_a)*np.exp(-1j*f/lisa.fstar) + (y21_a - y31_a)
        
        #XX_TDI = 8./5*np.abs(X_TDI)**2
            XX_TDI = 0.5*np.abs(X_TDI)**2
        
#            y12_a = Trans12*Aeff*0.5*dc12/2
#            y21_a = Trans21*Aeff*0.5*dc21/2
#            y13_a = Trans13*Aeff*0.5*dc13/2
#            y31_a = Trans31*Aeff*0.5*dc31/2

            y12_a = 0.5*Trans12*Aeff*(-0.5*(1 + np.cos(iota)**2)*dc12 + 1j*np.cos(iota)*dp12)
            y21_a = 0.5*Trans21*Aeff*(-0.5*(1 + np.cos(iota)**2)*dc21 + 1j*np.cos(iota)*dp21)
            y13_a = 0.5*Trans13*Aeff*(-0.5*(1 + np.cos(iota)**2)*dc13 + 1j*np.cos(iota)*dp13)
            y31_a = 0.5*Trans31*Aeff*(-0.5*(1 + np.cos(iota)**2)*dc31 + 1j*np.cos(iota)*dp31)

            X_TDI = (y12_a - y13_a)*np.exp(-1j*f/lisa.fstar) + (y21_a - y31_a)
        
#XX_TDI += 8./5*np.abs(X_TDI)**2
            XX_TDI += 0.5*np.abs(X_TDI)**2
        
            X_char = np.sqrt(4*XX_TDI*np.sqrt(f)*(self.f_end - f))
    
    return f, X_char
    
def CalcSNR(self, f, X_char, lisa):
    """ """

    if (self.Figure_Type == 'track'):
        N = len(f) # number of frequency samples

        d_logf       = np.log(f[1:]) - np.log(f[:N-1])
        
        term_i   = X_char[1:]**2/lisa.Sn(f[1:])
        term_im1 = X_char[:N-1]**2/lisa.Sn(f[:N-1])

        snrSQ = np.sum(0.5*(term_i + term_im1)*d_logf) 
        
    elif (self.Figure_Type == 'track_sky_dependent'):
        N = len(f) # number of frequency samples

        d_logf       = np.log(f[1:]) - np.log(f[:N-1])
        
        term_i   = X_char[1:]**2/lisa.PnC(f[1:])*lisa.NC
        term_im1 = X_char[:N-1]**2/lisa.PnC(f[:N-1])*lisa.NC
        
        snrSQ = np.sum(0.5*(term_i + term_im1)*d_logf)
        
    elif (self.Figure_Type == 'point'):

        snrSQ = (X_char**2/np.sqrt(f)/lisa.Sn(f))[0]
        
    elif (self.Figure_Type == 'point_sky_dependent'):

        snrSQ = (X_char**2/np.sqrt(f)/lisa.PnC(f)*lisa.NC)[0]
        
    return np.sqrt(snrSQ)
    
def PlotStrain(self, freqs, X_char, lisa):
    """ """

    if (self.Figure_Type == 'track'):
        fig, ax = plt.subplots(1, figsize=(8,6))
        plt.tight_layout()

        ax.set_xlabel(r'f [Hz]', fontsize=20, labelpad=10)
        ax.set_ylabel(r'Characteristic Strain', fontsize=20, labelpad=10)
        ax.tick_params(axis='both', which='major', labelsize=20)
        
        ax.set_xlim(1.0e-5, 1.0e0)
        ax.set_ylim(3.0e-22, 1.0e-15)

        ax.loglog(freqs, np.sqrt(freqs)*X_char) 

        f = np.logspace(np.log10(1.0e-5), np.log10(1.0e0), 1000)
        ax.loglog(f, np.sqrt(f*lisa.Sn(f)))
        
    elif (self.Figure_Type == 'track_sky_dependent'):
        fig, ax = plt.subplots(1, figsize=(8,6))
        plt.tight_layout()

        ax.set_xlabel(r'f [Hz]', fontsize=20, labelpad=10)
        ax.set_ylabel(r'Characteristic Strain', fontsize=20, labelpad=10)
        ax.tick_params(axis='both', which='major', labelsize=20)
        
        ax.set_xlim(1.0e-5, 1.0e0)
        ax.set_ylim(3.0e-22, 1.0e-15)

        ax.loglog(freqs, np.sqrt(freqs)*X_char) 

        f = np.logspace(np.log10(1.0e-5), np.log10(1.0e0), 1000)
        ax.loglog(f, np.sqrt(f*lisa.PnC(f)))    
        
    elif (self.Figure_Type == 'point_sky_dependent'):
        fig, ax = plt.subplots(1, figsize=(8,6))
        plt.tight_layout()

        ax.set_xlabel(r'f [Hz]', fontsize=20, labelpad=10)
        ax.set_ylabel(r'Characteristic Strain', fontsize=20, labelpad=10)
        ax.tick_params(axis='both', which='major', labelsize=20)
        
        ax.set_xlim(1.0e-5, 1.0e0)
        ax.set_ylim(3.0e-22, 1.0e-15)

        ax.loglog(freqs, np.sqrt(freqs)*X_char, c='xkcd:bright orange', marker='o', markeredgecolor='xkcd:black', alpha=0.9, markersize=20) 
        ax.loglog(freqs, np.sqrt(freqs)*6.5e-19, c='xkcd:poo brown', marker='o', markeredgecolor='xkcd:black', alpha=0.6, markersize=20) 


        f = np.logspace(np.log10(1.0e-5), np.log10(1.0e0), 1000)
        ax.loglog(f, np.sqrt(f*lisa.PnC(f)), color='xkcd:chocolate')   

    elif (self.Figure_Type == 'point'):
        fig, ax = plt.subplots(1, figsize=(8,6))
        plt.tight_layout()

        ax.set_xlabel(r'f [Hz]', fontsize=20, labelpad=10)
        ax.set_ylabel(r'Characteristic Strain', fontsize=20, labelpad=10)
        ax.tick_params(axis='both', which='major', labelsize=20)
        
        ax.set_xlim(1.0e-5, 1.0e0)
        ax.set_ylim(3.0e-22, 1.0e-15)

        ax.loglog(freqs, np.sqrt(freqs)*X_char, c='xkcd:bright orange', marker='o', markeredgecolor='xkcd:black', alpha=0.8, markersize=20) 
        ax.loglog(freqs, np.sqrt(freqs)*6.5e-19, c='xkcd:mushroom', marker='o', markeredgecolor='xkcd:black', alpha=0.8, markersize=20) 


        f = np.logspace(np.log10(1.0e-5), np.log10(1.0e0), 1000)
        ax.loglog(f, np.sqrt(f*lisa.Sn(f)), color='xkcd:chocolate')    

    return


class Binary():
    """ 
    Binary Class
    -------------------------------------------
    Inputs:
        Specify source-frame masses: m1, m2
        Specify a distance parameter: z, Dl (redshift, luminosity distance IN SECONDS)
        Specify an initial condition parameter: T_merge, f_start
                    (note that an upper limit of 4 years will be set on the 
                     observation period)
    
    Methods:
        CalcStrain: Calculate the characteristic strain of the binary. If (the optional
                    arguments) sky angles are provided use the stataionary phase approximation
                    signal generator, else use PhenomA amplitude exclusively
                    
        CalcSNR: Calculate the SNR averaged over polarization, inclination,
                  and sky angles. Theta, phi (spherical polar) are optional arguments
                  allowing the user to calculate the SNR at a specific sky location
                  averaged over only polarization and inclination angles

        PlotStrain: Plot the characteristic strain
    
    """
    
    def __init__(self, m1, m2, z=None, Dl=None):
        # source-frame component masses
        self.m1 = m1
        self.m2 = m2
        
        # Store distance parameters
        if (Dl == None): # convert redshift into luminosity distance
            self.z = z # TODO: check that one of these is provided
            self.Dl = get_Dl(self.z, Omega_m, H0) # Dl returned in seconds (i.e. G=c=1, geometric units)
            print("Redshift provided. \n\tLuminosity Distance........... {} Mpc".format(self.Dl/MPC))

        else: # convert luminosity distance to redshift
            self.Dl = Dl # TODO: check that one of these is provided
            self.z = optimize.root(get_z, 1., args=(self.Dl, Omega_m, H0)).x[0]
            print("Luminosity Distance provided. \n\tredshift........... {}".format(self.z))
            
        # adjust source-frame masses to detector-frame masses
        self.m1 *= 1. + self.z 
        self.m2 *= 1. + self.z
        
        # calculate relevant mass parameters
        self.M   = self.m1 + self.m2 # total mass
        self.eta = self.m1*self.m2/self.M**2 # symmetric mass ratio
        
        self.f_start = None
        self.f_end   = None
        
        
    # Methods
    SetFreqBounds = SetFreqBounds
    CalcStrain    = CalcStrain
    CalcSNR       = CalcSNR
    PlotStrain    = PlotStrain








