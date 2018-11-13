# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import math

#Simple Lotka-Volterra

#dN/dt = r*N - a*N*P
#dP/dt = a*b*N*P

def dn_dt(N, P, r, a):
    return(N+(N*r-a*N*P))
def dp_dt(N, P, r, a, b):
    return(P+(a*N*P*b - a*N*P))

def sim_pop(N_init, P_init, r, a, b, nsteps):
    #Simulates the community for given parameters & number of steps
    track = np.zeros((nsteps, 2))
    N = N_init
    P = P_init
    for i in range(0, nsteps):
        track[i, 0] = N
        track[i, 1] = P
        N = dn_dt(N, P, r, a)
        P = dp_dt(track[i, 0], P, r, a, b)
    return(track)

#nsteps = 10000
#ratio = 100
#multiplier = 0.00005
#
#for r in (0.1, 1, 10):
#    for a in (0.1, 1, 10):
#        for b in (1, 10, 100):
#            p_init_vals = np.linspace(0, 10, num = 50)
#            peak_track = np.zeros((len(p_init_vals), 2))
#            for i in range(len(p_init_vals)):
#                track = sim_pop(N_init = ratio*p_init_vals[i], 
#                                P_init = p_init_vals[i], 
#                                r = r*multiplier, 
#                                a = a*multiplier, 
#                                b = b, nsteps = nsteps)
                
#                #Make plot showing this population w/ 2 y axes but shared x axis
#                #Define plot & axis 1 (for bacteria)
#                fig, ax1 = plt.subplots()
#                #Plot bacteria
#                ax1.plot(range(0, nsteps), track[:, 0], label = "N", color = 'blue')
#                #Set up bacterial y axis 
#                ax1.set_ylabel('N', color='b')
#                #Make the phage x axis a twin of the bacterial x axis
#                ax2 = ax1.twinx()
#                #Plot phage
#                ax2.plot(range(0, nsteps), track[:, 1], label = "P", color = 'red')
#                #Set up phage y axis
#                ax2.set_ylabel('P', color='r')
#                plt.title("r="+str(r)+", a="+str(a)+", b="+str(b)+", p_init="+
#                          str(p_init_vals[i]))
#                plt.show()
                
                
#                if max(track[:, 0]) == ratio*p_init_vals[i]:
#                    peak_track[i, 0] = 0 #the bact pop never grew
#                else:
#                    peak_track[i, 0] = max(track[:, 0])
#            
#            plt.plot(p_init_vals, peak_track[:, 0])
#            plt.title("r="+str(r)+", a="+str(a)+", b="+str(b))
#            plt.show()


#Modified from Turner et al 2012 & Wang

#dN/dt = r*N - a*N*P
#dI/dt = a*N*P - I(t-l)
#dP/dt = b*I(t-l) - a*N*P
            
#Note: this is just a general case of the Lotka-Volterra above
#Since that can be arrived at by setting l = 0
    
def dn_dt2(N, P, r, a):
    return(max(0, N+(r*N - a*N*P)))
    
def di_dt2(N, P, I, I_t_l, a):
    return(max(0, I+(a*N*P - I_t_l)))
    
def dp_dt2(N, P, I_t_l, a, b):
    return(max(0, P+(b*I_t_l - a*N*P)))

def sim_pop2(N_init, I_init, P_init, r, a, b, l, nsteps):
    #Simulates the community for given parameters & number of steps
    track = np.zeros((nsteps, 3))
    N = N_init
    I = I_init
    P = P_init
    for i in range(0, nsteps):
        track[i, 0] = N
        track[i, 1] = I
        track[i, 2] = P
        if i >= l:
            I_t_l = track[int(i-l), 1]
        else:
            I_t_l = 0
        N = dn_dt2(N, P, r, a)
        I = di_dt2(track[i, 0], P, I, I_t_l, a)
        P = dp_dt2(track[i, 0], P, I_t_l, a, b)
    return(track)

#Let's say that one timestep is 1 minute
#Then we can subdivide from there to get nsteps

#That way burst size & latent time can be real numbers

run_time = 12*60 #in minutes
step_time = 0.1 #in minutes

nsteps = int(run_time/step_time)

init_bact_phage_ratio = 100000
moi = 1/init_bact_phage_ratio

doub_time = 20 #minutes
r = math.log(2)/doub_time

a = 1/10000000000 #adsorption parameter
l = 100 #latent time in minutes
b = 100 #burst size

big_track = []
p_init_vals = np.linspace(0, 1000, num = 00)
peak_track = np.zeros((len(p_init_vals), 3))
for j in range(len(p_init_vals)):
    track = sim_pop2(N_init = p_init_vals[j]*1/moi,
                     I_init = 0,
                     P_init = p_init_vals[j],
                     r = r*step_time,
                     a = a*step_time,
                     b = b,
                     l = l/step_time,
                     nsteps = nsteps)
     
    big_track.append(track)
    
    if max(track[:, 0]) == p_init_vals[j]*1/moi:
        peak_track[j, 0] = 0 #the bact pop never grew
    else:
        peak_track[j, 0] = max(track[:, 0])

#plt.subplot(3, 3, cntr)
plt.plot(p_init_vals, peak_track[:, 0])
plt.title("r="+str(r)+", a="+str(a)+", b="+str(b))
plt.show()

j = 10
#Make plot showing this population w/ 2 y axes but shared x axis
#Define plot & axis 1 (for bacteria)
fig, ax1 = plt.subplots()
#Plot bacteria
ax1.plot(range(0, nsteps), big_track[j][:, 0], label = "N", color = 'blue')
#Set up bacterial y axis 
ax1.set_ylabel('N', color='b')
ax1.set_yscale("log")
#Make the phage x axis a twin of the bacterial x axis
ax2 = ax1.twinx()
#Plot phage
ax2.plot(range(0, nsteps), big_track[j][:, 2], label = "P", color = 'red')
#Set up phage y axis
ax2.set_ylabel('P', color='r')
ax2.set_yscale("log")
plt.show()

#if cntr == 9:
#    plt.show()
#    cntr = 1
#else:
#    cntr += 1