'''
Title: ASL Grid Search
Desc: Shrinking grid search for best fit seismic source location using a synthetic
    source and amplitude ratios with noise.
Author: Bradford Mack
Date: 3 Dec 19
Last modified: 9 Dec 19
'''

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas
        

#%% Functions

def nCr(n,r):
    return math.factorial(n)/(math.factorial(r)*math.factorial(n-r))


def dist(x1, y1, x2, y2):
    return np.sqrt((x1-x2)**2+(y1-y2)**2)


#%% Input variables

invars = pandas.read_csv('./GridSearchInputs.csv')
gspace = invars['Grid Space (m)'][0]
glen = invars['Grid Length (m)'][0]
shrink = invars['Shrink'][0]
f = invars['Frequency (Hz)'][0]
Q = invars['Q'][0]
v = invars['Velocity (m/s)'][0]
gspace_min = invars['Min Grid Space (m)'][0]
amp_noise = invars['Noise SD'][0]
glen0 = glen.copy()
gspace0 = gspace.copy()

source = pandas.read_csv('./SourceLocs.csv')
sox = source['Source X (m)']
soy = source['Source Y (m)']

#sox = 500 #synthetic source x loc
#soy = 505 #synthetic source y loc
stax = np.array([300,200,550,750,800]) #station x locs
stay = np.array([100,800,800,700,200]) #station y locs

x = (-math.pi*f)/(v*Q) #constant for amplitude ratio formula, do not change

plt.close('all')



#%% Synthetic amplitude ratios at each station

r = np.empty([len(sox),len(stax)]) #radii for each station from synthetic source
amprat = np.empty([len(sox),int(nCr(len(stax),2))]) #amplitude ratios between stations with synthetic source
ampstr = amprat.copy() 
ampstr = ampstr.astype(str) #string array of station pairs
count = 0 #counter for amplitude ratio pairs
noise = np.random.normal(0,amp_noise,len(amprat)) #noise to be added (multiplied by amprat) to amplitude ratios

for t in range(0,len(sox),1):
    for i in range(0,len(stax),1): #calculate radii from synthetic source
        r[t][i] = dist(stax[i], stay[i], sox[t], soy[t])

for t in range(0,len(sox),1): 
    for i in range(0,len(stax)-1,1): #calculate amplitude ratios with synthetic source
        for j in range(i+1,len(stax),1):
            amprat[t][count] = np.sqrt(r[t][j]/r[t][i])*np.exp(x*(r[t][i]-r[t][j]))+noise[count]
            ampstr[:,count] = str(i)+ "/" + str(j)
            count +=1
    count = 0
        

#%% Grid Search
        
ga = np.empty([int(nCr(len(stax),2)),int(np.ceil(glen/gspace)+1),int(np.ceil(glen/gspace)+1)]) #amplitude ratios for source at each grid loc
gx = np.empty(len(ga[0])) #x values for grid
E = ga.copy() #error array for finding best fit loc
rads = np.empty([len(stax),int(np.ceil(glen/gspace)+1),int(np.ceil(glen/gspace)+1)]) #3D radii array for each station vs grid locs
count = 0 #counter for amplitude ratio pairs
e = 0 #summed error for each station pair at each grid loc
re_best = 9999 #finding best RMS error at each grid loc
runs = 1 #counter for number of runs
E_dist = np.empty(len(sox)) #array of distance between best fit and synthetic source in m

for i in range(0,len(gx),1): #set grid x locs for 1st grid
    gx[i] = i*gspace
    
gy = gx.copy() #y values for 1st grid
x_best = np.empty(len(sox))
y_best = x_best.copy()

for t in range (0,len(sox),1): #t is time window for moving source
    while gspace > gspace_min:
        for i in range(0,len(gx),1): #calculate radii from each station to grid locs
            for j in range(0,len(gy),1):
                for k in range(0,len(stax),1):
                    rads[k][j][i] = dist(stax[k],stay[k],gx[i],gy[j]) #station, y, x
                    if rads[k][j][i] == 0:
                        rads[k][j][i] = 0.000001
                
        for i in range(0,len(gx),1): #calculate amplitude ratios and error at each grid loc
            for j in range(0,len(gy),1):
                for k in range(0,len(stax)-1,1):
                    for l in range(k+1,len(stax),1):
                        #print(count)
                        ga[count][j][i] = np.sqrt(rads[l][j][i]/rads[k][j][i])*np.exp(x*(rads[k][j][i]-rads[l][j][i])) #pair, y, x
                        E[count][j][i] = ga[count][j][i]-amprat[t][count] #error at each grid point for each pair
                        count +=1
                count = 0
          
        for i in range(0,len(E[0]),1): #find best fit source loc
            for j in range(0,len(E[0]),1):
                for k in range(0,len(amprat[0])):
                    e += E[k][j][i]**2 #sum of errors
                re = np.sqrt(e) #RMS error
                #print(re)
                if (re<re_best):
                    re_best = re #best RMS error
                    i_best = i #best fit x index
                    j_best = j #best fit y index
                    x_best[t] = gx[i_best] #best fit x value
                    y_best[t] = gy[j_best] #best fit y value
                e = 0 #reset error sum for next set of pairs at grid loc
        print(str(runs) + ' m\370\370se, ' + 't = ' + str(t) + ', Grid Space = ' + str(gspace)[0:4] + ' m')
        gspace *= shrink #shrink grid spacing
        glen *= shrink #shrink grid length
        runs += 1
        :q

        for i in range(0,len(gx),1): #set grid x,y locs
            gx[i] = i*gspace+(x_best[t]-glen/2+gspace/2)
            gy[i] = i*gspace+(y_best[t]-glen/2+gspace/2)
            
    stations, = plt.plot(stax,stay,'bo',label='Stations')
    synsource, = plt.plot(sox[t],soy[t],'ro',label='Synthetic Source')
    plt.axis([0,glen0,0,glen0])
    best_fit, = plt.plot(x_best[t],y_best[t],'k*',label='Best Fit Source Loc')
    plt.legend([stations,synsource,best_fit], ["Stations","Synthetic Source","Best Fit Source Loc"])
    E_dist[t] = dist(sox[t],soy[t],x_best[t],y_best[t])
    
    e = 0
    re_best = 9999
    runs = 1
    gspace = gspace0.copy()
    glen = glen0.copy()

plt.show()
print('We apologise for the fault in the subtitles. Those responsible have been sacked.')