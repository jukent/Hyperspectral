# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 11:31:52 2018

@author: julia
"""
rho_i = 0.917*10**6 #g/m^3
Z = [10,9,8,7]
DZ = (Z[0]-Z[-1])*10**3

#Values from Twomey and Cocks 1989
tau = [1,1.5,2,3,4,6,8,12,16,24,32,48,64,96,128,192,20,28,36,40,44,52,56 \
       ,60,68,92,100,104,108,112,116,120,124,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188]
reff = [5,6,8,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45, \
        47.5,50,52.5,55,57.5,60] #um (LRT supported range 5-60)

##Write the Ice Files
files = []    
for t in tau:
    for r in reff:
        IWC = 2/3*t*(r*10**-6)*rho_i/DZ
        
        namestr = 'IC'+str(IWC)+'_r'+str(r)+'_t'+str(t)+'.DAT'
        files.append(namestr)
        f = open(namestr,"w+")
        line1 = str("\t".join(['#Z (km)','IWC (g/m^3)','r_eff (um)']))
        line2 = str("\t".join([str(Z[0]),'0','0']))
        f.write('{}\n{}\n'.format(line1,line2))
        for d in range(3):
            line = str("\t".join([str(Z[d+1]),str(IWC),str(r)]))
            f.write('{}\n'.format(line))
        f.close()

##Call up the Ice Files in UVSPEC input file
habits=['ghm']#['ghm','solid-column','rough-aggregate']
inputs = []
for file in files:
    for habit in habits:
        namestr='Clouds1_'+file[0:-4]+'_'+habit+'.INP'
        inputs.append(namestr)
        f = open(namestr,"w+")
        
        line1 = 'include ../juke6049/Inputs/UVSPEC_Clouds1.INP'
#        line1 = 'include ../juke6049/Inputs/CERES_23h29m.INP'
        f.write('{}\n'.format(line1))
        line2 = 'ic_properties baum_v36 interpolate'
        line3 = 'ic_habit '+str(habit)
        line4 = 'ic_file 1D Data/'+str(file)
        f.write('\n{}\n{}\n{}\n'.format(line2,line3,line4))
        f.close()
        
#Create Bash scripting file
outputs = []
f = open('IceCloud.sh',"w+")
shebang = "#!/bin/bash"
f.write('{}\n'.format(shebang))
for i in inputs:
    line = str("\t".join(['../bin/uvspec','< Inputs/'+str(i)+' >','Outputs/'+str(i[0:-4])+'_baum.dat']))
    f.write('{}\n'.format(line))
    outputs.append(str(i[0:-4])+'.dat')
f.close()