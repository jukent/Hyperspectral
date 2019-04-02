# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 11:31:52 2018

@author: julia
"""
rho_w = 1*10**6 #g/m^3
Z = [10,9,8,7]
DZ = (Z[0]-Z[-1])*10**3

#Values from Twomey and Cocks, 1989
tau = [1,1.5,2,3,4,6,8,12,16,24,32,48,64,96,128,192]
reff = [5,6,8,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40] #um (LRT supported range 5-60)

#Write the Water Files
files = []    
for t in tau:
    for r in reff:
        LWC = 2/3*t*(r*10**-6)*rho_w/DZ
        
        namestr = 'WC'+str(LWC)+'_r'+str(r)+'_t'+str(t)+'.DAT'
        files.append(namestr)
        f = open(namestr,"w+")
        line1 = str("\t".join(['#Z','LWC','R_eff']))
        line2 = str("\t".join(['#(km)','(g/m^3)','(um)']))
        line3 = str("\t".join([str(Z[0]),'0','0']))
        f.write('{}\n{}\n{}\n'.format(line1,line2,line3))
        for d in range(0,len(Z[0:-1])):
            line = str("\t".join([str(Z[d+1]),str(LWC),str(r)]))
            f.write('{}\n'.format(line))
        f.close()

#Call up the Ice Files in UVSPEC input file
inputs = []
for file in files:
    namestr='CERES_21h30m_water_'+file[0:-4]+'_550nm.INP'
    inputs.append(namestr)
    f = open(namestr,"w+")
    
    line1 = 'include ../juke6049/Inputs/CERES_21h30m.INP'
#    line1 = 'include ../juke6049/Inputs/UVSPEC_Desert_Vegetation_Clouds_550nm.INP'
    f.write('{}\n'.format(line1))
    f.write('\n')
    line2 = 'wc_file 1D Data/'+str(file)
    line3 = 'verbose'
    f.write('{}\n{}\n'.format(line2,line3))
    f.close()
        
#Create Bash scripting file
outputs = []
f = open('CERES_21h30m_Bash.sh',"w+")
shebang = "#!/bin/bash"
f.write('{}\n'.format(shebang))
for i in inputs:
    line = str("\t".join(['(../bin/uvspec','< Inputs/'+str(i)+' >','Outputs/'+str(i[0:-4])+'.dat)','>& Verbose/'+str(i[0:-4])+'_verbose.txt']))
    f.write('{}\n'.format(line))
    outputs.append(str(i[0:-4])+'.dat')
f.close()