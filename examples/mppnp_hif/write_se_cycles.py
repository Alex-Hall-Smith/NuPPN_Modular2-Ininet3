
##############Creates se files:
##This script writes the data of cycle cycle_all from file run_H5 into
##the sefile outputfile for a number of cycles given by cycles variable

import nugridse as mp
import sewrite as sw
import numpy as np


###############Read se file to get header information


run_H5='e2D14.0077501.se.h5'
cycle_all=77991

sefiles=mp.se('.',run_H5,rewrite=True)

hattrs=sefiles.se.hattrs
hattrs_data=[]
for k in range(len(hattrs)):
	if hattrs[k]=='age_unit':
		hattrs_data.append(1.)
		continue
	hattrs_data.append(sefiles.get(hattrs[k]))

cattrs=[]
for k in range(len(sefiles.se.cattrs)):
	if ((sefiles.se.cattrs[k] == 'age') or (sefiles.se.cattrs[k] == 'deltat')) or (sefiles.se.cattrs[k]=='shellnb'):
		continue
	cattrs.append(str(sefiles.se.cattrs[k]))
cattrs_data=[]
for k in range(len(cattrs)):
	cattrs_data.append(sefiles.get(cycle_all,str(cattrs[k])))

dcols=["dcoeff","radius","rho","temperature"]

dcols_data=[]
for k in range(len(dcols)):
	#print dcols[k]
        dcols_data.append(sefiles.get(cycle_all,dcols[k]))

mass2=sefiles.get(cycle_all,'mass')

age=sefiles.get(cycle_all,'age') # unit in years
deltat=sefiles.get(cycle_all,'deltat') #unit in s


#use the radius below
#plt.figure()
#plt.plot(mass,radius,label='radius special')
#plt.legend()
#plt.figure()
########Grid changes/mesh refinement:
#mass2 represents the se restart grid, see first section in file
#xm in fortran with maximum poosibble grid size entries/empty

#dxm from ppn_frame.input
dxm=1e-4
xmrmaxi=0.5987
xmrmin=0.57556

xmsg=[0.5880, 0.5890 ,0.5979]
dxmsg=[4.e-4,6.e-5,4.e-4]

xm = [xmrmin]
dq=[7.e-4] #mass grid se, set first cell size to this value
dq=[6.e-4]
dq=[7.5e-4]
#dq between 0 and 2*4e-4, for dq=7e-4 shows no negative dq anymore
#dq between 0 and 2*4e-4

xmrmax=max(mass2)
xmrmax = min(xmrmaxi,xmrmax)

ddxm=dxm
ddxm1=dxm
i=0


while(True):
        ddxm=dxm
        for j in range(len(xmsg))[::-1]:
                if xm[i]+dxmsg[j] < xmsg[j]:
                        ddxm=dxmsg[j]
        #print 'change xm at mass',xm[i],ddxm,i
        i+=1
        xm.append(round(xm[i-1] + ddxm,10))
	dq.append(round(2* (ddxm - dq[-1]/2.),10 ))
        if xm[i] >=xmrmax:
                xm[i] = xmrmax
                break
#result: xm,xmm, xm center, xmm, outer boundary
xm1=xm
shellnb=len(xm)

#Checck ingestion arad: how many cells?
#if(xm(j).ge.xmrmax-4.d-4)then
counter=0
for k in range(len(xm)):
	if xm[k]>=(xm[-1]-4e-4 ):
		counter+=1

print 'Number of cells in which H will be ingested:'
print counter

##below for findig the right start dq which gives the for all other dq>0

#plt.figure()
#plt.plot(xm1,dq,label='dq/delta mass')

#diff=(np.array(xm[1:])-np.array(xm[:-1]))
#plt.plot(xm[1:],diff,label='xm diff',marker='o',markevery=1)


############ Interpolate other parameter (rho,T) on new grid
#adapt those parameter below on new grid/except mass, 
#mass2: outer boundary of cell, from e2D14.0077501.se.h5

#the original routine l. 1862 needs the grid of the se input file xmm
#+ the new mppnp grid xm (with center mass coord)
# therefore take mass2 as xmm, and xm from above

dcols_name=["dcoeff","radius","rho","temperature"]
#data in dcols_data

#number of shells mi, 

#mi=len(mass2)
xmm=mass2 #from se file

from operator import itemgetter

xfind = 0.
ifindpos = 0
rhoppg=[]
t9ppg=[]
dppg = []
rppg = []

#carries the modified diffusion coeff
dppgLagr=[]
dppgLagr_split=[]
dppg_split=[]

#print hattrs
dcoeff_unit=hattrs_data[hattrs.index('dcoeff_unit')][0]
rho_unit=hattrs_data[hattrs.index('rho_unit')][0]
temperature_unit=hattrs_data[hattrs.index('temperature_unit')][0]
radius_unit=hattrs_data[hattrs.index('radius_unit')][0]

#print test
#print 'age unit',hattrs_data[hattrs.index('age_unit')][0] #in years


#print dcoeff_unit,rho_unit,temperature_unit,radius_unit

for j in range(len(xm)):
	xfind = xm[j]
        #ifindpos = MINLOC (xmm(1:mi), MASK = xmm(1:mi) .GT. xfind)
	
	ifindpos=list(xmm).index(min(x for x in xmm if x > xfind))	
	
	x=xfind
	x1=xmm[ifindpos]
	x2=xmm[ifindpos+1]
	y1=dcols_data[2][ifindpos] *rho_unit #dcol_data[2] is rho
	y2=dcols_data[2][ifindpos+1] *rho_unit
	rhoppg.append(y1+((y2-y1)/(x2-x1))*(x-x1) )	
        # if (ifindpos(1).eq.0) print *,"j, xfind, ifindpos: ",j, xfind,
     #$        ifindpos,(xmm(i),i=1,2)
         #rhoppg(j)= ylin(xfind, xmm(ifindpos(1)) ,xmm(ifindpos(1)+1)
     #$        ,rhose(ifindpos(1)),rhose(ifindpos(1)+1))

	y1=dcols_data[3][ifindpos] *temperature_unit #dcol_data[3] temperature
	y2=dcols_data[3][ifindpos+1] *temperature_unit
	t9ppg.append(  ( y1+((y2-y1)/(x2-x1))*(x-x1)))
	#print 't9ppg'
	#print ( y1+((y2-y1)/(x2-x1))*(x-x1))
     #    t9ppg(j)= ylin(xfind, xmm(ifindpos(1)) ,xmm(ifindpos(1)+1)
     #$        ,t8se(ifindpos(1)),t8se(ifindpos(1)+1))


		#	
	y1=dcols_data[0][ifindpos]*dcoeff_unit  # [d] = cm^2/s
	y2=dcols_data[0][ifindpos+1]*dcoeff_unit
	dppg.append(  ( y1+((y2-y1)/(x2-x1))*(x-x1)) )

     #    dppg(j)= ylin(xfind, xmm(ifindpos(1)) ,xmm(ifindpos(1)+1)
     #$        ,dse(ifindpos(1)),dse(ifindpos(1)+1))

	
	y1=dcols_data[1][ifindpos]  *radius_unit#radius
	y2=dcols_data[1][ifindpos+1] *radius_unit
	rppg.append(  ( y1+((y2-y1)/(x2-x1))*(x-x1)))
     #    rppg(j)= ylin(xfind, xmm(ifindpos(1)) ,xmm(ifindpos(1)+1)
     #$        ,rse(ifindpos(1)),rse(ifindpos(1)+1))
	
	###here in loop adapt diffusion coefficient	
		
        dconst = 3.0589703419947736e-11
	dconst=6.314572864321608e-33	
#C *** transform to Lagrangian coordinate: D_Lagr = (4*PI*R**2*RHO)**2*D_Eularian
#C *** and also delta_xm in mixdiffnet is in solar mass units, and since there is 
#C *** a second order spatial derivative:
#C     D_Lagr,mppn = (4*PI*R^2*RHO/Msun)**2*D_Eularian
#C solar constants: radius 6.96E10cm, mass 1.99E33 g
#C units from reading HPF files: [R]=Rsun, [rho]=cgs         
#C     4*Pi*Rsun^2/Msun = 3.0589703419947736e-11



	dppgLagr.append(0)
	dppgLagr_split.append(0)

        if (dppg[j]>0):
        	if (j==1):
        		RH = np.log(rhoppg[j]) #rhoppg 
        		HLOG = np.log(rppg[j]) #rppg in cm
        	else:
               		RH = (np.log(rhoppg[j]) + np.log(rhoppg[j-1]))/2.
               		HLOG = (np.log(rppg[j]) + np.log(rppg[j-1]))/2.
            	#endif   ! dppgLagr: Lagrangian coordinates assuming solar mass units
            #if (ksubc.ge.isplitstart .and. isplit.eq.1) then
		a1=1.e4           
		a2=1.e7
		xmsplit=0.5885
		#print xm, xmsplit
		dsplit = dppg[j]  / (1. + a2*np.exp(-(a1*(xm[j]-xmsplit))**2))
		    #else
		    #   dsplit = dppg(J)
		    #end if

		dppg_split.append(dppg[j]  / (1. + a2*np.exp(-(a1*(xm[j]-xmsplit))**2)))	

		dppgLagr_split[j] = dsplit * (dconst**2)*np.exp(2*(2*HLOG+RH)) 
		dsplit = dppg[j]
		dppgLagr[j] = dsplit * (dconst**2)*np.exp(2*(2*HLOG+RH))
		#ppgLagr(j) = dsplit * dconst**2*exp(2*(2*HLOG+RH))


#dppgLagr(j) = dppg(J) * dconst**2.*exp(2.d0*(2.d0*HLOG+RH))


##########Do some test polotting
if (False):
	'''
	plt.plot(xm1,t9ppg,label='T9',marker='o',markevery=1)
	plt.xlabel('mass coordinate')
	plt.ylabel('temperature')
	plt.figure()
	plt.legend()
	plt.plot(xm1,rhoppg,label='Rho',marker='o',markevery=1)
	plt.legend()
	plt.xlabel('mass coordinate')
	plt.ylabel('Rho')
	plt.figure()
	plt.plot(xm1,np.log10(dppgLagr),label='se diffusion')
	plt.plot(xm1,np.log10(dppgLagr_split),label='Modification of se diffusion')
	#plt.yscale('log')
        plt.xlabel('mass coordinate')
        plt.ylabel('Log Diffusion coefficient')
	plt.legend()
	plt.figure()
	plt.plot(xm1,dppg,label='dppg')
	plt.plot(xm1,dppg_split,label='dppg split')
	plt.legend()
	plt.figure()
	plt.plot(xm1,dppg,label='dppg')
	plt.plot(xm1,dppg_split,label='dppg split')
	plt.legend()
        plt.xlabel('mass coordinate')
        plt.ylabel('D [cm^2/s]')
	'''
	plt.figure('test')
	plt.plot(xm1,dq,label='dq '+str(dq[0]),marker='o',markevery=1)
	#diff=(np.array(xm1[1:])-np.array(xm1[:-1]))
	#plt.plot(xm1[1:],diff,label='xm1 diff',marker='o',markevery=1)
	#plt.legend()
	#plt.ylim(-1e-5,1e-3)
#print 'XM ...'
#print xm1[:5]
###############Write out cycle_all for each cycle of cycles
##Important: change 'age','deltat' so that timesteps of 63s

deltat=63.  #s

cycles=range(77991,(77991+3002))

cycle_bndy=[77991,78001,78501,79001,79501,80001,80501,81001]

outputfiles=['e2D14_hif.0077501.se.h5','e2D14_hif.0078001.se.h5','e2D14_hif.0078501.se.h5','e2D14_hif.0079001.se.h5','e2D14_hif.0079501.se.h5','e2D14_hif.0080001.se.h5','e2D14_hif.0080501.se.h5']
i=0



#for se files inverted mass, surface is first element
dppg=dppg[::-1]
dppgsplit=dppg_split[::-1]
dppgLagr=dppgLagr[::-1]
rppg=rppg[::-1]
xm1=xm1[::-1]
dq=dq[::-1]
rhoppg=rhoppg[::-1]

t9ppg=t9ppg[::-1]
print 'compare types'
print type(dppgLagr),type(dppg),type(dppg_split)
for cyc in cycles:
	cyc = int(cyc)
	if cyc==cycle_bndy[i]:
		print 'start writing in file ',outputfiles[i]
		f=sw.startfile(outputfiles[i])
		f.write_hattr(hattr_name=hattrs,hattr_val=hattrs_data)
		#Add specific header info about first cycle
		f.write_hattr(hattr_name=['firstcycle'],hattr_val=[np.array([cyc], dtype=np.int32)])
		i+=1

	#f.write_dcol(cyc,dcol_name=dcols,dcol_val=dcols_data)	

	#start with split
	#here the Mesa se mass coord for 'mass' and delta mass were chosen
	if cyc>=(951+77991):
		f.write_dcol(cyc,dcol_name=["dcoeff","radius","mass","delta_mass","rho","temperature"],dcol_val=[dppgsplit,rppg,xm1,dq,rhoppg,t9ppg]) 
	else:
		f.write_dcol(cyc,dcol_name=["dcoeff","radius","mass","delta_mass","rho","temperature"],dcol_val=[dppg,rppg,xm1,dq,rhoppg,t9ppg])

	f.write_cattr(cyc,cattr_name=cattrs,cattr_val=cattrs_data)
	#due to timesteps of 63s
	f.write_cattr(cyc,cattr_name=['deltat','age','shellnb'],cattr_val=[deltat,age,shellnb])	
	age+=(deltat) #age als oin years #*3.1689e-8)
