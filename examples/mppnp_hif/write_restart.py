
####Script to create .restart.h5 file (with one cycle)
   #based on the abundance of the ABUPP###.DAT (mppnp output first cycle)
   #and the header infos from the se file. Do a mass grid refinement that
   #of the se grid to match the ABuPP grid

############Read restart0077991.check (from Falks RUN48 dir) file to get abundances,Z,A,(isomeric_state)

#guess its the same structure as for the ABUPP00779910000.DAT files (inside to surface)

import utils as u

f1=open('restart0077991.check')
lines=f1.readlines()
f1.close()

shells_massfrac=[]
A=[]
Z=[]
mass=[]
for k in range(len(lines)):
	#skip header
	if k<2:
		continue
	if (k==2):
		#column titles
		isotopes=[]
		line=lines[k]#[70:-1]#200]#150]
		#print line
		idx_iso=[]
		isotopes1=[]
		for i in xrange(0,len(line),11):
			isotopes1.append(line[i:i+11].strip())
			#idx_iso.append(line[i:i+11][:5])
		isotopes1=isotopes1[:-1]
		for iso in range(len(isotopes1)):		
			if 'NEUT' in isotopes1[iso]:
				A.append(1)
				isoname=isotopes1[iso][:-3].strip()
				isotopes.append( isoname[0]+isoname[1:].lower()+'-'+str(A[-1]))
			elif 'PROT' in isotopes1[iso]:
				A.append(1)
				isotopes.append( 'H'+'-'+str(A[-1])) 
			else:
				A.append(int(isotopes1[iso][-3:]))
				isoname=isotopes1[iso][:-3].strip()
				isotopes.append( isoname[0]+isoname[1:].lower()+'-'+str(A[-1]))
		#print len(isotopes)
		#print len(isotopes)
		#print len(idx_iso)
		print isotopes
		#print A
		continue
	massfracs=lines[k].split()#[7:]
	#if k==4:
		#print lines[k].split()	
		#print len(massfracs)
		#print massfracs
	shells_massfrac.append(massfracs)

#print 'shells massffrac'
#print shells_massfrac
u.convert_specie_naming_from_h5_to_ppn(isotopes[2:])

Z=[0,1]+[int(k) for k in list(u.znum_int)]
#A1=[1,1]+[int(k) for k in list(u.amass_int)]


#Assume no isomeric states!!!:
isomeric_state=[1]*len(isotopes)


###############Read se file to get header information

import nugridse as mp
import sewrite as sw


run_H5='e2D14.0077501.se.h5'
cycle_all=77991

sefiles=mp.se('.',run_H5,rewrite=True)

hattrs=sefiles.se.hattrs
hattrs=[]
notneeded_header=['HDF5_version','modname','dcoeff_unit','radius_unit','rho_unit','temperature_unit']
for k in  range(len(sefiles.se.hattrs)):
	if sefiles.se.hattrs[k] in notneeded_header:
		continue
	hattrs.append(sefiles.se.hattrs[k])
hattrs_data=[]
for k in range(len(hattrs)):
        hattrs_data.append(sefiles.get(hattrs[k]))

mass2=[]
dcols=[]
for k in range(len(sefiles.se.dcols)):
        if sefiles.se.dcols[k] == 'yps':
                continue
        dcols.append(sefiles.se.dcols[k])
dcols_data=[]
for k in range(len(dcols)):
	if dcols[k]=='mass':
		mass2=sefiles.get(cycle_all,dcols[k])
        dcols_data.append(sefiles.get(cycle_all,dcols[k]))

cattrs=[]
for k in range(len(sefiles.se.cattrs)):
        cattrs.append(str(sefiles.se.cattrs[k]))
cattrs_data=[]
for k in range(len(cattrs)):
	if sefiles.se.cattrs[k]=='shellnb':
		cattrs_data.append(len(shells_massfrac))
		continue
        cattrs_data.append(sefiles.get(cycle_all,str(cattrs[k])))

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



######## Do test plotting


if(False):

	#for se files: index: 0: highest m, last: 0
	#therefore:
	mass2=mass2[::-1]
	plt.plot(range(len(mass)),mass,marker='d',markevery=10,label='ABUPP DAT')
	plt.plot(range(len(mass2)),mass2)
	plt.plot(range(len(xm)),xm,label='xm',marker='o',markevery=10)
	plt.plot([0,1000],[xmrmaxi,xmrmaxi],label='--')
	plt.plot([0,1000],[xmrmin,xmrmin],label='--')

	plt.legend()




######Write restart cycle here#########


name='e2D14_hif.0077501.restart.h5'

f=sw.startfile(name)

#Same header as in se files
f.write_hattr(hattr_name=hattrs+['zisnb'],hattr_val=hattrs_data+[len(shells_massfrac[0])])
#f.write_hattr(hattr_name='zisnb',hattr_val=len(shells_massfrac[0]))#for number of isotopes

#A,Z, isomeric state, mass,iso_massf is all needed in restart file (except header)

f.write_table('A',A)
f.write_table('Z',Z)
f.write_table('isomeric_state',isomeric_state)
print 'last test#######################################################################'
print shells_massfrac[0]
print shells_massfrac[1]


#right units? since im writing out info header_unit + hattrs from se > consistent 
f.write_dcol(cycle_all,dcol_name=['mass']+['iso_massf'],dcol_val=[xm1]+[shells_massfrac])

#write column attr,really needed: "shellnb", but write out more,e.g. age
f.write_cattr(cycle_all,cattr_name=cattrs,cattr_val=cattrs_data)










