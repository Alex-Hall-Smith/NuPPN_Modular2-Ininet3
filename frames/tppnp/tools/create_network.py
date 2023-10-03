Zwidth = 12 # allowed range in Z for given A
maxA = 140

# read full network from isotopedatabase_all.txt
d = np.loadtxt("isotopedatabase_all.txt",usecols=[0,1],skiprows=3)
a = d[:,1].copy(); z = d[:,0].copy()
n = a-z

@np.vectorize
def zopt(a):
    """return Z for which binding energy is minimized for the given A, based on the SEMF"""
    a_s = 17.8; a_v = 15.75; a_c = 0.711; a_a = 23.7
    return int(np.round(2*a_a*(a_c/a**(1./3.) + 4*a_a/a)**(-1)))


# trim network based on width around valley of stability and max. mass number
newn = []; newz = []
for i in range(len(a)):
    zo = zopt(a[i])
    if z[i] < zo + Zwidth/2 and z[i] > zo - Zwidth/2 and a[i] <= maxA:
        newz.append(z[i])
        newn.append(a[i] - z[i])

# plot
clf()
abest = np.arange(1,301)
zbest = zopt(abest)
plot(a-z,z,'.',ms=3,label='full network: %s' % str(len(a)))
plot(newn,newz,'.',ms=3,label='trim network: %s' % str(len(newn)))
plot(abest-zbest,zbest,label='min BE',c='k')
legend(loc='best',fontsize=6)
ylim(-5,90)
ylabel("proton number")
xlabel("neutron number")
minorticks_on()

# write new network file
newn = np.array(newn,np.int32); newz = np.array(newz,np.int32)
newisos = [(newn[i]+newz[i],newz[i]) for i in range(len(newn))]
header = open("isotopedatabase_all.txt","r").readlines()[:3]
body = open("isotopedatabase_all.txt","r").readlines()[3:]

f = open("output_network_%s.txt" % str(len(newn)),"w")
for line in header: f.write(line)
for line in body:
    z = int(line.split()[0])
    a = int(line.split()[1])
    if (a,z) in newisos or "nn" in line:  # don't forget the neutrons!!
        f.write(line)
f.close()

print("Done. Don't forget to change NNN_max in your network file")
