import leafs as l

model = "J1a"
run = "J01a"
stride = 1

# load tracer data
if True:
    t = l.leafs_readtracer(model=model, snappath="./sncode_ONe."+run+"/output/")
    rt = t.loadalltracers()

# extract the following information, taking time-evolutions only from 1 time
# step before the level-set crosses the particle
if True:
    ntimesteps, all_times = t.getTimes()
    # ini/fin positions, mass, time, temperature, density, particle id, level-set
    # number of time steps, number of tracer particles
    posini = []; posfin = []; dm = []; time = []; temp = []; dens = []; ids = []; ls=[];
    nstp = []; ntracers = 0

    for i in range(0, t.npart, stride):
        ntracers += 1
        ids.append(i+1)
        dm.append(t.tmass[i])
        posini.append(np.array([rt.data[0,i,0:2]], np.float64))
        posfin.append(np.array([rt.data[-1,i,0:2]], np.float64))
        burning = np.where(rt.data[:,i,14] > 0.)
        if size(burning) == 0:
            nstp.append(2)
            time.append(np.array([all_times[0],all_times[-1]],np.float64))
            temp.append(np.array([5.e5,5.e5],np.float64))
            dens.append(np.array([rt.data[0,i,3],rt.data[-1,i,3]],np.float64))
            ls.append(np.array([-1e6,-1e6],np.float64))
        else:
            start = max(np.min(burning) - 1,0)
            time.append(all_times[start:])
            temp.append(rt.data[start:,i,4])
            dens.append(rt.data[start:,i,3])
            ls.append(rt.data[start:,i,14])
            nstp.append(len(all_times[start:]))

# write file for tppnp
if True:
    fb = open(run+"_tp.lsb","wb")

    for i in range(ntracers):
        fb.write(np.int32(ids[i]))
        fb.write(np.int32(nstp[i]))
        fb.write(np.float64(posini[i]))
        fb.write(np.float64(posfin[i]))
        fb.write(np.float64(dm[i]))
        fb.write(np.ascontiguousarray(time[i]))
        fb.write(np.ascontiguousarray(temp[i]))
        fb.write(np.ascontiguousarray(dens[i]))

    fb.close()
