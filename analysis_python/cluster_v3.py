#
# TIME EVOLUTION OF AGGREGATION
# Use as input files all simulation dcd
# Version as python3 script
# David Malaspina and Jordi Faraudo
#

# Imports
import MDAnalysis as md
import numpy as np
import scipy as sp
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster
import matplotlib.pyplot as plt

#List of dcd files
dcd_files = ['../../equilibration/MDequil.dcd','../MD.dcd','../MDs2.dcd','../MDs3.dcd','../MDs4.dcd','../MDs5.dcd','../MDs6.dcd','../MDs7.dcd','../MDs8.dcd','../MDs9.dcd','../MDs10.dcd','../../MDr3cont/MDs10.dcd','../../MDr3cont/MDs11.dcd','../../MDr3cont/MDs12.dcd','../../MDr3cont/MDs13.dcd','../../MDr3cont/MDs14.dcd','../../MDr3cont/MDs15.dcd','../../MDr3cont/MDs16.dcd','../../MDr3cont/MDs17.dcd','../../MDr3cont/MDs18.dcd','../../MDr3cont/MDs19.dcd','../../MDr3cont/MDs20.dcd'] #Add all dcd files


#Import structure and add dcd files
for dcd_file in dcd_files:
    U = md.Universe('../../input/system.psf', dcd_files)

#create output data files
file1= open('nclust-all.dat', 'w')
file2= open('distri3d-all.dat', 'w') 
fileheader="# time(ps) total_number_of_clusters number_of_clusters_of_sizes_1_to_14 \n"
file2.write(f'{fileheader}')
file3= open('monomers.dat', 'w')
file4= open('meansize-all.dat', 'w')

#
#Do the calculation
#
sel1=U.select_atoms("name P and resname HPO4") #selection

#distri=np.zeros(20)#cluster size distribution
distri3d=np.zeros((15,U.trajectory.n_frames+1))
numdistri=np.zeros((15,U.trajectory.n_frames+1))
dbins=np.arange(0.5,.5)

j=0

for ts in U.trajectory:
    time_ps=round(U.trajectory.time)
    time_ns=U.trajectory.time/1000.0
    frame_number=ts.frame
    #print("Frame: %5d, Time: %8.0f ps" % (ts.frame, U.trajectory.time))
    print("Frame: %5d, Time: %8.2f ns" % (frame_number, time_ns))
    distan=md.lib.distances.self_distance_array(sel1.positions,box=U.dimensions)#distance matrix
    clust=sp.cluster.hierarchy.linkage(distan,'single') #clustering using complete linkage method
    max_d = 6.25 #cutoff based on rdf
    clustdist = fcluster(clust, max_d, criterion='distance')#clasifica los clusters con distancias menores a cutoff
    #fcluster devuelve una matriz donde cada elemento tiene un numero asignado que indica a que cluster pertenece

    nclust=np.amax(clustdist) # amax devuelve el numero mas grande de clustdist que equivale a number of clusters
    #file1.write('{:d}\n'.format(nclust))#print max number of clusters
    file1.write(f'{time_ps} {nclust}\n') #print max number of clusters

        
    #distribucion de clusters
    binss=nclust #-1
    clustsize=np.histogram(clustdist,bins=binss,range=(0.5,nclust+0.5))#la funcion histograma es un poco complicada, devuelve dos arrays
    #el primero con la acumulacion y el segundo con los bins. pero el array de bins siempre es bins+1
    #Por eso usamos el maxclust-1 mas arriba
    #Ej: si pusimos bins=2 devuelve un array (1,2,3), pero solo dos valores de acumulacion
    # pues los bins seran 1-2 y 2-3 y los bins incluiran los numeros acumulados en 3
    
    #np.savetxt(file2, clustsize[0])
    
    for i in range(clustsize[0].size):
        #distri[clustsize[0][i]]=distri[clustsize[0][i]]+(1./clustsize[0].size)

        numdistri[clustsize[0][i],j]=numdistri[clustsize[0][i],j]+1.0
        
        distri3d[clustsize[0][i],j]=distri3d[clustsize[0][i],j]+(1./clustsize[0].size)
    
   
    #screen print of instantaneous distribution only for debugging purposes
    print(numdistri[:,j])    

    #number of monomers
    monomers=numdistri[1,j]
    file3.write(f'{time_ps} {monomers}\n') #save number of monomers

    #average aggregate size
    meanclust=np.mean(clustsize[0])   
    file4.write(f'{time_ps} {meanclust}\n')      
    #file4.write('{:f}\n'.format(meanclust))    
    
    #save full distribution (time, total number of clusters and number of clusters of each type)
    a= numdistri[:,j]
    file2.write(f'{time_ps} {nclust} {a[1]} {a[2]} {a[3]} {a[4]} {a[5]} {a[6]} {a[7]} {a[8]} {a[9]} {a[10]} {a[11]} {a[12]} {a[13]} {a[14]}  \n') #save distribution
    
    #update frame number
    j=j+1
    
#distri=distri/U.trajectory.frame
#distri3d=distri3d/U.trajectory.frame
#file2 = "./distri3d-all.dat"

#np.savetxt(file3, np.transpose([dbins, distri]))
#np.savetxt(file2,distri3d)

#plot results
#plt.plot(distri, label='sys3-all')
#plt.legend(frameon=False, loc='upper right')
#plt.xlim(0,15)
#plt.xlabel('Average Aggregate size')
#plt.ylabel('Probability density')
#plt.savefig('./meansize-all.png', dpi=600)
#fig = plt.figure(figsize=(10, 8), label='sys3-all') 
#plt.contourf(distri3d)
#plt.ylim(0,13)
#plt.xlabel('Frames')
#plt.ylabel('Distribution of HPO4 in clusters')
#plt.savefig('./distri3d-all.png')

