from prody import *
from pylab import *
ion()

pca = loadModel('/home/daya/Documentos/NORMALCOV/pca/no2/SPIKE_X-Ray.pca.npz')
anm= loadModel('/home/daya/Documentos/NORMALCOV/pca/no2/spike.anm.npz')
ensemble = loadEnsemble('/home/daya/Documentos/NORMALCOV/pca/no2/open_x-ray.ens.npz')

for mode in pca[:10]:
	var = calcFractVariance(mode)*100
	print(mode, var)

for mode in pca[:3]: # Print PCA mode collectivity
	coll = calcCollectivity(mode)
	print(mode, coll)

pca_coords, anm_coords = calcCrossProjection(ensemble, pca[0], anm[2])

print(np.corrcoef(pca_coords, anm_coords))

for i in range(len(anm)):
	if vals[i] != 0:
		print(i)
		slowest_mode = anm[i]
		print(slowest_mode.getEigval().round(3))
		print(slowest_mode.getEigvec().round(3))
		writeNMD('ensemble.nmd', anm[:15], calphas)
		break


color_list = ['red', 'red', 'red', 'red', 'yellow', 'red', 'yellow', 'yellow', 'red', 'red', 'blue', 'blue', 'red', 'blue' , 'blue', 'red', 'blue', 'red', 'green', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'blue', 'red', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'green', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'blue', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red',  'red', 'red', 'red']


color2label = {'red': 'Glucoside bound', 'blue':'Inhibitor bound', 'yellow':'Peptide/Protein Bound with ACE2', 'green':'Unbound'}

label_list = [color2label[color] for color in color_list]


