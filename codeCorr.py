from prody import *
from pylab import *
ion()
pdbids = ['6vyb']

open_sequence = '''MGILPSPGMPALLSLVSLLSVLLMGCVAETGTQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNEVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPSGAGSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGA'''

blast_record = blastPDB(open_sequence)

pdbids = blast_record.getHits()

pdbfiles = fetchPDB(*pdbids, compressed=False)



len(pdbids)

ref_structure = parsePDB('/home/daya/Documentos/NORMALCOV/pca/no2/6vyb.pdb')

ref_selection = ref_structure.select('resnum 27 to 1127 and calpha')
ref_chain = ref_selection.getHierView().getChain('B')
#ref_chainB = ref_selection.getHierView().getChain('B')
#ref_chainC = ref_selection.getHierView().getChain('C')
repr(ref_chain)

startLogfile('open_pca')

ensemble = PDBEnsemble('open x-ray')
ensemble.setAtoms(ref_chain)
ensemble.setCoords(ref_chain)

a=[]

for pdbfile in pdbids:
	structure = parsePDB(pdbfile, subset='calpha')
	mappings = mapOntoChain(structure, ref_chain)
	a.append(mappings)

for i in range(len(a)):
	if len(a[i]) !=0:
		atommap = a[i][0][0]
		ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))


ensemble.iterpose()
closeLogfile('open_pca')

writePDB('open_spike_ensemblei.pdb', ensemble)


pca = PCA('SPIKE X-Ray')

pca.buildCovariance(ensemble)

pca.calcModes('all')

pca_svd = PCA('spike svd')

pca_svd.performSVD(ensemble)

print(abs(pca_svd.getEigvals()[:1] - pca.getEigvals()).max())

print(abs(calcOverlap(pca, pca_svd).diagonal()[:1]).min())

anm = ANM('spike')

anm.buildHessian(ref_chain)

anm.calcModes()

saveModel(pca)

saveModel(anm)

saveEnsemble(ensemble)

writePDB('p38_ref_chain.pdb', ref_chain)

####################################################

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




