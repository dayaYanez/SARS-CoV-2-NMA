from prody import *
from matplotlib.pylab import *
import os

ion()

openS = parsePDB('/home/daya/Documentos/NORMALCOV/Correccion/Close state mutation/N501Y/Elastic Netwoek Models/ANM/closeN.pdb')

calphas = openS.select('protein and name CA')
calphas2 = openS.select('calpha')

if calphas==calphas2:
	anm = ANM('Close State N501Y ANM analysis')
	anm.buildHessian(calphas)
	anm.getHessian()
	anm.calcModes(300, zeros=True)
else: 
	print("No same Atoms selection")

vals = anm.getEigvals().round(4)

fileO = open("EigenValuesANM.txt", "w")
for i in range(len(vals)):
	fileO.write(str(vals[i]) + os.linesep)
fileO.close()

vecs = anm.getEigvecs().round(4)

file1 = open("EigenVectorsANM.txt", "w")
for i in range(len(vecs)):
	file1.write(str(vecs[i]) + os.linesep)
file1.close()

cov = anm.getCovariance().round(3)

file2 = open("CovarianceANM.txt", "w")
for i in range(len(cov)):
	file2.write(str(cov[i]) + os.linesep)
file2.close()

for i in range(len(anm)):
	if vals[i] != 0:
		print(i)
		slowest_mode = anm[i]
		print(slowest_mode.getEigval().round(3))
		print(slowest_mode.getEigvec().round(3))
		writeNMD('opnS_anm_modes.nmd', anm[:15], calphas)
		break
saveModel(anm)
showCrossCorr(anm)
