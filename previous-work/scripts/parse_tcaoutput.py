import pandas
import glob
import numpy as np

rankedlist = 'data_cleaned_v2.refactor.rankedlist.txt'
ranked = []

with open(rankedlist, 'r') as r:
	for line in r:
		cpg = line.strip('\n')
		ranked.append(cpg)

# create a list of the top 100 discriminatory cpgs
#top = ranked[:100]
top = ranked[:50]

filelist = glob.glob('ct[1-6].csv')


celltypes = []


for file in filelist:
	final = []
	with open(file, 'r') as matrix:
		counter = 0
		for line in matrix:
			if counter != 0:
				line_split = line.strip('\n').split(',')
				data = line_split[1:]
				cpg = line_split[0].strip('"')
				data = [float(i) for i in data]
				mean = np.mean(data)
				if cpg in top:
					final.append([cpg, mean])
			counter += 1
	celltypes.append(final)

finalframes = []
for entry in celltypes:
	newdf = pandas.DataFrame(entry)
	newdf = newdf.set_index(0)
	finalframes.append(newdf)

counter = 1
for entry in finalframes:
	entry.to_csv("tcacelltype"+str(counter)+'top100.csv')
	counter += 1



cellt = []
cpgs = []
for entry in celltypes:
	data = []
	for cpg in entry:
		data.append(cpg[1])
		if cpg[0] not in cpgs:
			cpgs.append(cpg[0])
	cellt.append(data)	

# create joint matrix
joint_df = pandas.DataFrame(cellt)
joint_df = joint_df.transpose()
joint_df['cpgs'] = cpgs
joint_df = joint_df.set_index('cpgs')
joint_df.columns=['celltype1', 'celltype2', 'celltype3', 'celltype4', 'celltype5', 'celltype6']
# export to csv
joint_df.to_csv("tca_allcelltypes_top100discrim.csv")