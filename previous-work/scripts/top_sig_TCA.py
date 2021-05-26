import pandas 
import glob
import numpy as np

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
				final.append([cpg, mean])
			counter += 1
	celltypes.append(final)


final = []
for celltype in celltypes:
	mylist = celltype
	mylist.sort(key=lambda x: x[1], reverse=True)
	finalcelltype = mylist[:100]
	final.append(finalcelltype)

finalframes = []
for entry in final:
	newdf = pandas.DataFrame(entry)
	newdf = newdf.set_index(0)
	finalframes.append(newdf)

counter = 1
for entry in finalframes:
	entry.to_csv("celltype"+str(counter)+'top100significant_TCA.csv')
	cpgs = entry.index get
	with open("celltype"+str(counter)+'top100significant_TCA_aslist.csv','w') as out:
		out.write(','.join(cpgs))
	out.close()
	counter += 1
