import pandas 
import glob

df = pandas.read_csv("RA_HIRE_pvalues.csv")

# add row names
matrix = pandas.read_csv('RAmeth_matrix_linesremoved_10Kmostvar.csv')

celltypes = [df[df.columns[0]], df[df.columns[1]], df[df.columns[2]], df[df.columns[3]],
			 df[df.columns[4]], df[df.columns[5]]]

# convert dfs to csv files
counter = 1
for entry in celltypes:
	entry.to_csv("celltype"+str(counter)+'.csv')
	counter += 1

# get cpg names
cpgs = list(matrix[matrix.columns[0]])

# load csv files back in 

filelist = glob.glob('celltype[1-6].csv')

finaldata = []
for entry in filelist:
	datalines = []
	with open(entry, 'r') as data:
		for line in data:
			split_line = line.strip('\n').split(',')
			datalines.append(float(split_line[1]))
	finaldata.append(datalines)

datalines = []

for celltype in finaldata:
	cell_type = []
	for cpg in celltype:
		cpg_name = cpgs[celltype.index(cpg)]
		print cpg_name, cpg
		cell_type.append([cpg_name, cpg])
	datalines.append(cell_type)


final = []
for celltype in datalines:
	mylist = celltype
	mylist.sort(key=lambda x: x[1])
	finalcelltype = mylist[:100]
	final.append(finalcelltype)

# convert back to dfs and export to csvs
finalframes = []
for entry in final:
	newdf = pandas.DataFrame(entry)
	newdf = newdf.set_index(0)
	finalframes.append(newdf)

counter = 1
for entry in finalframes:
	entry.to_csv("celltype"+str(counter)+'top100significant_HIRE.csv')
	cpgs = entry.index
	with open("celltype"+str(counter)+'top100significant_HIRE_aslist.csv','w') as out:
		out.write(','.join(cpgs))
	out.close()
	counter += 1

