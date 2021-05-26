# using test data

python glint.py --datafile tutorial_files/datafile.txt --covarfile tutorial_files/covariates.txt --phenofile tutorial_files/phenotypes.txt --gsave

# look for outliers in the data
python glint.py --datafile datafile.glint --plot --plotpcs --numpcs 2 --out pcs_plot


# remove outliers
python glint.py --datafile datafile.glint --maxpcstd 1 4 --gsave --out data_cleaned

# capture cell type composition
python glint.py --datafile data_cleaned.glint --refactor --k 6 --covar age gender chip1 chip2 chip3 chip4 chip5 chip6 chip7 chip8 --gsave --out data_cleaned_v2

# infer population structure
python glint.py --datafile data_cleaned_v2.glint --epi --covar rc1 rc2 rc3 rc4 rc5 rc6 --gsave --out data_final




## USING REAL DATA ## 
# had to remove "" marks from dataset 
''' Had to change the code - to change position to str in order to account for
	nan values '''
python glint.py --datafile real_data/methdata.txt --covarfile real_data/covariates.txt --phenofile real_data/ra_status.txt --gsave

# remove outliers
python glint.py --datafile methdata.glint --maxpcstd 1 4 --gsave --out data_cleaned

# capture cell type composition as houseman estimates
python glint.py --datafile data_cleaned.glint --refactor --k 6 --covar age gender smoking_occasional smoking_ex smoking_current --houseman --gsave --out data_cleaned_v2



# and now to format the output to make it easier to load into R... 

matrix = []
colnames = []
rownames = []

counter = 1

with open('data_cleaned_v2.houseman_estimates.txt', 'r') as output:
	for line in output:
			splitline = line.strip('\n').split(' ')
			if counter == 1:
				colnames = splitline
			else:
				matrix_nums = splitline[1:]
				to_add = []
				for entry in matrix_nums:
					to_add.append(float(entry))
				matrix.append(to_add)
				rownames.append(splitline[0])
		counter += 1

colnames = colnames[1:]

df = pandas.DataFrame(matrix, columns=colnames, index=rownames)
df.to_csv('cell_proportions_houseman.csv')
