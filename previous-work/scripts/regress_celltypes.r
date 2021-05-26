library(devtools)

hireoutput <- data.matrix(read.csv("hire_allcelltypes_top100discrim.csv"))
hireoutput <- hireoutput[,-1]
tcaoutput <- data.matrix(read.csv("tca_allcelltypes_top100discrim.csv"))
tcaoutput <- tcaoutput[,-1]
test_tca <- data.matrix(read.csv('tcacelltype1top100.csv'))
test_tca <- test_tca[,2]
test_hire <- data.matrix(read.csv('celltype1top100.csv'))
test_hire <- test_hire[,2]

# ref analysis
ref <- data.matrix(read.csv('10kref.csv'))

# define celltypes
hire_1 <- data.matrix(read.csv('celltype1top100.csv'))[,2]
hire_2 <- data.matrix(read.csv('celltype2top100.csv'))[,2]
hire3 <- data.matrix(read.csv('celltype3top100.csv'))[,2]
hire4 <- data.matrix(read.csv('celltype4top100.csv'))[,2]
hire5 <- data.matrix(read.csv('celltype5top100.csv'))[,2]
hire6 <- data.matrix(read.csv('celltype6top100.csv'))[,2]
tca1 <- data.matrix(read.csv('tcacelltype1top100.csv'))[,2]
tca2 <- data.matrix(read.csv('tcacelltype2top100.csv'))[,2]
tca3 <- data.matrix(read.csv('tcacelltype3top100.csv'))[,2]
tca4 <- data.matrix(read.csv('tcacelltype4top100.csv'))[,2]
tca5 <- data.matrix(read.csv('tcacelltype5top100.csv'))[,2]
tca6 <- data.matrix(read.csv('tcacelltype6top100.csv'))[,2]

# proof of principal
fit_hire <- lm(test_hire~., data=data.frame(hireoutput))
fit_tca <- lm(test_tca~., data=data.frame(tcaoutput))
fit_ref <- lm(test_hire~., data=data.frame(ref))
summary(fit_hire)
summary(fit_tca)


sink("test.txt")

# matching HIRE celltype 1
fit_hire1 <- lm(hire_1~., data=data.frame(tcaoutput))
summary(fit_hire1)  # cell type 3 - 0.00655 (3*)

# matching HIRE celltype 2
fit_hire2 <- lm(hire_2~., data=data.frame(tcaoutput))
summary(fit_hire2) # cell type 4 - 0.0171(1*)


# matching HIRE celltype 3
fit_hire3 <- lm(hire3~., data=data.frame(tcaoutput))
summary(fit_hire3) # cell type 3 - 8.09e-10 (3*)

# matching HIRE celltype 4
fit_hire4 <- lm(hire4~., data=data.frame(tcaoutput))
summary(fit_hire4) # celltype 3 - 0.00523 (2*)

# matching HIRE celltype 5
fit_hire5 <- lm(hire5~., data=data.frame(tcaoutput))
summary(fit_hire5) # cell type 5 0.0877 (.)

# matching HIRE celltype 6
fit_hire6 <- lm(hire6~., data=data.frame(tcaoutput))
summary(fit_hire6) # cell type 3 0.00032 (3*)


# matching TCA celltype 1
fit_tc1 <- lm(tca1~., data=data.frame(hireoutput))
summary(fit_tc1) # cell type 1 - 0.00173 (2*)

# matching TCA celltype 2
fit_tc2 <- lm(tca2~., data=data.frame(hireoutput))
summary(fit_tc2) # cell type 3 - 0.000305 (3*)

# matching TCA celltype 3
fit_tc3 <- lm(tca3~., data=data.frame(hireoutput))
summary(fit_tc3) # cell type 3 - 0.0206 (*)

# matching TCA celltype 4
fit_tc4 <- lm(tca4~., data=data.frame(hireoutput))
summary(fit_tc4) # cell type 3 1.64e-05 (3*)

# matching TCA celltype 5
fit_tc5 <- lm(tca5~., data=data.frame(hireoutput))
summary(fit_tc5) # cell type 1 0.00114 (1*)

# matching TCA celltype 6
fit_tc6 <- lm(tca6~., data=data.frame(hireoutput))
summary(fit_tc6) # cell type 1 0.00468 (1*)

sink()


sink("comparetoref.txt")

# FOR COMPARING THE REFERENCE #
ref <- data.matrix(read.csv('10kref.csv'))
ref <- ref[,-1]
ref <- ref[,-1]


hire1 <- data.matrix(read.csv('celltype1forref.csv'))[,2]
fith1 <- lm(hire1~., data=data.frame(ref))
summary(fith1)


hire2 <- data.matrix(read.csv('celltype2forref.csv'))[,2]
fith2 <- lm(hire2~., data=data.frame(ref))
summary(fith2)

hire3 <- data.matrix(read.csv('celltype3forref.csv'))[,2]
fith3 <- lm(hire3~., data=data.frame(ref))
summary(fith3)


hire4 <- data.matrix(read.csv('celltype4forref.csv'))[,2]
fith4 <- lm(hire4~., data=data.frame(ref))
summary(fith4)

hire5 <- data.matrix(read.csv('celltype5forref.csv'))[,2]
fith5 <- lm(hire5~., data=data.frame(ref))
summary(fith5)


hire6 <- data.matrix(read.csv('celltype6forref.csv'))[,2]
fith6 <- lm(hire6~., data=data.frame(ref))
summary(fith6)


tca1 <- data.matrix(read.csv('tcacelltype1forref.csv'))[,2]
fitt1 <- lm(tca1~., data=data.frame(ref))
summary(fitt1)


tca2 <- data.matrix(read.csv('tcacelltype2forref.csv'))[,2]
fitt2 <- lm(tca2~., data=data.frame(ref))
summary(fitt2)

tca3 <- data.matrix(read.csv('tcacelltype3forref.csv'))[,2]
fitt3 <- lm(tca3~., data=data.frame(ref))
summary(fitt3)


tca4 <- data.matrix(read.csv('tcacelltype4forref.csv'))[,2]
fitt4 <- lm(tca4~., data=data.frame(ref))
summary(fitt4)


tca5 <- data.matrix(read.csv('tcacelltype5forref.csv'))[,2]
fitt5 <- lm(tca5~., data=data.frame(ref))
summary(fitt5)


tca6 <- data.matrix(read.csv('tcacelltype6forref.csv'))[,2]
fitt6 <- lm(tca6~., data=data.frame(ref))
summary(fitt6)



sink()

sink("ref_summaries.txt")
summary(fith1)
summary(fith2)
summary(fith3)
summary(fith4)
summary(fith5)
summary(fith6)
summary(fitt1)
summary(fitt2)
summary(fitt3)
summary(fitt4)
summary(fitt5)
summary(fitt6)
sink()











