import pandas as pd
import glob
import numpy as np
# from scipy import stats
# from statsmodels.stats import multitest
# import seaborn as sns
# import matplotlib.pyplot as plt


cpath ="*.gz"
samples = pd.DataFrame()
for i, filename in enumerate(glob.glob(cpath)):
    colV = filename.strip('.txt.gz')
    samples2 = pd.DataFrame()

    df = pd.read_csv(filename, compression='gzip',sep='\t')
    if 'GeneID' not in samples.columns:
        samples = df[['GeneID','RPKM']]
        samples.columns=['GeneID',colV]
    else:
        samples2 = df[['GeneID','RPKM']]
        samples2.columns=['GeneID',colV]

        samples = pd.merge(samples,samples2, on='GeneID', )

print("\n",'List of file names read to dataframe:','\n')
samples.set_index("GeneID", inplace = True) 
print(samples.head())

samples.to_csv('samples_count_matrix.csv',)
samples_annotated = samples
#############################################################
import rpy2.robjects as robjects 
from rpy2.robjects.packages import importr 
ALL = importr('ALL') 
limma = importr('limma')
exprs = robjects.r['exprs']
summary = robjects.r['summary']
matrix = robjects.r['as.matrix']
new = robjects.r['new']
robjects.r('data("ALL")')
data = robjects.globalenv['ALL'] 
featureNames = robjects.r['featureNames'] 
ExpressionSet = robjects.r['ExpressionSet']
character = robjects.r['as.character']
pas = robjects.r['paste']
fac = robjects.r['as.factor']
mmax = robjects.r['model.matrix']


from rpy2.robjects import pandas2ri
pandas2ri.activate()

print(samples.head())
samples_log2=samples.apply(np.log2)
print(samples_log2.head())

import rpy2.robjects as ro

# print(asda)
R = ro.r

control_sample = input('Please enter column numbers of control samples (ex:0,2,4): ')
control_samples = control_sample.split(',')
control_samples = [int(i) for i in control_samples]
Test_sample =  input('Please enter column numbers of test samples (ex:3,5,7): ')
Test_samples = Test_sample.split(',')
Test_samples = [int(i) for i in Test_samples]

f_samples = pd.DataFrame() 
f_samples2 = samples_log2.iloc[:,Test_samples] 
f_samples1 = samples_log2.iloc[:,control_samples]
f_samples = f_samples1.join(f_samples2)

r_samples = pandas2ri.py2ri(f_samples)

asda = new("ExpressionSet", exprs=matrix(r_samples))


a_labels=control_samples + Test_samples
# print(a_labels)
a_lab=list()
for i, t in enumerate(a_labels):
    if i < len(control_samples):
        a_lab.append("1")
    else:
        a_lab.append("0")
sml = pas("G",a_lab,sep="")
fl = fac(sml)
# print(fl)
R.assign('fl',fl)
R.assign('r_samples',asda)
R('design<-model.matrix(~ fl + 0,r_samples)')
R('print(design)')
R('colnames(design)<-levels(fl)')
R('print(design)')
R('fit<-lmFit(r_samples,design)')
R('cont.matrix <- makeContrasts(G1-G0, levels=design)')
R('fit2 <- contrasts.fit(fit, cont.matrix)')
R('fit2 <- eBayes(fit2, 0.01)')
R('tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)')
R('tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC"))')
# R('print(tT)')

R('write.csv(tT, file="F_matrix.csv", sep=",")')
