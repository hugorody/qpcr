#!/usr/bin/python3

from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
import matplotlib.pyplot as plt
from scipy import stats
import itertools as it
import pandas as pd
import numpy as np
import string
import math

############################################################### PARSE INPUT FILE
DCTdict = {} #dict {'Sample1':[DCTfloat1,DCTfloat2,DCTfloatn...]}
qpcr = pd.read_csv("anova_r.csv")
for i in range(1,len(qpcr)):
    sample = qpcr.iloc[i][0]
    sampleDCTs = []
    for jj in range(1,len(qpcr.columns)):
        endogenous_control = qpcr.iloc[0][jj] #endogenous must be the index 0 of input
        gene_test = qpcr.iloc[i][jj]
        DCTval = 2 ** (endogenous_control - gene_test) #Delta Ct = Ct gene test – Ct endogenous control
        sampleDCTs.append(DCTval)
    DCTdict[sample] = sampleDCTs

################################################################################
ANOVA = {} # {'AP18':[f,p],'Cmu1':[f,p]} ...
MEANS = {} # {'AP18':['mean1','mean2','meann']} n = number of treatments
SDs = {} # {'AP18':['std1','std2','stdn']} n = number of treatments
SEs = {} # {'AP18':['error1','error2','errorn']} n = number of treatments
barletters = {} # {'Sample1': {1: ['a'], 2: ['a', 'b'], 3: ['b'], 4: ['d'],n...},'Sample2':...}
for i in DCTdict.items():
    sample = i[0]
    DCTvalues = i[1] #[float1,float2,float3,float4,n...]
    listanova = [] #[[float1,float2,float3],[float1,float2,float3],n...] n = number of treatments
    for j in range(0,len(DCTvalues),3):
        addlist = []
        addlist.append(DCTvalues[j])
        addlist.append(DCTvalues[j+1])
        addlist.append(DCTvalues[j+2])
        listanova.append(addlist)
    ###################################################################### ANOVA
    f_value, p_value = stats.f_oneway(*listanova)
    ANOVA[sample] = [str(f_value),str(p_value)]

    ############################################ FEED STATISTICS TO DICTIONARIES
    listmedia = []
    listdesvio = []
    listerro = []
    for m in listanova:
        listmedia.append(str(np.mean(m))) #Mean
        listdesvio.append(str(np.std(m))) #Standard deviation (SD)
        listerro.append(str(np.std(m) / math.sqrt(3))) #Standard error (SE), 3 because of three replicates
    MEANS[sample] = listmedia
    SDs[sample] = listdesvio
    SEs[sample] = listerro

    ################################################################# TUKEY TEST
    countt = 0
    tukeydata = [] #[float1,float2,float3,float4,float5,n...]
    tukeygroups = [] #[1,1,1,2,2,n...]
    treatments = [] #list of treatments [1,2,3,4,n...]
    for ix in listanova:
        countt += 1
        treatments.append(countt)
        for ixx in ix:
            tukeygroups.append(countt)
            tukeydata.append(ixx)

    mc = MultiComparison(tukeydata, tukeygroups)
    result = mc.tukeyhsd()
    #print (sample,"\n",result)
    ############################################################# LETTER DANCING
    #The idea is to set equal letters for groups found not significantly different
    #by the TUKEY test
    myletters = list(string.ascii_lowercase)
    groupPAIR = list(it.combinations(treatments, 2)) #all possible treatments pairs
    groupRELATION = list(result.reject) # True, the two group's means ARE significantly different
                                        # False, the two group's means ARE NOT significantly different
    elemletters = {}
    cl = 0
    for u in range(len(groupPAIR)):
        ele1 = groupPAIR[u][0]
        ele2 = groupPAIR[u][1]
        if str(groupRELATION[u]) == "False":
            if ele1 in elemletters and ele2 in elemletters:
                equal_letter = any(elem in elemletters[ele1] for elem in elemletters[ele2])
                #print (ele1,elemletters[ele1],ele2,elemletters[ele2],groupRELATION[u],equal_letter)
                if equal_letter == False:
                    cl += 1
                    ele12l = myletters[cl]
                    add1 = elemletters[ele1]
                    add2 = elemletters[ele2]
                    add1.append(ele12l)
                    add2.append(ele12l)
                    elemletters[ele1] = add1
                    elemletters[ele2] = add2
            elif ele1 not in elemletters and ele2 not in elemletters:
                elemletters[ele1] = [myletters[cl]]
                elemletters[ele2] = [myletters[cl]]
            elif ele1 in elemletters and ele2 not in elemletters:
                ele2l = elemletters[ele1][0]
                elemletters[ele2] = [ele2l]
        else: #if relation is True
            if ele1 in elemletters and ele2 in elemletters:
                if elemletters[ele1] == elemletters[ele2]: #se forem diferentes nao precisa add
                    cl += 1
                    elemletters[ele1] = [myletters[cl]]
                    for uu in range(len(groupPAIR)):
                        if str(groupRELATION[uu]) == "False" and ele1 in groupPAIR[uu]:
                            for elem in groupPAIR[uu]:
                                if elem != ele1:
                                    add1 = elemletters[elem]
                                    add1.append(myletters[cl])
                                    elemletters[elem] = add1
            if ele1 in elemletters and ele2 not in elemletters:
                cl += 1
                elemletters[ele2] = [myletters[cl]]
            if ele1 not in elemletters and ele2 in elemletters:
                cl += 1
                elemletters[ele1] = [myletters[cl]]
            if ele1 not in elemletters and ele2 not in elemletters:
                elemletters[ele1] = [myletters[cl]]
                cl += 1
                elemletters[ele2] = [myletters[cl]]
    barletters[sample] = elemletters
#print (barletters)
############################################################## OUTPUT STATISTICS
outputQPCR = open("calculoqpcr.csv","w")
writeindex = 0
for i in MEANS.items():
    mykey = i[0]
    myvalues = i[1]
    if writeindex == 0:
        DCTMean = "DCTMean\t" * len(myvalues)
        DCTstd = "DCTstd\t" * len(myvalues)
        DCTerror = "DCTError\t" * len(myvalues)
        outputQPCR.write("Sample\t" + DCTMean + DCTstd + DCTerror + "F\tpval\n")
        writeindex += 1
    outputQPCR.write(mykey + "\t" + "\t".join(myvalues) + "\t" +
        "\t".join(SDs[mykey]) + "\t" + "\t".join(SEs[mykey]) + "\t" +
        "\t".join(ANOVA[mykey]) + "\n")
outputQPCR.close()

######################################################################## BARPLOT
N = len(treatments) #number of bars in barplot
ind = np.arange(N) #the x locations for the bars
width = 0.35 # the width of the bars
factor = 0.3

fig = plt.figure(figsize=(16,10))
fig.subplots_adjust(hspace=0.4, wspace=0.4, left=0.07, bottom=0.04, right=0.95, top=0.96)
labele = []
for i in treatments:
    labele.append("T" + str(i))
labels = tuple(labele)

addn = 0
for i in MEANS.items():
    xval = tuple(map(float, i[1]))
    xyerr = tuple(map(float, SEs[i[0]]))
    ylim = max(xval) + (0.5 * max(xval))
    sample = i[0]
    if float(ANOVA[sample][1]) <= 0.05: #if ANOVA is significant, print * on title
        anoval = "*"
    else:
        anoval = ""
    addn += 1
    ax = fig.add_subplot(2, 3, addn) #two first numbers: number of rows and columns respectively
    rects1 = ax.bar(ind + factor, xval, width, color='darkseagreen', yerr=xyerr, ecolor='black') #control

    ax.set_alpha(1)
    ax.set_xticks(ind + factor + width/2)
    ax.set_xticklabels(labels,rotation=20,fontsize=12)
    ax.set_ylabel('DCT value',fontsize=12)
    ax.set_ylim([0,ylim])
    ax.set_title(sample + anoval, horizontalalignment='center', fontname='Verdana',
                        fontsize=14, fontstyle='normal', fontweight='bold')

    count = 0
    for i in rects1:
        height = i.get_height()
        x = i.get_x()
        ax.text(x + width/2, (xyerr[count] + height + (0.02 * (ylim - xyerr[count] + height))),
            str(round(height, 3)) + "\n" + "".join(barletters[sample][count + 1]), horizontalalignment ='center',fontsize=12)
        count +=1

plt.savefig('barplot.png',dpi=100,bbox_inches='tight')
#plt.show()
