# -*- coding: utf-8 -*-
import xlrd
import numpy
import numpy as np
import math
import numpy.linalg as LA
from sklearn.cluster import KMeans
import random
from sklearn import tree
a=open(r'.\data\Known disease-miRNA association number.xlsx') 
b=open(r'.\data\Disease semantic similarity matrix 1.xlsx')
c=open(r'.\data\Disease semantic similarity matrix 2.xlsx')
d=open(r'.\data\Disease semantic similarity weighting matrix.xlsx')
e=open(r'.\data\miRNA functional similarity matrix.xlsx')
f=open(r'.\data\miRNA functional similarity weighting matrix.xlsx')
g=open(r'.\data\miRNA number.xlsx')
h=open(r'.\data\disease number.xlsx')
SS=numpy.zeros((383,383))                 
A=numpy.zeros((383,495))                     
KD=numpy.zeros((383,383))                    
SD=numpy.zeros((383,383))                  
DWM=numpy.zeros((383,383))                
MWM=numpy.zeros((495,495))              
FS=numpy.zeros((495,495))                   
SM=numpy.zeros((495,495))            
KM=numpy.zeros((495,495))                    
D1=numpy.zeros(10950)   
D2=numpy.zeros(10950)  
D3=numpy.zeros(10950)  
D4=numpy.zeros(10950)    
D5=numpy.zeros(10950)  
D6=numpy.zeros(10950)  
D7=numpy.zeros(10950)  
D8=numpy.zeros(10950)  
D9=numpy.zeros(10950)  
D10=numpy.zeros(10950)  
D11=numpy.zeros(10950)  
D12=numpy.zeros(10950)  
D13=numpy.zeros(10950)  
D14=numpy.zeros(10950)   
D15=numpy.zeros(10950) 
D16=numpy.zeros(10950)
D17=numpy.zeros(10950)  
D18=numpy.zeros(10950)   
D19=numpy.zeros(10950)  
xlsx1=xlrd.open_workbook(r'.\data\Disease semantic similarity matrix 1.xlsx')
xlsx2=xlrd.open_workbook(r'.\data\Disease semantic similarity matrix 2.xlsx')
sheet1=xlsx1.sheets()[0]                    
sheet2=xlsx2.sheets()[0]            
for i in range(383):
    for j in range(383):
        s1=sheet1.row_values(i)            
        s2=sheet2.row_values(i)
        m=s1[j]                             
        n=s2[j]
        SS[i,j]=float(m+n)/2                                                   #Obtain disease semantic similarity SS
xlsx3=xlrd.open_workbook(r'.\data\Known disease-miRNA association number.xlsx')
sheet3=xlsx3.sheets()[0]                    
for i in range(5430):
    s3=sheet3.row_values(i)                
    m=int(s3[0])
    n=int(s3[1])
    A[n-1,m-1]=1                                                               #Obtain adjacency matrix A
xlsx4=xlrd.open_workbook(r'.\data\Disease semantic similarity weighting matrix.xlsx')
sheet4=xlsx4.sheets()[0]                  
for i in range(383):
    for j in range(383):
        s4=sheet4.row_values(i)            
        DWM[i,j]=s4[j]                                                         #Get disease semantic weighting matrix DJQ
xlsx5=xlrd.open_workbook(r'.\data\miRNA functional similarity weighting matrix.xlsx')
sheet5=xlsx5.sheets()[0]                  
for i in range(495):
    for j in range(495):
        s5=sheet5.row_values(i)              
        MWM[i,j]=s5[j]                                                         #Get miRNA functional similarity weighting matrix MJQ
xlsx6=xlrd.open_workbook(r'.\data\miRNA functional similarity matrix.xlsx')
sheet6=xlsx6.sheets()[0]                 
for i in range(495):
    for j in range(495):
        s6=sheet6.row_values(i)              
        FS[i,j]=s6[j]                                                          #Get miRNA functional similarity matrix FS
C=np.asmatrix(A)
gamd=383/(LA.norm(C,'fro')**2);
kd=np.mat(np.zeros((383,383)))
km=np.mat(np.zeros((495,495)))
D=C*C.T;
for i in range(383):
        for j in range(i,383):
            kd[j,i]=np.exp(-gamd*(D[i,i]+D[j,j]-2*D[i,j]))
kd=kd+kd.T-np.diag(np.diag(kd))
KD=np.asarray(kd)                                                              # Obtain Gaussian interaction profile kernel similarity for disease SD
kd=[]                            
SD = np.multiply(SS,DWM)+np.multiply(KD,(1-DWM))
SD=np.asarray(SD)
gamam = 495/(LA.norm(C,'fro')**2);
E=C.T*C;
for i in range(495):
    for j in range(i,495):
        km[i,j]=np.exp(-gamam*(E[i,i]+E[j,j]-2*E[i,j]))
km=km+km.T-np.diag(np.diag(km))
KM=np.asarray(km)
km=[]
SM = np.multiply(FS,MWM)+np.multiply(KM,(1-MWM))
SM=np.asarray(SM)                                                              #Obtain Gaussian interaction profile kernel similarity for miRNA SM
#K_mean clustering
major=[]                          
minor=[]                            
for x in range(383):                 
    for y in range(495):
        if A[x,y]==0:                 
            major.append((x,y))
        else:
            minor.append((x,y))                                                #Divide the samples into Major and Minor
kmeans=KMeans(n_clusters=23, random_state=0).fit(major)      
center=kmeans.cluster_centers_
center_x=[]
center_y=[]
for j in range(len(center)):
    center_x.append(center[j][0])
    center_y.append(center[j][1])
labels=kmeans.labels_
type1_x=[]
type1_y=[]
type2_x=[]
type2_y=[]
type3_x=[]
type3_y=[]
type4_x=[]
type4_y=[]
type5_x=[]
type5_y=[]
type6_x=[]
type6_y=[]
type7_x=[]
type7_y=[]
type8_x=[]
type8_y=[]
type9_x=[]
type9_y=[]
type10_x=[]
type10_y=[]
type11_x=[]
type11_y=[]
type12_x=[]
type12_y=[]
type13_x=[]
type13_y=[]
type14_x=[]
type14_y=[]
type15_x=[]
type15_y=[]
type16_x=[]
type16_y=[]
type17_x=[]
type17_y=[]
type18_x=[]
type18_y=[]
type19_x=[]
type19_y=[]
type20_x=[]
type20_y=[]
type21_x=[]
type21_y=[]
type22_x=[]
type22_y=[]
type23_x=[]
type23_y=[]
for i in range(len(labels)):
    if labels[i]==0:
        type1_x.append(major[i][0])
        type1_y.append(major[i][1])
    if labels[i]==1:
        type2_x.append(major[i][0])
        type2_y.append(major[i][1])
    if labels[i]==2:
        type3_x.append(major[i][0])
        type3_y.append(major[i][1])
    if labels[i]==3:
        type4_x.append(major[i][0])
        type4_y.append(major[i][1])
    if labels[i]==4:
        type5_x.append(major[i][0])
        type5_y.append(major[i][1])
    if labels[i]==5:
        type6_x.append(major[i][0])
        type6_y.append(major[i][1])
    if labels[i]==6:
        type7_x.append(major[i][0])
        type7_y.append(major[i][1])
    if labels[i]==7:
        type8_x.append(major[i][0])
        type8_y.append(major[i][1])
    if labels[i]==9:
        type10_x.append(major[i][0])
        type10_y.append(major[i][1])
    if labels[i]==10:
        type11_x.append(major[i][0])
        type11_y.append(major[i][1])
    if labels[i]==11:
        type12_x.append(major[i][0])
        type12_y.append(major[i][1])
    if labels[i]==12:
        type13_x.append(major[i][0])
        type13_y.append(major[i][1])
    if labels[i]==13:
        type14_x.append(major[i][0])
        type14_y.append(major[i][1])
    if labels[i]==14:
        type15_x.append(major[i][0])
        type15_y.append(major[i][1])
    if labels[i]==15:
        type16_x.append(major[i][0])
        type16_y.append(major[i][1])
    if labels[i]==16:
        type17_x.append(major[i][0])
        type17_y.append(major[i][1])
    if labels[i]==17:
        type18_x.append(major[i][0])
        type18_y.append(major[i][1])
    if labels[i]==18:
        type19_x.append(major[i][0])
        type19_y.append(major[i][1])
    if labels[i]==19:
        type20_x.append(major[i][0])
        type20_y.append(major[i][1])
    if labels[i]==20:
        type21_x.append(major[i][0])
        type21_y.append(major[i][1])
    if labels[i]==21:
        type22_x.append(major[i][0])
        type22_y.append(major[i][1])
    if labels[i]==22:
        type23_x.append(major[i][0])
        type23_y.append(major[i][1])
type=[[]]*23                                
mtype=[[]]*23                                      
dataSet=[]                                
for k1 in range(len(type1_x)):
    type[0].append((type1_x[k1],type1_y[k1]))      
for k2 in range(len(type2_x)):
    type[1].append((type2_x[k2],type2_y[k2]))
for k3 in range(len(type3_x)):
    type[2].append((type3_x[k3],type3_y[k3]))
for k4 in range(len(type4_x)):
    type[3].append((type4_x[k4],type4_y[k4]))
for k5 in range(len(type5_x)):
    type[4].append((type5_x[k5],type5_y[k5]))
for k6 in range(len(type6_x)):
    type[5].append((type6_x[k6],type6_y[k6]))
for k7 in range(len(type7_x)):
    type[6].append((type7_x[k7],type7_y[k7]))
for k8 in range(len(type8_x)):
    type[7].append((type8_x[k8],type8_y[k8]))     
for k9 in range(len(type9_x)):
    type[8].append((type9_x[k9],type9_y[k9]))
for k10 in range(len(type10_x)):
    type[9].append((type10_x[k10],type10_y[k10]))
for k11 in range(len(type11_x)):
    type[10].append((type11_x[k11],type11_y[k11]))
for k12 in range(len(type12_x)):
    type[11].append((type12_x[k12],type12_y[k12]))
for k13 in range(len(type13_x)):
    type[12].append((type13_x[k13],type13_y[k13]))
for k14 in range(len(type14_x)):
    type[13].append((type14_x[k14],type14_y[k14]))
for k15 in range(len(type15_x)):
    type[14].append((type15_x[k15],type15_y[k15]))     
for k16 in range(len(type16_x)):
    type[15].append((type16_x[k16],type16_y[k16]))
for k17 in range(len(type17_x)):
    type[16].append((type17_x[k17],type17_y[k17]))
for k18 in range(len(type18_x)):
    type[17].append((type18_x[k18],type18_y[k18]))
for k19 in range(len(type19_x)):
    type[18].append((type19_x[k19],type19_y[k19]))
for k20 in range(len(type20_x)):
    type[19].append((type20_x[k20],type20_y[k20]))
for k21 in range(len(type21_x)):
    type[20].append((type21_x[k21],type21_y[k21]))
for k22 in range(len(type22_x)):
    type[21].append((type22_x[k22],type22_y[k22]))
for k23 in range(len(type23_x)):
    type[22].append((type23_x[k23],type23_y[k23]))                             #Divide Major into 23 clusters by K-means clustering
for k in range(23):
    mtype[k]=random.sample(type[k],240)                                        #Randomly extract 240 samples from each cluster
for m2 in range(383):
    for n2 in range(495):
        for z2 in range(23):
            if (m2,n2) in mtype[z2]: 
                dataSet.append((m2,n2))                                        #Store the randomly extracted 23X240 samples in the dataSet
for m3 in range(383):
    for n3 in range(495):
        if A[m3,n3]==1:
            dataSet.append((m3,n3))                                            #Combine Major and Minor into training samples containing 10,950 samples
#Decision tree
sumy1=numpy.zeros(10950)
sumy2=numpy.zeros(10950)
sumy3=numpy.zeros(10950)
sumy4=numpy.zeros(10950)
sumy5=numpy.zeros(10950)
sumy6=numpy.zeros(10950)
sumy7=numpy.zeros(10950)
sumy8=numpy.zeros(10950)
sumy9=numpy.zeros(10950)
sumy10=numpy.zeros(10950)
sumy11=numpy.zeros(10950)
sumy12=numpy.zeros(10950)
sumy13=numpy.zeros(10950)
sumy14=numpy.zeros(10950)
sumy15=numpy.zeros(10950)
sumy16=numpy.zeros(10950) 
sumy17=numpy.zeros(10950)
sumy18=numpy.zeros(10950)
sumy19=numpy.zeros(10950)
sumy20=numpy.zeros(10950)
x=[]                                                  
x1=[]                                     
x2=[]
y=[]                                                                              
D=numpy.ones(10950)*1.0/10950.0                                                #Initialize the weight of the training sample
for xx in dataSet:
    q=SD[xx[0],:].tolist()+SM[xx[1],:].tolist()
    x.append(q)                                                         
    if (xx[0],xx[1]) in minor:
        y.append(1)
    else:
        y.append(0)                                     
ys=numpy.array(y)
clf1=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf1=clf1.fit(x,y,sample_weight=D)
v1=clf1.predict(x)                                                             #Predict the labels of training samples
vs1=numpy.array(v1)
sumy1[vs1!=ys]=1                                  
sumd1=(sumy1*D).sum()                                                          #Calculate the error function of the weak classifier
at1=math.log((1-sumd1)/sumd1)*0.5                                              #Calculate the weight of the weak classifier in the strong classifier
z1=2*((sumd1*(1-sumd1))**0.5)       
for n3 in range(10950):
    D1[n3]=D[n3]*numpy.exp(-at1*y[n3]*v1[n3])/z1                               #Update the weight of the training sample
clf2=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf2=clf2.fit(x,y,sample_weight=D1)
v2=clf2.predict(x)                                       
vs2=numpy.array(v2)
sumy2[vs2!=ys]=1
sumd2=(sumy2*D1).sum()
at2=math.log((1-sumd2)/sumd2)*0.5         
z2=2*((sumd2*(1-sumd2))**0.5)       
for n3 in range(10950):
    D2[n3]=D1[n3]*math.exp(-at2*y[n3]*v2[n3])/z2 
clf3=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)    
clf3=clf3.fit(x,y,sample_weight=D2)
v3=clf3.predict(x)              
vs3=numpy.array(v3)
sumy3[vs3!=ys]=1
sumd3=(sumy3*D2).sum()
at3=math.log((1-sumd3)/sumd3)*0.5         
z3=2*((sumd3*(1-sumd3))**0.5)       
for n3 in range(10950):
    D3[n3]=D2[n3]*math.exp(-at3*y[n3]*v3[n3])/z3
clf4=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf4=clf4.fit(x,y,sample_weight=D3)
v4=clf4.predict(x)                
vs4=numpy.array(v4)
sumy4[vs4!=ys]=1
sumd4=(sumy4*D3).sum()
at4=math.log((1-sumd4)/sumd4)*0.5        
z4=2*((sumd4*(1-sumd4))**0.5)       
for n3 in range(10950):
    D4[n3]=D3[n3]*math.exp(-at4*y[n3]*v4[n3])/z4 
clf5=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)    
clf5=clf5.fit(x,y,sample_weight=D4)
v5=clf5.predict(x)               
vs5=numpy.array(v5)
sumy5[vs5!=ys]=1
sumd5=(sumy5*D4).sum()
at5=math.log((1-sumd5)/sumd5)*0.5        
z5=2*((sumd5*(1-sumd5))**0.5)       
for n3 in range(10950):
    D5[n3]=D4[n3]*math.exp(-at5*y[n3]*v5[n3])/z5 
clf6=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf6=clf6.fit(x,y,sample_weight=D5)
v6=clf6.predict(x)                
vs6=numpy.array(v6)
sumy6[vs6!=ys]=1
sumd6=(sumy6*D5).sum()
at6=math.log((1-sumd6)/sumd6)*0.5 
z6=2*((sumd6*(1-sumd6))**0.5)       
for n3 in range(10950):
    D6[n3]=D5[n3]*math.exp(-at6*y[n3]*v6[n3])/z6 
clf7=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf7=clf7.fit(x,y,sample_weight=D6)
v7=clf7.predict(x)               
vs7=numpy.array(v7)
sumy7[vs7!=ys]=1
sumd7=(sumy7*D6).sum()
at7=math.log((1-sumd7)/sumd7)*0.5        
z7=2*((sumd7*(1-sumd7))**0.5)       
for n3 in range(10950):
    D7[n3]=D6[n3]*math.exp(-at7*y[n3]*v7[n3])/z7 
clf8=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf8=clf8.fit(x,y,sample_weight=D7)
v8=clf8.predict(x)                
vs8=numpy.array(v8)
sumy8[vs8!=ys]=1
sumd8=(sumy8*D7).sum()
at8=math.log((1-sumd8)/sumd8)*0.5
z8=2*((sumd8*(1-sumd8))**0.5)
for n3 in range(10950):
    D8[n3]=D7[n3]*math.exp(-at8*y[n3]*v8[n3])/z8
clf9=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf9=clf9.fit(x,y,sample_weight=D8)
v9=clf9.predict(x)                                     
vs9=numpy.array(v9)
sumy9[vs9!=ys]=1
sumd9=(sumy9*D8).sum()
at9=math.log((1-sumd9)/sumd9)*0.5            
z9=2*((sumd9*(1-sumd9))**0.5)       
for n3 in range(10950):
    D9[n3]=D8[n3]*math.exp(-at9*y[n3]*v9[n3])/z9           
clf10=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf10=clf10.fit(x,y,sample_weight=D9)
v10=clf10.predict(x)                                       
vs10=numpy.array(v10)
sumy10[vs10!=ys]=1
sumd10=(sumy10*D9).sum()
at10=math.log((1-sumd10)/sumd10)*0.5         
z10=2*((sumd10*(1-sumd10))**0.5)       
for n3 in range(10950):
    D10[n3]=D9[n3]*math.exp(-at10*y[n3]*v10[n3])/z10     
clf11=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf11=clf11.fit(x,y,sample_weight=D10)
v11=clf11.predict(x)              
vs11=numpy.array(v11)
sumy11[vs11!=ys]=1
sumd11=(sumy11*D10).sum()
at11=math.log((1-sumd11)/sumd11)*0.5         
z11=2*((sumd11*(1-sumd11))**0.5)       
for n3 in range(10950):
    D11[n3]=D10[n3]*math.exp(-at11*y[n3]*v11[n3])/z11 
clf12=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf12=clf12.fit(x,y,sample_weight=D11)
v12=clf12.predict(x)                
vs12=numpy.array(v12)
sumy12[vs12!=ys]=1
sumd12=(sumy12*D11).sum()
at12=math.log((1-sumd12)/sumd12)*0.5        
z12=2*((sumd12*(1-sumd12))**0.5)       
for n3 in range(10950):
    D12[n3]=D11[n3]*math.exp(-at12*y[n3]*v12[n3])/z12 
clf13=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf13=clf13.fit(x,y,sample_weight=D12)
v13=clf13.predict(x)               
vs13=numpy.array(v13)
sumy13[vs13!=ys]=1
sumd13=(sumy13*D12).sum()
at13=math.log((1-sumd13)/sumd13)*0.5        
z13=2*((sumd13*(1-sumd13))**0.5)       
for n3 in range(10950):
    D13[n3]=D12[n3]*math.exp(-at13*y[n3]*v13[n3])/z13    
clf14=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf14=clf14.fit(x,y,sample_weight=D13)
v14=clf14.predict(x)                
vs14=numpy.array(v14)
sumy14[vs14!=ys]=1
sumd14=(sumy14*D13).sum()
at14=math.log((1-sumd14)/sumd14)*0.5 
z14=2*((sumd14*(1-sumd14))**0.5)       
for n3 in range(10950):
    D14[n3]=D13[n3]*math.exp(-at14*y[n3]*v14[n3])/z14
clf15=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf15=clf15.fit(x,y,sample_weight=D14)
v15=clf15.predict(x)               
vs15=numpy.array(v15)
sumy15[vs15!=ys]=1
sumd15=(sumy15*D14).sum()
at15=math.log((1-sumd15)/sumd15)*0.5        
z15=2*((sumd15*(1-sumd15))**0.5)       
for n3 in range(10950):
    D15[n3]=D14[n3]*math.exp(-at15*y[n3]*v15[n3])/z15 
clf16=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf16=clf16.fit(x,y,sample_weight=D15)
v16=clf16.predict(x)                
vs16=numpy.array(v16)
sumy16[vs16!=ys]=1
sumd16=(sumy16*D15).sum()
at16=math.log((1-sumd16)/sumd16)*0.5                      
z16=2*((sumd16*(1-sumd16))**0.5)       
for n3 in range(10950):
    D16[n3]=D15[n3]*math.exp(-at16*y[n3]*v16[n3])/z16
clf17=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf17=clf17.fit(x,y,sample_weight=D16)
v17=clf17.predict(x)               
vs17=numpy.array(v17)
sumy17[vs17!=ys]=1
sumd17=(sumy17*D16).sum()
at17=math.log((1-sumd17)/sumd17)*0.5        
z17=2*((sumd17*(1-sumd17))**0.5)       
for n3 in range(10950):
    D17[n3]=D16[n3]*math.exp(-at17*y[n3]*v17[n3])/z17
clf18=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf18=clf18.fit(x,y,sample_weight=D17)
v18=clf18.predict(x)                
vs18=numpy.array(v18)
sumy18[vs18!=ys]=1
sumd18=(sumy18*D17).sum()
at18=math.log((1-sumd18)/sumd18)*0.5  
z18=2*((sumd18*(1-sumd18))**0.5)       
for n3 in range(10950):
    D18[n3]=D17[n3]*math.exp(-at18*y[n3]*v18[n3])/z18 
clf19=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf19=clf19.fit(x,y,sample_weight=D18)
v19=clf19.predict(x)                
vs19=numpy.array(v19)
sumy19[vs19!=ys]=1
sumd19=(sumy19*D18).sum()
at19=math.log((1-sumd19)/sumd19)*0.5  
z19=2*((sumd19*(1-sumd19))**0.5)       
for n3 in range(10950):
    D19[n3]=D18[n3]*math.exp(-at19*y[n3]*v19[n3])/z19
clf20=tree.DecisionTreeClassifier(max_depth=9,min_samples_leaf=5)
clf20=clf20.fit(x,y,sample_weight=D19)
v20=clf20.predict(x)                
vs20=numpy.array(v20)
sumy20[vs20!=ys]=1
sumd20=(sumy20*D19).sum()
at20=math.log((1-sumd20)/sumd20)*0.5                                           #End of training process
for yy in major:
    q1=SD[yy[0],:].tolist()+SM[yy[1],:].tolist()          
    x1.append(q1)  
fs=clf1.predict_proba(x1)*at1+clf2.predict_proba(x1)*at2+clf3.predict_proba(x1)*at3+clf4.predict_proba(x1)*at4+clf5.predict_proba(x1)*at5+clf6.predict_proba(x1)*at6+clf7.predict_proba(x1)*at7+clf8.predict_proba(x1)*at8+clf9.predict_proba(x1)*at9+clf10.predict_proba(x1)*at10+clf11.predict_proba(x1)*at11+clf12.predict_proba(x1)*at12+clf13.predict_proba(x1)*at13+clf14.predict_proba(x1)*at14+clf15.predict_proba(x1)*at15+clf16.predict_proba(x1)*at16+clf17.predict_proba(x1)*at17+clf18.predict_proba(x1)*at18+clf19.predict_proba(x1)*at19+clf20.predict_proba(x1)*at20
px1=fs[:,1].tolist()                                                           #Score all unknown samples                                                                                            
xlsx7=xlrd.open_workbook(r'.\data\disease number.xlsx')
xlsx8=xlrd.open_workbook(r'.\data\miRNA number.xlsx')
sheet7=xlsx7.sheets()[0]                                 
sheet8=xlsx8.sheets()[0]  
px1=numpy.matrix(px1)
Sampleindex=numpy.argsort(-px1).tolist()    
Sampleindex=Sampleindex[0]  
f=open(r'Prediction results for all unknown samples.txt','a+')
f.writelines(['disease','\t','miRNA','\t','Score','\n'])
f.close()
for i in range(184155):                     
    a=fs[:,1][Sampleindex[i]]
    s7=sheet7.row_values(major[Sampleindex[i]][0])                   
    s8=sheet8.row_values(major[Sampleindex[i]][1])                     
    f=open(r'Prediction results for all unknown samples.txt','a+')
    f.writelines([s7[1],'\t',s8[1],'\t',str(a),'\n'])
    f.close()                                                                  #Obtain the prediction results for all unknown samples
