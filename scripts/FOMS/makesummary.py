#!/usr/local/bin/python
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob, re, collections

'''create graphs from confidence intervals and partial aucs for all targets'''

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)    
mpl.rcParams['axes.color_cycle'] = tableau20[::2]

names = dict(
 cathg='CathG',
eralpha=r'ER$\alpha$',
eralpha_pot='ER$\\alpha$\nagonist',
erbeta=r'ER$\beta$',
fak=r'FAK',
fxia=r'FXIa',
hivrt=r'HIVrt',
hsp90=r'HSP90',
pka=r'PKA',
rho=r'Rho')

methodmap = dict(
foms='FOMS',
rdkit='RDKit',
vams='VAMS',
#fp="rdfp",
fp2='FP2')

values = collections.defaultdict(dict) #fpr - target - method - vals
fprs = set() #values partial aucs calculated at
for cfile in glob.glob('*/*.ci'):
    m = re.findall(r'(\S+)\/(\S+).\1\.f.*\.ci',cfile)
    if m:
        target = m[0][0]
        method = m[0][1]
        if method in methodmap:
            for line in open(cfile):
                (f, auc, l, u) = map(float,line.split())
                fprs.add(f)
                if target not in values[f]:
                    values[f][target] = dict()
                values[f][target][method] = (auc,l,u)
            
for f in fprs: #separate graph for each partial
    targnames = []
    data = collections.defaultdict(list) #indexed by method
    for targ in sorted(values[f].keys()):
        targnames.append(names[targ])
        for (method,vals) in values[f][targ].iteritems():
            data[method].append(vals)
    methodnames = sorted(data.keys(),key=lambda x: x if not x.startswith('fp') else 'zzz'+x)
    mn = len(methodnames)
    width = 1.0/(mn+1)
    ind = np.arange(len(targnames))
    for (i,method) in enumerate(methodnames):
        d = np.array(data[method])
        means = d[:,0]
        low = means-d[:,1]
        high = d[:,2]-means
        err = np.vstack((low,high))
        err[err < 0] = 0
        plt.bar(ind+i*width+0.5*width, d[:,0],width,yerr=err,color=tableau20[2*i],ecolor='k',label=methodmap[method],linewidth=0.5)
    plt.legend(loc='best',fontsize=18,ncol=2)
    plt.xticks(ind+0.5,targnames,fontsize=14)
    plt.yticks(fontsize=14)
    plt.axhline(y=0.5,color='gray',linestyle='--')
    plt.ylabel("Area Under the Curve",fontsize=18)    
    plt.savefig('aucs_%.2f.pdf' % f,bbox_inches='tight')
    plt.close()
