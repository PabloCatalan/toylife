#%%toyLIFE patterns figures for paper
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
import pylab
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mticker

#toyLIFE PARAMETERS
time=100
cells=31
ring=0
inputs=['zero', 'linear', 'linear2', 'exp', 'gauss', 'point']
morphogen='_'+inputs[5]
D=['only_P1', 'only_P2', 'only_P1_P2', 'also_dimer']
diff=D[0]

#COLOR PARAMETERS
color0='white'
color1='#30678d'
color2='#e06a00'
color3='#4d4d4d'

#READ ALL PATTERNS
with open('results/phen_id_'+str(cells)+'_'+str(time)+'_'+str(ring)+str(morphogen)+'_'+diff+'.txt', 'r') as f:
    read_data=f.read()
rows=len(read_data.split())
phens=int(rows/(time+2))
s=read_data.split()
dict1={}
for p in range(phens):
    name=s[p*(time+2)]
    size=s[p*(time+2)+1]
    z1=[s[p*(time+2)+i+2] for i in range(time)]     
    dict1[name]=[z1,size]
    
def plot_pattern_list(P, figName, dict1):
    color0='white'
    color1='#30678d'
    color2='#e06a00'
    color3='#4d4d4d'
    maxrows=1
    maxcols=len(P)
    fig=plt.figure(figsize=(maxcols*4,maxrows*6))
    Lab=['a)', 'b)', 'c)', 'd)']
    for i1, pat in enumerate(P):
        ax=fig.add_subplot(maxrows,maxcols,i1+1)
        name=str(pat)
        f1=dict1[name][0]
        data=[[0 for z in range(cells)] for y in range(time)]
        for i in range(time):
            f=f1[i]
            for j in range(cells):
                data[i][j]=float(f[j])/3.0
    
        cmap = ListedColormap([color0, color1, color2, color3], 'indexed')
        heatmap = ax.pcolormesh(data, cmap=cmap, vmin=0.0, vmax=1.0, edgecolor='face')
        ax.set_xlim([0,cells])
        ax.set_ylim([0,time])  
        ax.invert_yaxis()
        ax.text(-0.4, 1.10, Lab[i1], fontsize=40, transform = ax.transAxes)        #plt.xlabel('gene number', family='sans-serif', weight='bold', size='xx-large')
        ax.set_ylabel('time', fontsize=30)
        ax.set_xlabel('cells', fontsize=30)
        ax.set_yticklabels(ax.get_yticks(), fontsize=20)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))  
        ax.set_xticks([])
        fig.tight_layout()
        fig.savefig('results/'+figName, bbox_inches='tight')

#%%FIGURE 3
P=[11,108,21,3]
plot_pattern_list(P, 'fig3.pdf', dict1)

#FIGURE 5a
P=[113,109,170]
plot_pattern_list(P, 'fig5a.pdf', dict1)

#%%FIGURE 4E AND SUPP FIG10
from scipy.optimize import curve_fit

#GAUSSIAN
def exp(x, a, b, c):
    return a*np.exp(-b*(x-c)**2)

#PLOT
def plot_abundance(AB, yLab, figName):
    n, bins=np.histogram(AB, bins=50)
    mask=n>0
    n1=n[mask]
    w=bins[1]-bins[0]
    b1=bins[:-1][mask]+w/2

    #FIT TO GAUSSIAN
    popt, pcov=curve_fit(exp, b1, n1)
    c1='#30678dff'
    c2='#e06a00ff'
    
    #COMPUTE R SQUARED
    pred=exp(b1,*popt)
    eT=n1-n1.mean()
    ssT=np.sum(eT**2)
    eR=n1-pred
    ssE=np.sum(eR**2)
    R2=1-ssE/ssT
    R2=np.round(R2*100)/100

    #PLOT
    fig=plt.figure(figsize=(4.5,4))
    ax=plt.gca()
    ax.plot(10**b1, n1, 'o', color=c1)
    xp=np.linspace(0,12,100)
    ax.plot(10**xp, exp(xp, *popt), color=c2)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('S', size=15)
    ax.set_ylabel(yLab, size=15)
    ax.set_xticklabels(ax.get_xticks(), size=10)
    ax.set_yticklabels(ax.get_yticks(), size=10)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(g))
    ax.set_ylim([0.8*min(n1),1.2*max(n1)])
    ax.set_xlim([10**0.5,10**12.5])
    ax.text(0.05, 1.05, 'R$^2$='+str(R2), transform=ax.transAxes)
    fig.savefig('results/'+figName, bbox_inches='tight')
    

#READ GRN ABUNDANCE
grnAB=[]
with open('results/grn_abundance.txt', 'r') as f:
    for line in f:
        w=line.split()
        grnAB.append(np.log10(float(w[1])))
plot_abundance(grnAB, 'number of GRNs', 'fig4e.pdf')

#READ CELLULAR AUTOMATA ABUNDANCE
caAB=[]
with open('results/ca_abundance.txt', 'r') as f:
    for line in f:
        w=line.split()
        caAB.append(np.log10(float(w[1])))
plot_abundance(caAB, 'number of CA', 'suppfig10a.pdf')

phenAB=[]
with open('results/phen_size_31_100_0_point_only_P1.txt', 'r') as f:
    for line in f:
        w=line.split()
        phenAB.append(np.log10(float(w[1])))
plot_abundance(phenAB, 'number of patterns', 'suppfig10b.pdf')

#%%FIG5d
def par(x,a,b,c):
    return a+b*x+c*x**2

def sig(x,a,b,c):
    return 1/(1+a*np.exp(-b*(x-c)))

#READ DATA
import pandas as pd
DF=pd.read_csv('results/simulations_prob_of_arrival_patterns.txt',
               sep='\t',names=['Pattern', 'Pavg', 'Pstd', 'Tavg',
                               'Tstd', 'S'])
N=np.sqrt(1000)
#FIRST SUBPLOT
fig=plt.figure()
ax=fig.add_subplot(1,2,1)
x=np.array(np.log10(DF['S']).tolist())
y1=np.array(DF['Pavg'].tolist())
y1err=np.array((DF['Pstd']/N).tolist())
y2=np.array(DF['Tavg'].tolist())
y2err=np.array((DF['Tstd']/N).tolist())
#y2err2=np.array((DF['Tstd']).tolist())
y3=np.array(np.log10(DF['Tavg']).tolist())
y3err=np.array((np.log10(DF['Tstd'])/N).tolist())
#y3=np.array(y2)
c1='#30678dff'
c2='#e06a00ff'
ax.errorbar(10**x,y1, yerr=y1err, fmt='o', c=c1)
ax.set_xlabel('S', size=15)
ax.set_ylabel('fixation probability', size=15)
ax.set_xticks(10**np.array([2,4,6,8,10,12]))
ax.set_xticklabels(ax.get_xticks(), size=10)
ax.set_yticklabels(ax.get_yticks(), size=10)
ax.set_xscale('log')
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#FIT
mask = ~np.isnan(x) & ~np.isnan(y1)
popt, pcov=curve_fit(sig, x[mask], y1[mask], p0=[1,1,0])
xp=np.linspace(min(x),max(x),100)
yp=sig(xp,*popt)
#COMPUTE R SQUARED
pred=sig(x[mask],*popt)
eT=y1[mask]-y1[mask].mean()
ssT=np.sum(eT**2)
eR=y1[mask]-pred
ssE=np.sum(eR**2)
R2=1-ssE/ssT
R2=np.round(R2*100)/100
ax.plot(10**xp,yp,color=c2)
ax.text(0.05, 1.05, 'R2='+str(R2), transform=ax.transAxes)
ax.set_xlim([10**1.8,10**12.2])
ax.text(-0.1, 1.05, 'a)', transform=ax.transAxes, fontsize=15)

#SECOND SUBPLOT
ax=fig.add_subplot(1,2,2)
ax.errorbar(10**x,y2,yerr=y2err, fmt='o', c=c1)
ax.set_xlabel('S', size=15)
ax.set_ylabel('fixation time', size=15)
ax.set_xticks(10**np.array([2,4,6,8,10,12]))
ax.set_xticklabels(ax.get_xticks(), size=10)
ax.set_yticklabels(ax.get_yticks(), size=10)
ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.set_xlim([10**1.8,10**12.2])
ax.set_ylim([1,100010])
ax.set_xscale('log')
ax.set_yscale('log')
ax.text(-0.1, 1.05, 'b)', transform=ax.transAxes, fontsize=15)
#FIT
mask = ~np.isnan(x) & ~np.isnan(y3)
popt, pcov=curve_fit(par, x[mask], y3[mask], #sigma=y2err[mask], 
                     p0=[3.8, 0.1, -0.1])
xp=np.linspace(min(x),max(x),100)
yp=par(xp,*popt)
#COMPUTE R SQUARED
pred=par(x[mask],*popt)
eT=y3[mask]-y3[mask].mean()
ssT=np.sum(eT**2)
eR=y3[mask]-pred
ssE=np.sum(eR**2)
R2=1-ssE/ssT
R2=np.round(R2*100)/100
ax.plot(10**xp,10**yp,color=c2)
ax.text(0.05, 1.05, 'R2='+str(R2), transform=ax.transAxes)
fig.tight_layout()
fig.savefig('results/fig5d.pdf', bbox_inches='tight')

#%%SUPPFIG S8
import gc
import itertools

def plot_all_patterns(tissueS):
    time=100
    cells=tissueS
    ring=0
    inputs=['pointleft', 'linearP2', 'point', 'linear']
    inputs=[inputs[2]]
    D=['only_P1', 'only_P2', 'only_P1_P2', 'also_dimer']
    D=['only_P1']
    for (i1, diff) in itertools.product(inputs, D):
        morphogen='_'+i1
        with open('results/phen_id_'+str(cells)+'_'+str(time)+'_'+str(ring)+str(morphogen)+'_'+diff+'.txt', 'r') as f:
            read_data=f.read()
    blocksize=time+2
    rows=len(read_data.split())
    phens=int(rows/(blocksize))
    s=read_data.split()
    maxrows=10
    maxcols=10
    itph=10*10
    L=phens//itph+1
    for figura in range(L):
        plt.clf()
        fig, axarr = plt.subplots(maxrows, maxcols, figsize=(maxcols*2,maxrows*3))
        rows=0
        columns=0
        for p in range(figura*itph, min([figura*itph+itph, phens])):
            data=[[0 for z in range(cells)] for y in range(time)]
            name=s[p*(blocksize)]
            dsize=s[p*(blocksize)+1]
            for i in range(time):
                f=s[p*(blocksize)+i+2]
                for j in range(cells):
                    data[i][j]=float(f[j])/3.0
            a1=rows%maxrows
            a2=columns%maxcols
            ax=axarr[a1][a2]
            columns+=1
            if (columns%maxcols)==0:
                rows+=1  
                
            color0='white'
            color1='#30678d'
            color2='#e06a00'
            color3='#4d4d4d'
            
            cmap = ListedColormap([color0, color1, color2, color3], 'indexed')
            heatmap = ax.pcolormesh(data, cmap=cmap, vmin=0.0, vmax=1.0, rasterized=True)
            xticks = ax.xaxis.get_major_ticks()
            for i in range(len(xticks)):
                xticks[i].label1.set_visible(False)
            yticks = ax.yaxis.get_major_ticks()
            for i in range(len(yticks)):
               yticks[i].label1.set_visible(False)
            ax.set_xlim([0,cells])
            ax.set_ylim([0,time])  
            ax.invert_yaxis()
            ax.text(0, 1.01, name+'_'+dsize, fontsize=10, transform = ax.transAxes)        #plt.xlabel('gene number', family='sans-serif', weight='bold', size='xx-large')
        fig.savefig('results/supffigS8_tissuesize_'+str(cells)+'_'+str(figura)+'.pdf', bbox_inches='tight')
        fig.clf()
        plt.close('all')#pylab.close()
    del read_data
    del s
    gc.collect()
    
plot_all_patterns(31)
plot_all_patterns(51)

#%%SUPPFIG S11
def line(x, a, b):
   return a+b*x

def plot_complexity(phen):
    fig=plt.figure(figsize=(3,3))
    color1='#30678dff'
    color2='#e06a00ff'
    DF=pd.read_csv('results/complexity_'+phen+'.txt',
               sep='\t',names=['C', 'S'])
    ax=plt.gca()
    x=np.array(DF['C'])
    y=np.array(DF['S'])/DF['S'].sum()
    ax.plot(x,y,'o',color=color1)
    #THEORETICAL UPPER BOUND
    N0=len(DF)
    maxK=x.max()
    aslope=np.log2(N0)/maxK
    xp=np.linspace(x.min(), x.max(), 100)
    yp=2**(-aslope*xp)
    ax.plot(xp,yp,color=color2)
    ax.set_xlabel('$\~K$', size=20)
    ax.set_ylabel('S', size=20)
    ax.set_xticklabels(ax.get_xticks(), size=10)
    ax.set_yticklabels(ax.get_yticks(), size=10)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.set_yscale('log')
    ax.text(-0.1, 1.05, phen, transform=ax.transAxes, fontsize=15)
    fig.savefig('results/suppfigS11_'+phen+'.pdf', bbox_inches='tight')
    
plot_complexity('grn')
plot_complexity('ca')
plot_complexity('pattern')
#%%SUPPFIG S12

#THEORY
T={}
with open('results/simulations_patternFig3_theory.txt', 'r') as f:
    for line in f:
        w=line.split()
        T[w[0]]=float(w[2])
#SIMULATIONS
S={}
with open('results/simulations_patternFig3.txt', 'r') as f:
    for line in f:
        w=line.split()
        S[w[0]]=float(w[1])

x=[]
y=[]
for key, value in S.items():
    if key in T:
        x.append(T[key])
        y.append(value)
Labels=['GRN I', 'GRN III', 'GRN VII', 'GRN IX', 
        'GRN XI', 'GRN XV', 'GRN XVI']
    
def line(x, a, b):
   return a+b*x

fig=plt.figure()
color1='#30678dff'
color2='#e06a00ff'    
ax=plt.gca()
ax.plot(x,y,'o', color=color1)
#LABELS
for i,l in enumerate(Labels):
    ax.text(x[i]*10**(0.05),y[i]/10**(0.1),l,size=10)
ax.set_xlabel('relative abundance', size=20)
ax.set_ylabel('observed frequency', size=20)
ax.set_xticklabels(ax.get_xticks(), size=10)
ax.set_yticklabels(ax.get_yticks(), size=10)
#FIT
popt, pcov=curve_fit(line, x, y)
xp=np.linspace(0,1,10000)
yp=line(xp,*popt)
yp2=line(xp,0,1)
ax.plot(xp,yp,color=color2)
ax.set_yscale('log')
ax.set_xscale('log')
#COMPUTE R SQUARED
pred=line(np.array(x),*popt)
eT=y-np.average(y)
ssT=np.sum(eT**2)
eR=y-pred
ssE=np.sum(eR**2)
R2=1-ssE/ssT
R2=np.round(R2*100)/100
ax.text(0.05,1.05,'R2='+str(R2),transform=ax.transAxes)
fig.savefig('results/suppfigS12.pdf', bbox_inches='tight')