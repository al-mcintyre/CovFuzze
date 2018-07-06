#!/usr/bin/env python
#a program to plot gene coverage in multiple samples with peak regions
#Alexa McIntyre, 2018
import numpy as np
import pysam 
import pandas as pd
import seaborn as sns
import matplotlib
import re
import os
from matplotlib import pyplot as plt
from matplotlib import cm
from pybedtools import BedTool
import sys

def fill_coverage(coverage,strand,gene_cov,expn,normalize):
    if strand == '-':
        gene_cov = gene_cov[::-1]
    if expn not in coverage:
        coverage[expn] = {}
        coverage[expn]['count'] = []
        coverage[expn]['total'] = []
        coverage[expn]['gene_ind'] = []
    coverage[expn]['total'].append(float(sum(gene_cov)))
    if normalize:
        coverage[expn]['count'].append(np.array(gene_cov)*(len(gene_cov)/coverage[expn]['total'][-1]))
    else:
        coverage[expn]['count'].append(gene_cov)
    coverage[expn]['gene_ind'] = range(len(gene_cov))
    return coverage

def plot_gene(gene,outpref,bams,bedname,gtfname,peakname,labels,nsubplots,normalize):
    cds_indices = []
    transcript_end = 0
    if os.path.isfile(gtfname):
        gtf = open(gtfname,'r')
        for line in gtf.readlines():
            tabbed = line.split('\t')
            csome,feature = tabbed[0],tabbed[2]
            strand = tabbed[6]
            start = str(int(tabbed[3])-1)
            end = tabbed[4]
            if feature == 'CDS':
                if strand == '-':
                    inds = (int(end)-1,int(start))
                else:
                    inds = (int(start),int(end)-1)
                cds_indices.append(inds)
    else:
        strand = '+'
    print '{} strand'.format(strand)

    print cds_indices
    coverage = {}
    bed = BedTool(bedname)
    for label,bamname in zip(labels,bams): #add replicates here
        bam = BedTool(bamname) 
        cov = bed.coverage(bam, d=True,stream=True)
        indices = {}
        gene_cov = []
        for i,line in enumerate(cov):
            csome, start, end, intind, nreads = str(line).split('\t')
            indices[int(start)+int(intind)-1] = i
            gene_cov.append(int(nreads))

        if strand == '-':
            length_exon = len(indices)
            for ind in indices:
                indices[ind] = length_exon - indices[ind] - 1
        coverage = fill_coverage(coverage,strand,gene_cov,label,normalize)
        print "coverage calculated for {}: {}...".format(label,','.join([str(x) for x in gene_cov[:10]]))

    peaks = []
    try:
        peaksfi = open(peakname,'r')
        for line in peaksfi.readlines():
            if strand == '+':
                csome, start, end = line.split('\t')[:3]
            else:
                csome, end, start = line.split('\t')[:3]
            start,end = int(start), int(end)
            try:
                peaks.append((indices[start],indices[min(end,max(indices.keys()))]))
                #print start,end,indices[start],indices[min(end,max(indices.keys()))]
            except KeyError:
                ind_s, ind_e = 0,0
                found_s,found_e = False, False
                while not found_s and ind_s == 0 and start < max(indices.keys()):
                    try: 
                        ind_s = indices[start]
                        found_s = True
                    except KeyError:
                        start += 1 
                while not found_e and ind_e == 0 and end > 0:
                    try:
                        ind_e = indices[end]
                        found_e = True
                    except KeyError:
                        end += -1
                if found_s and found_e: 
                    #print start,end,ind_s,ind_e
                    peaks.append((ind_s,ind_e))
    except TypeError:
        pass
    print "peaks added"

    sns.set_style("ticks")
    figp,axes = plt.subplots(figsize = (5,2*nsubplots), nrows = nsubplots,sharey=True,sharex=True)
    if nsubplots == 1:
        axes = [axes]
    print "{} subplots".format(len(axes))
    #pal1 = cm.ScalarMappable(sns.light_palette("navy", as_cmap=True, reverse=True)).to_rgba(range(len(labels)/2))
    #pal2 = cm.ScalarMappable(sns.light_palette("orange", reverse=True, as_cmap=True)).to_rgba(range(len(labels)/2))
    #colours = [(0,0,0,1)]*len(labels)
    colours = ['#000000','#2294a3']*nsubplots #,'#006e90','#f18f01']*nsubplots #,'#adcad6','#c2724d','#9883e5','#f76b49','#3a208e','#f24236'] #don't like the colours? change them here!
    #colours[::2] = pal1 
    #colours[1::2] = pal2
    #print colours
    gene_len = coverage[labels[-1]]['gene_ind'][-1]
    maxpeak = 10 
    sp = 0

    def get_unique(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]
    uniq_labels = get_unique(labels)
    npersubplot = len(uniq_labels)/nsubplots

    reps = False
    for i,label,col in zip(range(len(uniq_labels)),uniq_labels,colours):
        if len(axes) > 1:
            ax = axes[sp]
        else:
            ax = axes[0]
        xs = coverage[label]['gene_ind']
        if len(coverage[label]['count']) > 1:
            ys = coverage[label]['count']
            std = np.std(ys,axis=0)
            ys = np.mean(ys,axis=0)
            reps = True
        else:
            ys = coverage[label]['count'][0]
        if len(xs) == 1:
            xs = xs[0]
        ax.plot(xs,ys,label=label,color=col)
        if i >= (sp+1)*(npersubplot)-1:
            sp += 1
        if reps:
            #print ys[:10],std[:10]
            #maxpeak = max(maxpeak,ax.get_ylim()[1]) #max(ys+std))
            ax.fill_between(xs, ys-std, ys+std ,alpha=0.3, facecolor=col)
        maxpeak = max(maxpeak,ax.get_ylim()[1])
        #else:
            #maxpeak = ax.get_ylim()[1] #max([maxpeak]+ax.get_ylim()[1])
    for i,ax in enumerate(axes):
        for (s,e) in cds_indices:
            print s,e 
            rect = matplotlib.patches.Rectangle((indices[s],0), indices[e]-indices[s], maxpeak, angle=0.0, alpha = 0.1, color = '#a5abaf')
            ax.add_patch(rect)
        for peak in peaks:
            rect = matplotlib.patches.Rectangle((peak[0],0), peak[1]-peak[0], maxpeak, angle=0.0, alpha = 0.3, color = '#ffff00')
            ax.add_patch(rect)
        #ax.set_xlabel(gene)
        ax.xaxis.set_ticks(np.arange(0,gene_len,max(100*np.round(gene_len/500.,0),500)))
        ax.set_xlim([0,gene_len+1])
        #print maxpeak, int(0.05*maxpeak)
        ax.set_ylim([0,maxpeak]) #+int(0.05*maxpeak)])
        if normalize:
            ax.set_ylabel('normalized\ncoverage')
        else:
            ax.set_ylabel('coverage')
        if i == nsubplots-1:
            ax.set_xlabel(gene)
        ax.legend(loc=2)
    #sns.despine(fig=figp)
    plt.tight_layout()
    plt.savefig(outpref+'_'+gene+'_coverage.pdf',dpi=500,bbox_inches='tight')

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Plot gene coverage for multiple experiments')
    parser.add_argument('-g','--gene',type=str,required=False,help='gene name (default = GeneDoe)',default="GeneDoe")
    parser.add_argument('-o','--out',type=str,required=False,help='output prefix (incl. dir)',default="gene_coverage_plots/Sample1")
    parser.add_argument('--bams',nargs='+',required=True,help='bam files (separated by spaces)')
    parser.add_argument('--bed',type=str,required=True,help='bed file with region(s) for which to calculate coverage')
    parser.add_argument('--gtf',type=str,required=False,help='gtf file for gene',default='imagtf')
    parser.add_argument('-p','--peaks',type=str,required=False,help='bed file with peaks',default=None)
    parser.add_argument('-l','--labels',nargs='+',required=True,help='labels associated with bams - if replicates, use same labels')
    parser.add_argument('-n','--nsubplots',type=int,required=False,help='number of subplots- bams will be split evenly based on the order given (default = 1)',default=1)
    parser.add_argument('--normalize',action='store_true',required=False,help='normalize by gene length/summed coverage (default = False)',default=False)
    parser.add_argument('-v','--version',action='version',version='%(prog)s (v0.1.2)') 

    args = parser.parse_args()
    outdir = (args.out).split('/')
    if len(outdir) > 1:
        outdir = outdir[0]
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
    assert len(args.labels) == len(args.bams), 'mismatch between length of bam list and length of conditions'
    for fi in args.bams + [args.bed]: 
        assert os.path.isfile(fi),'no file found at {}'.format([fi])
    print 'sample {}: plotting coverage for {}...'.format(args.out,args.gene)
    plot_gene(args.gene,args.out,args.bams,args.bed,args.gtf,args.peaks,args.labels,args.nsubplots,args.normalize)

if __name__ == "__main__":
    main()
