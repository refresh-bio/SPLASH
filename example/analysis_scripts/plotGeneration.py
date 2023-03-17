import numpy as np
import os
import glob
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator
import time
import scipy
from scipy.sparse import find
from bisect import bisect_left, bisect_right
import argparse
import seaborn as sns

import colorlover as cl




### command line inputs
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dsName",
        type=str,
        help='Dataset name to be used in title of plots.',
    )
    parser.add_argument(
        "--satcFolder",
        type=str,
        help='Folder containing .satc files.'
    )
    parser.add_argument(
        "--metadataPath",
        type=str,
        help='.tsv file containing two columns with names [\'sampleName\', \'metadata\'].'
    )
    parser.add_argument(
        "--outFolder",
        type=str,
        help='Output folder to save files to.'
    )
    parser.add_argument(
        "--pvDfPath",
        type=str,
        help='Path to p-value .tsv file.'
    )
    parser.add_argument(
        "--anchorFile",
        type=str,
        help='File of anchors to be plotted, one anchor per line, no header.'
    )
    parser.add_argument(
        "--skipSATC",
        default=False,
        action='store_true',
        help='Boolean: if rerunning, can skip SATC dump step, otherwise necessary. Flag is optional, if included will skip SATC generation. To run SATC generation, do not include --skipSATC flag'
    )
    
    parser.add_argument(
        "--satc_dump_file",
        type=str,
        help='Path to satc_dump utility file (within NOMAD2, /bin/satc_dump)'
    )

    args = parser.parse_args()
    return args



### parse satc files to dump files
def parseSATC():
    print('parsing SATC')
    Path(outFolder+'/satc_unpacked/all').mkdir(parents=True, exist_ok=True)
    
    for fname in tqdm(glob.glob(satcFolder+'/result.bin*.satc')):
        outFile = outFolder+'/satc_unpacked/'+fname.split('_')[-1]+'.dump'
        
        cmd=f"{satc_dump_file} --anchor_list {anchorFile} {fname} {outFile}"
        print(cmd)
        os.system(cmd)


### generate counts dataframe from dumped files
def generateCtsDf():
    dfArr = []
    for fname in tqdm(glob.glob(outFolder+'/satc_unpacked/all/*.dump'), desc='reading .dump files'):
        dfArr.append(pd.read_csv(fname,names=['id','anchor','target','counts'],sep='\t'))
    
    ctsDf = pd.concat(dfArr)

    id_to_sample_mapping = pd.read_csv(satcFolder+'/sample_name_to_id.mapping.txt',names=['sample','id'],sep=' ')

    ctsDf = pd.merge(ctsDf,id_to_sample_mapping).drop(columns='id')
    ctsDf = ctsDf.groupby(['anchor','target','sample']).counts.sum().reset_index()
    ctsDf.to_csv(outFolder+'/countsDf.tsv',sep='\t')
    return ctsDf
        
        
        
#### for target plotting
def bpToInt(x):
    if x=='A':
        return 0
    if x=='T':
        return 1
    if x=='C':
        return 2
    if x=='G':
        return 3
    return 4 #### for nan

dna_colors = ['#ffffff','#8dd3c7','#ffffb3','#fb8072','#80b1d3'] #["N","A","T","C","G"]
col_dict={0:dna_colors[1],
  1:dna_colors[2],
  2:dna_colors[3],
  3:dna_colors[4]}
labels = np.array(["A","T","C","G"])
len_lab = len(labels)
cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

norm_bins = np.sort([*col_dict.keys()]) + 0.5
norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

norm = mpl.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
diff = norm_bins[1:] - norm_bins[:-1]
tickz = norm_bins[:-1] + diff / 2
        
    
### main function for generating plots
def generatePlots(ctsDf):
    Path(outFolder+'/plots/').mkdir(parents=True, exist_ok=True)
    
    
    anchSet = set(ctsDf.anchor.to_list())
    ctsDf = ctsDf.set_index('anchor')
    
    metadataDf = pd.read_csv(metadataPath)
    
    ### if more than one anchor, perform clustering on anchors to more easily parse through list of outputted anchors
    if len(anchSet)<=1:
        clustMapping = {list(anchSet)[0]:0}
    else:
        clustMapping = clusterAnchors(list(anchSet))
    
    print('reading p-values')
    pvDf = pd.read_csv(pvDfPath,sep='\t').set_index('anchor')

    
    print('starting loop')
    for anch in tqdm(ctsDf.index.unique(),desc='plotting loop over anchors'):
        if anch not in anchSet: ### should never reach
            print(f'skipping {anch}')
            continue
            
        cont_table = (ctsDf.loc[anch]
                            .pivot(index='target', columns='sample', values='counts')
                            .fillna(0)
                            .reset_index()
                            )
        targets = cont_table.target.to_numpy()
        samples = cont_table.columns[1:]
        cont_mat = cont_table.drop(columns='target').to_numpy()

        ## reorder targets to go in descending order
        targOrdering = np.argsort(-cont_mat.sum(axis=1))
        cont_mat = cont_mat[targOrdering]
        targets = targets[targOrdering]

        ## filter down to fewer rows
        sortedTargCounts = cont_mat.sum(axis=1)
        thresh = max(sortedTargCounts[min(10,len(targOrdering)-1)],
                     5,sortedTargCounts[0]/20) ## for large datasets
        thresh = min(thresh,sortedTargCounts[1]) ### need at least 2 targets
        relevantRows = cont_mat.sum(axis=1)>=thresh

        ## filter out columns with fewer than 5 counts
        relevantCols = cont_mat.sum(axis=0)>=5
        cont_mat = cont_mat[:,relevantCols]
        samples = samples[relevantCols]
        colSums = cont_mat.sum(axis=0)


        ## filter out rows with fewer than thresh counts
        cont_mat = cont_mat[relevantRows]
        targets = targets[relevantRows]
        
        ## reorder columns
        toSort = list(zip(metadataDf.set_index('sampleName').loc[samples].metadata,-cont_mat[0]/colSums))
        col_ordering = [x[1] for x in sorted((e,i) for i,e in enumerate(toSort))]
        cont_mat = cont_mat[:,col_ordering]
        samples = samples[col_ordering]
        colSums = colSums[col_ordering]
        
        
        ### adjust target mat
        targMat = np.zeros((len(targets),len(targets[0])))
        for i,targ in enumerate(targets):
            targMat[i,:] = [bpToInt(x) for x in targ]


        fig,axs = plt.subplots(nrows=2,ncols=2)#,gridspec_kw={'width_ratios': [1, 2]})
        fig.set_size_inches(10, 10)
        fig.set_facecolor('white')
        

        paspect = 'auto'
        plt.sca(axs[1,0])

        im= plt.imshow(cont_mat/ colSums,aspect=paspect,interpolation='nearest')

        nj = colSums.astype('int')


        cellType = metadataDf.set_index('sampleName').loc[samples].metadata.values

        labels = [cellType[i] + f' ({nj[i]})' for i in range(len(nj))]

        
        ##### if there are not too many samples, can include labels
        ### option 1: cellType + nj
#         plt.xticks(ticks=np.arange(cont_mat.shape[1])
#                    ,labels=labels,
#                    rotation=-80 ,ha="left", rotation_mode="anchor")
        ### option 2: nj
#         plt.xticks(ticks=np.arange(cont_mat.shape[1])
#                    ,labels=nj,
#                    rotation=-80 ,ha="left", rotation_mode="anchor")

        plt.ylabel('Target')
        plt.yticks(ticks=range(cont_mat.shape[0]), labels=cont_mat.sum(axis=1).astype(int))

        plt.colorbar(orientation="horizontal",pad=.2)
        
        
        tmpDf = pd.DataFrame({'sampName':samples, 't1Dominant':(cont_mat/colSums)[0]>.5, 'celltype': cellType})
        perfectConcordance = tmpDf.groupby('celltype').t1Dominant.nunique().max()==1
        


        plt.sca(axs[1,1])
        im = plt.imshow(targMat, aspect=paspect,cmap=cm,norm=norm, interpolation='nearest')
        plt.gca().yaxis.tick_right()
        plt.yticks(ticks=range(len(targets)),labels=targets)
        plt.xlabel('basepair')
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

        cb = plt.colorbar(im,orientation="horizontal") ### just to get cb.ax
        cb = plt.colorbar(im, format=fmt, ticks=tickz,
                              orientation="horizontal",
                             cax=cb.ax)
        
        
        #### metadata colorbar
        plt.sca(axs[0,0])

        #### get metadata labels sorted by number of samples
        labels =metadataDf.set_index('sampleName').loc[samples].metadata.value_counts().index.to_numpy() 
        metadataColorBar = sns.color_palette(n_colors=len(labels))
        np.sort(labels)
        col_dict = {}
        mapping_dict = {}
        for i,lab in enumerate(labels):
            col_dict[i] = metadataColorBar[i]
            mapping_dict[lab]=i
            
        toPlotMetadata = []
        for lab in metadataDf.set_index('sampleName').loc[samples].metadata:
            toPlotMetadata.append(mapping_dict[lab])
        
        len_lab = len(labels)
        metadata_cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

        norm_bins = np.sort([*col_dict.keys()]) + 0.5
        norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

        metadata_norm = mpl.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
        metadata_fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[metadata_norm(x)])
        diff = norm_bins[1:] - norm_bins[:-1]
        metadata_tickz = norm_bins[:-1] + diff / 2
        ####### end metadataColorbar
        
        im=plt.imshow(np.array(toPlotMetadata).reshape(1,-1), aspect='auto', cmap=metadata_cm,norm = metadata_norm,interpolation='nearest')
        plt.yticks([])
        plt.xlabel('Sample metadata')
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        cb = plt.colorbar(im,orientation="horizontal",pad=.5) ### just to get cb.ax
        cb = plt.colorbar(im, format=metadata_fmt, ticks=metadata_tickz,
                      orientation="horizontal",
                     cax=cb.ax)
        cb.ax.tick_params(labelrotation=-80)
        ######## end metadata plotting

        
        ### plot counts
        plt.sca(axs[0,1])
        im = plt.imshow(nj.reshape(1,-1),aspect='auto')
        plt.yticks([])
        plt.xlabel(r'Number of counts per column ($n_j$)')
        plt.colorbar(im,orientation="horizontal",pad=.5)
        
        ### title
        pv=np.nan
        esize=np.nan
        M=np.nan
        if anch in pvDf.index:
            anchRow = pvDf.loc[anch]
            pv = anchRow.pval_rand_init_alt_max_corrected
            esize=anchRow.effect_size_bin
            M=int(anchRow.M)
            
        fname ="{}_clust{}_partition-{}_pv{:.1E}_eSize{:.2F}_M{}_numSamples{}_{}.png".format(
            dsName,
             clustMapping[anch],
            perfectConcordance,
            pv, 
         esize,
            M,
            len(nj),
         anch
        )

        plt.suptitle(fname[:-4])
        plt.tight_layout()
        plt.savefig(outFolder+f'/plots/{fname}',bbox_inches='tight')
        plt.close()
    with open(outFolder+f'/runDetails.txt', "w") as text_file:
        print(statusStr, file=text_file)
    
        
        
        
### for cluster labelling
def clusterAnchors(anchLst):

    def is_partition(partition, n):
        # create a set of integers from 0 to n-1
        full_set = set(range(n))
        # create an empty set to store elements in the partition
        partition_set = set()
        for subset in partition:
            # convert each subset to a set
            subset_set = set(subset)
            # check if the subset is non-empty
            if not subset_set:
                return False
            # check if the subset is disjoint with the partition set
            if subset_set.intersection(partition_set):
                print('duplicate',subset_set)
                return False
            # add the subset to the partition set
            partition_set.update(subset_set)
        # check if the union of subsets is equal to the full set

        if partition_set != full_set:
            print('not complete',full_set-partition_set)

        return partition_set == full_set



    ### i,j-th entry is 1 if i-th row of A equals j-th row of B
    def compare_rows_hash(A, B):
        n, k = A.shape
        hash_A = np.array([hash(tuple(row)) for row in A])
        hash_B = np.array([hash(tuple(row)) for row in B])
        return hash_A[:, None] == hash_B



    #### main function
    def noiselessAssembly(anchLst,maxShiftDist=5):
        bpArr = np.array([[bpToInt(x) for x in s] for s in anchLst])
        ### generate similarity matrix, where i,j=1 indicates that
        #### read i falls ahead of read j within shift distance maxShiftDist
        n,k=bpArr.shape
        simMat = np.zeros((n,n),dtype=bool)
        for shift in range(1,maxShiftDist+1):
            simMatUpdate = compare_rows_hash(bpArr[:,shift:],bpArr[:,:-shift])
            simMat = np.logical_or(simMat,simMatUpdate)

        ### from this similarity matrix, generate clusters
        assemblies = clusterFromSimMat(simMat)

        assert is_partition(assemblies,n)


        return [[anchLst[i] for i in lst[::-1]] for lst in assemblies] ### output sequences
    
    
    def clusterFromSimMat(simMat):
        n = simMat.shape[0]

        sparseMat = scipy.sparse.lil_matrix(simMat)

        row_idx, col_idx, _ =find(sparseMat+sparseMat.T)


        followers = []

        for i in range(simMat.shape[0]):
            # Find the start index of the column indices for row i
            start = bisect_left(col_idx, i)
            # Find the end index of the column indices for row i
            end = bisect_right(col_idx, i)
            # Extract the column indices for row i
            followers.append(row_idx[start:end])



        ### list of finalized assemblies
        finalizedAssemblies = {}

        for idx in np.where([len(x)==0 for x in followers])[0]:
            finalizedAssemblies[idx]=[idx]

        isRoot = np.ones(n,dtype=bool)

        for idx in np.where([len(x)>0 for x in followers])[0]:

            if not isRoot[idx]:
                continue

            ### find the locations that follow this index        
            followerLst = findFollowers(followers,idx,finalizedAssemblies)

            ### all followers of this node cannot be a "Root"
            tmpVal = isRoot[idx]
            isRoot[idx] = False
            loopBreaker = np.any(isRoot[followerLst])
            isRoot[followerLst] = False
            isRoot[idx]=tmpVal or loopBreaker

            finalizedAssemblies[idx]=followerLst


        return [finalizedAssemblies[x] for x in np.where(isRoot)[0]]



    def findFollowers(followers,idx,finalizedAssemblies):
        retLst = set()
        to_visit = [idx]
        while to_visit:
            curr = to_visit.pop(0)
            if curr in finalizedAssemblies:
                retLst = retLst | set(finalizedAssemblies[curr])
                continue
            retLst.add(curr)
            for follower in followers[curr]:
                if follower not in retLst:
                    to_visit.append(follower)
        return list(retLst)
    
    
    outStr = ''
    assemblies = noiselessAssembly(anchLst, 2)
    clustMapping = {}
    for i,a in enumerate(assemblies):
        outStr+= f'Cluster {i}, length {len(a)} \n'
        for anch in a:
            outStr += f'{anch} \n'
            clustMapping[anch] = i
            
        outStr += '\n'
            
    ### can save fastaStr
    with open(outFolder+f'/plots/anch_lst.txt', "w") as text_file:
        print(outStr[:-2], file=text_file)
    
    return clustMapping

        
        
        
print('running')
args = get_args()
dsName = args.dsName
anchorFile = args.anchorFile
satcFolder = args.satcFolder
metadataPath = args.metadataPath
outFolder = args.outFolder
skipSATC = args.skipSATC
pvDfPath = args.pvDfPath
satc_dump_file = args.satc_dump_file

statusStr = "Dataset: {}\nAnchor file: {}\nSATC folder: {}\nmetadata: {}\npvDf: {}\nOut folder: {}\nskipSATC: {}".format(
dsName, anchorFile,satcFolder, metadataPath, pvDfPath,outFolder,skipSATC
)
print(statusStr)

if not skipSATC:
    parseSATC()
ctsDf = generateCtsDf()
generatePlots(ctsDf)

