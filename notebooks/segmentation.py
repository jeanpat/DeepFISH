# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 16:29:23 2014

@author: jeanpat
"""
import numpy as np
from scipy import ndimage as nd

import mahotas as mh


def extractParticles_2(greyIm, LabIm):
    #print 'start', greyIm.dtype, LabIm.dtype
    LabelImg= LabIm
    GreyImg = greyIm
    locations = nd.find_objects(LabelImg)
    
    #print locations
    i=1
    extracted_images=[]
    for loc in locations:
        
        lab_image = np.copy(LabelImg[loc])
        grey_image = np.copy(GreyImg[loc])
        
        lab_image[lab_image<>i]=0
        grey_image[lab_image <>i]=0
        extracted_images.append(grey_image)
        i=i+1
    return extracted_images
    
def branchedPoints(skel, showSE=True):
    X=[]
    #cross X
    X0 = np.array([[0, 1, 0], 
                   [1, 1, 1], 
                   [0, 1, 0]])
    X1 = np.array([[1, 0, 1], 
                   [0, 1, 0], 
                   [1, 0, 1]])
    X.append(X0)
    X.append(X1)
    #T like
    T=[]
    #T0 contains X0
    T0=np.array([[2, 1, 2], 
                 [1, 1, 1], 
                 [2, 2, 2]])
            
    T1=np.array([[1, 2, 1], 
                 [2, 1, 2],
                 [1, 2, 2]])  # contains X1
  
    T2=np.array([[2, 1, 2], 
                 [1, 1, 2],
                 [2, 1, 2]])
    
    T3=np.array([[1, 2, 2],
                 [2, 1, 2],
                 [1, 2, 1]])
    
    T4=np.array([[2, 2, 2],
                 [1, 1, 1],
                 [2, 1, 2]])
    
    T5=np.array([[2, 2, 1], 
                 [2, 1, 2],
                 [1, 2, 1]])
    
    T6=np.array([[2, 1, 2],
                 [2, 1, 1],
                 [2, 1, 2]])
    
    T7=np.array([[1, 2, 1],
                 [2, 1, 2],
                 [2, 2, 1]])
    T.append(T0)
    T.append(T1)
    T.append(T2)
    T.append(T3)
    T.append(T4)
    T.append(T5)
    T.append(T6)
    T.append(T7)
    #Y like
    Y=[]
    Y0=np.array([[1, 0, 1], 
                 [0, 1, 0], 
                 [2, 1, 2]])
    
    Y1=np.array([[0, 1, 0], 
                 [1, 1, 2], 
                 [0, 2, 1]])
    
    Y2=np.array([[1, 0, 2], 
                 [0, 1, 1], 
                 [1, 0, 2]])
    
    Y2=np.array([[1, 0, 2], 
                 [0, 1, 1], 
                 [1, 0, 2]])
    
    Y3=np.array([[0, 2, 1], 
                 [1, 1, 2], 
                 [0, 1, 0]])
    
    Y4=np.array([[2, 1, 2], 
                 [0, 1, 0], 
                 [1, 0, 1]])
    Y5=np.rot90(Y3)
    Y6 = np.rot90(Y4)
    Y7 = np.rot90(Y5)
    Y.append(Y0)
    Y.append(Y1)
    Y.append(Y2)
    Y.append(Y3)
    Y.append(Y4)
    Y.append(Y5)
    Y.append(Y6)
    Y.append(Y7)
    
    bp = np.zeros(skel.shape, dtype=int)
    for x in X:
        bp = bp + mh.morph.hitmiss(skel,x)
    for y in Y:
        bp = bp + mh.morph.hitmiss(skel,y)
    for t in T:
        bp = bp + mh.morph.hitmiss(skel,t)
        
    if showSE==True:
        fig = plt.figure(figsize=(4,5))
        tX =['X0','X1']
        tY =['Y'+str(i) for i in range(0,8)]
        tT =['T'+str(i) for i in range(0,8)]
        ti= tX+tY+tT
        SE=X+Y+T
        print len(SE), len(ti)
        n = 1
        ti = iter(ti)
        for se in SE:
            #print next(ti)
            #print se
            mycmap = mpl.colors.ListedColormap(['black','blue','red'])
            ax = fig.add_subplot(4,5,n,frameon=False, xticks=[], yticks=[])
            title(str(next(ti)))
            imshow(se, interpolation='nearest',vmin=0,vmax=2,cmap=mycmap)
            n = n+1
        fig.subplots_adjust(hspace=0.1,wspace=0.08)
        #ax_cb = fig.add_axes([.9,.25,.1,.3])#
        color_vals=[0,1,2]
        #cb = mpl.colorbar.ColorbarBase(ax_cb,cmap=mycmap, ticks=color_vals)
        #cb.set_ticklabels(['back', 'hit', 'don\'t care'])
        
        plt.show()
    return bp
    
def chromatids_elements (TopHatedChromosome):
    '''Take a High pass filtered (or top hat) image of a chromosome and label the chromatids elements
    '''
    threshed = TopHatedChromosome > 0
    #threshed = mh.open(threshed)
    labthres, _ =mh.label(threshed)
    labsz = mh.labeled.labeled_size(labthres)
    mh.labeled.remove_regions_where(labthres,labsz<2, inplace=True)
    threshed = labthres>0
    
    skel2 = mh.thin(threshed)
    bp2 = branchedPoints(skel2, showSE=False)>0
    rem = np.logical_and(skel2,np.logical_not(bp2))
    labskel, _ = mh.labeled.label(rem)
    #print labskel.dtype
    size_sk= mh.labeled.labeled_size(labskel)
    #print size_sk
    skelem = mh.labeled.remove_regions_where(labskel, size_sk < 4)
    
    distances = mh.stretch(mh.distance(threshed))
    surface = (distances.max() - distances)
    chr_label = mh.cwatershed(surface, skelem)
    #print chr_label.dtype, type(chr_label)
    chr_label *= threshed
    
    #This convertion is important !!
    chr_label = chr_label.astype(np.intc)
    #-------------------------------
    mh.labeled.relabel(chr_label, inplace=True)
    labsize2 = mh.labeled.labeled_size(chr_label)
    cleaned = mh.labeled.remove_regions_where(chr_label, labsize2 < 8)
    mh.labeled.relabel(cleaned, inplace=True)
    return cleaned