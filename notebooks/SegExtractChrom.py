from scipy import ndimage as nd
import scipy
import numpy as np
import readmagick
import mahotas
import pymorph
import pylab
import os
def distanceTranform(bIm):
    #from pythonvision.org
    dist = nd.distance_transform_edt(bIm)
    dist = dist.max() - dist
    dist -= dist.min()
    dist = dist/float(dist.ptp()) * 255
    dist = dist.astype(np.uint8)
    return dist
def BorderKill(imlabel):
    '''remove labelled objects touching the image border'''
    #from pythonvision.org
    whole = mahotas.segmentation.gvoronoi(imlabel)
    borders = np.zeros(imlabel.shape, np.bool)
    borders[ 0,:] = 1
    borders[-1,:] = 1
    borders[:, 0] = 1
    borders[:,-1] = 1
    at_border = np.unique(imlabel[borders])
    for obj in at_border:
        whole[whole == obj] = 0
    return whole
def gray12_to8(im):
    i=0.062271062*im
    return pymorph.to_uint8(i)
def gray16_to8(im):
    i=0.00390625*im
    return pymorph.to_uint8(i)

def GradBasedSegmentation(im):
    blur=nd.gaussian_filter(im, 16)
    rmax = pymorph.regmax(blur)
    T = mahotas.thresholding.otsu(blur)
    bImg0=im>T
    #bImg01=nd.binary_closing(bImg0,iterations=2)
    bImg01=pymorph.close(bImg0, pymorph.sedisk(3))
    bImg=pymorph.open(bImg01, pymorph.sedisk(4))
    #bImg=nd.binary_opening(bImg01,iterations=3)
    b=pymorph.edgeoff(bImg)
    d=distanceTranform(b)
    seeds,nr_nuclei = nd.label(rmax)
    lab=mahotas.cwatershed(d,seeds)
    return lab
def ModalValue(image):
    '''look for the modal value of an image'''
    #print image.dtype
    if image.dtype=="uint8":
        depthmax=255
        print "8bits"
    if image.dtype=="uint16":
        depthmax=65535
        print "16bits"
    histo=mahotas.fullhistogram(image)
    countmax=histo.max()
    print "countmax:",countmax
    print "image max",image.max()
    mig=image.min()#image min graylevel
    mag=image.max()#image max gray level
    mode=0
    countmax=0#occurence of a given grayscale
    print "mig=",mig,"  mag=",mag
    for i in range(mig,mag-1,1):
        test=histo[i]>countmax
        #print "test:",test,"histo(",i,")=", histo[i],"max",countmax
        if  test:
            countmax=histo[i]
            mode=i
            #print "mode",mode
    return mode
def RemoveModalBackground(image):
    mode=ModalValue(image)
    back = np.zeros(image.shape, image.dtype)
    back.fill(mode)
    #print "def background:",back.mean()
    im=pymorph.subm(image,back)
    return im
def LowResSegmentation(image):
    '''Perform a simple threshold after DoG filtering'''
    noBack=RemoveModalBackground(image)
    #print "from Segmeta noBack:",noBack.min(),noBack.mean()
    blurLowRes=nd.filters.gaussian_filter(noBack,13)
    blurHiRes=nd.filters.gaussian_filter(noBack,1)
    midPass=pymorph.subm(blurHiRes,0.70*blurLowRes)
    bin=(midPass>1.5*midPass.mean())
    binLowRes=pymorph.open(bin,pymorph.sedisk(4))
    return binLowRes
def extractParticles(grayIm,labIm):
    ''' give a grayscaled and a labelled image, extract the segmented particles
    ,returns a list of flattened particles'''
    #grayIm and labIm should have the same size
    def unflattenParticles(flatParticleList):
        '''take a list of flat particles and unflat them to yield an image'''
        unflatList=[]
        lenFlatList=len(flatParticleList)
        for i in range(0,lenFlatList):
            #get the i particle:current Particle
            curPart=flatParticleList[i]#current particle
            #x values(col) are stored in the third col (3-1)
            colmax=curPart[:,2].max()
            colmin=curPart[:,2].min()
            #y values(li) are stored in the fourth col (4-1)
            limax=curPart[:,3].max()
            limin=curPart[:,3].min()
            unflatIm=np.zeros((limax-limin+1,colmax-colmin+1),np.int16)
            #number of pixels in the particle
            nbPixel=len(curPart[:,1])#count how many lines at col=1
            for line in range(0,nbPixel):
                col=curPart[line,2]
                li=curPart[line,3]
                pixVal=curPart[line,1]
                unflatIm[li-limin,col-colmin]=pixVal
            unflatList.append(unflatIm)
        return unflatList
                   
    sx=grayIm.shape[0]
    sy=grayIm.shape[1]
    #flatten grayIm
    fg=grayIm.flatten()
    fl=labIm.flatten()
    labmax=fl.max()
    #print fg
    #print fl
    #build two 2D array containing x and y
    #of each pixel of the grayIm
    ax=np.zeros((sx,sy),np.int16)
    ay=np.zeros((sx,sy),np.int16)
    #vectorization with numpy may be 
    #more efficient than two loops
    for j in range(0,sy):
        for i in range(0,sx):
            ax[i,j]=j#filling ax with x=col
            ay[i,j]=i#filling ay with y values y=li
    #flat arrays of coordinates
    fax=ax.flatten()
    fay=ay.flatten()
    #1D merge graylevel, label and coordinates 
    #in one 1D array of 4-uplet
    extract=np.vstack((fl,fg,fax,fay))
    #transpose to watch it easily
    eT=extract.T
    #create a list of flatten particles
    #labIndex takes the value from 1 (the first particle to labmax the\
    #label of the last particle
    flatParticleList=[]#from Matthieu Brucher
    for labIndex in range(1,labmax+1):
        flatParticleList.append(eT[eT[:,0]==labIndex])#from Matthieu Brucher
    return unflattenParticles(flatParticleList)
#
#Modify your path to your images here. The script works with 16bits images
#
user=os.path.expanduser("~")
workdir=os.path.join(user,"Applications","ImagesTest","jp","Jpp48","2","DAPI")
file="1.tif"
complete_path=os.path.join(workdir,file)
#
#
if __name__ == "__main__":
    dapi=readmagick.readimg(complete_path)
    print "shape dapi",dapi.shape
    os.mkdir(os.path.join(workdir,"particules"))
    im1=RemoveModalBackground(dapi)
    #8bits dapi image
    d8=gray12_to8(im1)
    #try a simple segmentation procedure
    print "segmenting..."
    #imlabel=GradBasedSegmentation(im1)
    imlabel,npart=nd.label(LowResSegmentation(dapi))
    print "showing..."
    #pylab.imshow(imlabel)
    #pylab.show()
    #Particles is a list of images
    Particles=extractParticles(im1,imlabel)
    li=1
    ColNum=10#ten columns and len(Particles)/ColNum lines
    il=0
    iw=0
    for i in range(0,len(Particles)):
        #first convert and save particle as 8bits png image
        file='part'+str(i)+'.png'
        saveim=gray12_to8(Particles[i])
        pathtoparticles=os.path.join(workdir,"particules",file)
        savedfile=saveim.astype(np.uint8)
        readmagick.writeimg(savedfile, pathtoparticles)
        print "saving:"+file
        #loading a 8bits png
        #(next thing to do:building a mosaic directly from Particles list)
        #
        #print pathtoparticles
        #impng=readmagick.readimg(pathtoparticles)
        #col=i%ColNum+1
        #print "im shape:",impng.shape," i:",i," li:",li," col:",col
        #there should be len(Particles)+1 subplots
        #pylab.subplot(1,i+1,i+1)
        #pylab.imshow(impng)
        #if col == 1 :
        #    li=li+1
    #pylab.show()