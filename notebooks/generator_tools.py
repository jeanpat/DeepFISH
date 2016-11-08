import itertools

from scipy import ndimage as nd
import numpy as np
import skimage as sk
print sk.__version__
from skimage import io
from skimage import segmentation as skg

def ResizeImages(ImList):
        '''Find the largest width and height of images belonging to a list.
        Return a list of images of same width/height
        '''
        maxwidth=0
        maxheight=0
        if len(np.shape(ImList[0]))==3:
            components = np.shape(ImList[0])[2]
        imtype = ImList[0].dtype
        for i in range(len(ImList)):
            width=np.shape(ImList[i])[1]#width=column
            height=np.shape(ImList[i])[0]#height=line
            #print "width:height",width,":",height
            if width>maxwidth:maxwidth=width
            if height>maxheight:maxheight=height
        #print "maxwidth:maxheight",maxwidth,":",maxheight
        NewList=[]
        for i in range(0,len(ImList)):
            width=np.shape(ImList[i])[1]
            height=np.shape(ImList[i])[0]

            diffw=maxwidth-width
            startw=round(diffw/2)
            diffh=maxheight-height
            starth=int(round(diffh/2))
            startw=int(round(diffw/2))
            if len(np.shape(ImList[0]))==3:
                newIm=np.zeros((maxheight,maxwidth,components), dtype=imtype)
                newIm[starth:starth+height,startw:startw+width,:]=ImList[i][:,:,:]
                NewList.append(newIm)
            if len(np.shape(ImList[0]))==2:
                newIm=np.zeros((maxheight,maxwidth), dtype=imtype)
                newIm[starth:starth+height,startw:startw+width]=ImList[i][:,:]
                NewList.append(newIm)
        return NewList
def clip_img_to_bounding_box(img):
    bb = nd.find_objects(img[:,:,-1]>0)
    slice0 = bb[0][0]
    slice1= bb[0][1]
    clip = img[slice0,slice1]
    return clip

def patch_to_square(image, seepatch=False):
    if seepatch==True:
        s1=2
        s2=3
    else:
        s1=0
        s2=0
    row = image.shape[0]
    col = image.shape[1]
    Hyp = int(np.ceil(np.sqrt(row**2+col**2)))+1
    drow = int(np.ceil(Hyp-row)/2)
    dcol = int(np.ceil(Hyp-col)/2)
    patch_h= s1*np.ones((row,dcol),dtype=int)
    patch_v= s2*np.ones((drow,col+2*(dcol)), dtype=int)
    e2 = np.hstack((patch_h, image, patch_h))
    return np.vstack((patch_v,e2,patch_v))

def rotated_images(image, step_angle, half_turn = True):
    if half_turn == True:
        angle_max = 180
    else:
        angle_max = 360
    angles = np.arange(0, angle_max, step_angle)
    return [clip_img_to_bounding_box(nd.rotate(image, rotation)) for rotation in angles]

def collection_of_pairs_of_rotated_images(image1, image2, step_angle=10, half_turn=True):
    '''Take two images, rotate them all by "step_angle".
    Make all possible pairs by cartesian product
    then resize them such same shape(size)
    '''

    r_images1 = rotated_images(image1, step_angle, half_turn)
    r_images2 = rotated_images(image2, step_angle, half_turn)

    pairs = itertools.product(r_images1, r_images2)
    #print type(next(pairs))
    #print next(pairs)
    r_pairs = [ResizeImages([p[0],p[1]]) for p in pairs]
    return r_pairs

def translate(image, vector):
    #print 'vector',vector, vector[0]
    image = np.roll(image, int(vector[0]), axis = 0)
    image = np.roll(image, int(vector[1]), axis = 1)
    return image

def add_mask_to_image(image, mask_value = 1):
    image = np.dstack((image, mask_value*(image>0)))
    return image

def merge_and_roll(still_img, move_img, row, col, clip=True):
    '''row, col: assume two numbers in [0,1]
    images:last component supposed to be a mask
    '''
    u=row
    v=col

    target = np.copy(still_img)
    source = np.copy(move_img)
    #print target.shape, source.shape

    #bounding boxes
    bb0 = nd.find_objects(target[:,:,-1]>0)
    bb1 = nd.find_objects(source[:,:,-1]>0)
    #won't work if more than two components C0 and C1
    C_0 = target[bb0[0][0],bb0[0][1]]
    C_1 = source[bb1[0][0],bb1[0][1]]
    #col, row
    c1 = C_1.shape[1]
    r1 = C_1.shape[0]
    c0 = C_0.shape[1]
    r0 = C_0.shape[0]
    comp = C_0.shape[-1]
    #print 'c1,r1,c0,r0, comp:',c1,r1,c0,r0, comp

    still = np.zeros((r1+2*r0,c1+2*c0,comp),dtype=int)
    still[r0+1:r0+r1+1,c0+1:c0+c1+1]= C_1
    move = np.zeros(still.shape, dtype=int)
    move[:r0,:c0]=C_0
    vector = u*(still.shape[0]-r0),v*(still.shape[1]-c0)
    #print r1, c1, vector
    move = translate(move, vector)
    #merge_shape =
    #figsize(10,8)
    merge = move+still
    if clip==False:
        return merge
    else:
        return clip_img_to_bounding_box(merge)

    
def overlapping_generator(image1, image2, mask1 = 1, mask2 = 2, rotation_step = 30, translations_number = 9):
    '''
    * This function takes two greyscaled images of chromosomes on black backaground
    * add a mask to each image:
        + same value for each mask (1 and 1) so that the overlapping pixels are set to 2.
        + different values for each mask (1 and 2) so that the overlapping pixels are set to 3.
    * rotate one chromosome of a multiple of a fixed angle in degree and make a pair of chromosomes
    (rotated, fixed position).
    *perform relative translations of one chromosome with a given number
    *returns a list of merged pair of images.
    '''
    c0 = clip_img_to_bounding_box(add_mask_to_image(image1, mask_value = mask1))
    c1 = clip_img_to_bounding_box(add_mask_to_image(image2,mask_value = mask2))
    
    CR = collection_of_pairs_of_rotated_images(c0, c1, step_angle = rotation_step)
    
    #Prepare translation vectors u (same on rows and lines)
    u = np.linspace(0,1,num = translations_number, endpoint = True)[1:-1]
    #v = np.arange(0,1,0.2)[1:]
    P = [t for t in itertools.product(u,u)]
    
    overlappings = []
    #Take a pair of images
    for pair_of_images in CR:# CR[::3] Just play with a subset
        im1 = pair_of_images[0]
        im2 = pair_of_images[1]
        #Now translate one image relative to the other one
        for t in P:#P[::3]Just play with a subset
            u = t[0]
            v = t[1]
            #translations = [merge_and_roll(im1, im2, u, v) for t in P]
            #overlappings.append(translations)
            overlappings.append(merge_and_roll(im1, im2, u, v))

    overlapping_chrom_learning_set = ResizeImages(overlappings)
    return overlapping_chrom_learning_set
