# DeepFISH
Why Deep? because of deep-learning. Why FISH? Because of Fluorescent in-situ Hybridization.

This repository is intented to share data and code for resolving some problems met in cytogenetics imaging such overlapping chromosomes.

### dataset:
   * **DAPI.tif and Cy3.tif** : 12 bits images of metaphasic chromosomes. The telomeres marking the end of the chromosomes are visible in the Cy3.tif image.The metaphase doesn't contain overlapping chromosomes.

   * **LowRes_13434_overlapping_pairs.h5** : 13434 pairs of overlapping chromosomes generated from the two previous images. This dataset is intented to train a supervised learning algorithm to resolve overlapping chromosomes. The dataset is stored as a numpy array and saved in a hdf5 file. Compared to the DAPI and Cy3 images,the resolution was decreased by two.
   * **overlapping_chromosomes_examples.h5**: smaller dataset (~2000 images). The resolution of the images is the same than the DAPI/Cy3 images. 
 