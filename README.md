# DeepFISH
Why Deep? Because of deep-learning.

Why FISH? Because of **F**luorescent **i**n-**s**itu **H**ybridization.

This repository is intended to share data and code for resolving some problems met in cytogenetics imaging such overlapping chromosomes.

### Dataset
   * **DAPI.tif and Cy3.tif** : 12 bits images of metaphasic chromosomes. The telomeres marking the end of the chromosomes are visible in the Cy3.tif image.The metaphase doesn't contain overlapping chromosomes.

   * **LowRes_13434_overlapping_pairs.h5** : 13434 pairs of overlapping chromosomes generated from the two previous images. This dataset is intended to train a supervised learning algorithm to resolve overlapping chromosomes. The dataset is stored as a numpy array and saved in a hdf5 file. Compared to the DAPI and Cy3 images,the resolution was decreased by two.
   * **overlapping_chromosomes_examples.h5**: smaller dataset (~2000 images). The resolution of the images is the same than the DAPI/Cy3 images.

### References
  * [A Geometric Approach To Fully Automatic Chromosome Segmentation](https://arxiv.org/abs/1112.4164)
  * [Automated Discrimination of Dicentric and Monocentric Chromosomes by Machine Learning-based Image Processing](http://biorxiv.org/content/biorxiv/early/2016/01/19/037309.full.pdf)
  * [An Efficient Segmentation Method for Overlapping Chromosome Images](http://research.ijcaonline.org/volume95/number1/pxc3894861.pdf)
  * [A Review of Cytogenetics and its Automation](http://www.scialert.net/qredirect.php?doi=jms.2007.1.18&linkid=pdf)
