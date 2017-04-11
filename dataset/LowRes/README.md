This directory (*DeepFISH/dataset/LowRes*) contains:

   * **train/grey and train/groundtruth**: 100085 pairs of 8bits images (grey+groundtruth) at low resolution (16 times lower than the original raw tiff images) in png format. Those pairs can be use to train a semantix segmentation network
   
   * **validation/grey and validation/groundtruth**: contains 111123 pairs (grey + groundtruth label) of images (8bits, png) with the same resolution than those in the training dataset.


The training and the validation datasets were generated with different chromosomes belonging to the same metaphase (*DeepFISH/Raw images/jpp21/3*) from human normal lymphocytes couterstained with DAPI and hybridized with a telomeric Cy3-labelled PNA probe.
