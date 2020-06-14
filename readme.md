Riemannian kernel based nystrom method for approximate infinite-dimensional covariance descriptors (AidCovDs) with application to image set classification.

Written by Kai-Xuan Chen (e-mail: kaixuan_chen_jsh@163.com)
version 2.0 -- April/2019 
version 1.0 -- December/2017 

Please cite the following paper (more theoretical and technical details) if your are using this code:

Chen K X, Wu X J, Wang R, et al. Riemannian kernel based Nyström method for approximate infinite-dimensional covariance descriptors with application to image set classification[C]//2018 24th International conference on pattern recognition (ICPR). IEEE, 2018: 651-656.  

BibTex : 
```
@inproceedings{chen2018riemannian,
  title={Riemannian kernel based Nystr{\"o}m method for approximate infinite-dimensional covariance descriptors with application to image set classification},
  author={Chen, Kai-Xuan and Wu, Xiao-Jun and Wang, Rui and Kittler, Josef},
  booktitle={2018 24th International Conference on Pattern Recognition (ICPR)},
  pages={651--656},
  year={2018},
  organization={IEEE}
}
```


The ETH-80 dataset is needed to be downloaded(https://github.com/Kai-Xuan/ETH-80/),  and put 8 filefolders(visual image sets from 8 different categories) into filefolder '.\ETH-80\'.  

Please run 'read_ETH.m' to generate AidCovDs. Then run 'run_ETH.m' for image set classification.  

For Sift feature extraction, we employ VLFeat package.  Usage of VLFeat package, please cite the following paper:
Vedaldi A, Fulkerson B. VLFeat: An open and portable library of computer vision algorithms[C]//Proceedings of the 18th ACM international conference on Multimedia. ACM, 2010: 1469-1472.

For classification, we employ four NN classifiers and two discriminative classifiers in this source code(Version 2.0).  COV-LDA/COV-PLS:  This method was proposed in the paper:  
R. Wang, H. Guo, L. S. Davis, and Q. Dai. Covariance discriminative learning: A natural and efficient approach to image set classification. In Computer Vision and Pattern Recognition (CVPR), 2012 IEEE Conference on, pages 2496-2503. IEEE, 2012. 


For more experiment, you can test on Virus dataset (https://github.com/Kai-Xuan/Virus/)


