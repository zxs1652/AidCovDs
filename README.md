Riemannian kernel based nystrom method for approximate infinite-dimensional covariance descriptors (AidCovDs) with application to image set classification.

Written by Kai-Xuan Chen (e-mail: kaixuan_chen_jsh@163.com)
version 2.0 -- April/2019 
version 1.0 -- December/2017 

Please cite the following paper (more theoretical and technical details) if your are using this code:

Kai-Xuan Chen, Xiao-Jun Wu, Rui Wang, Josef Kittler. Riemannian kernel based nystrom method for approximate infinite-dimensional covariance  descriptors with application to image set classification. In 2018 24th International Conference on Pattern Recognition (ICPR), pages 651–656. IEEE, 2018.

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


For classification, we employ two classifiers in this source code(Version 2.0).  

COV-LDA/COV-PLS :  This method was proposed in the paper:  
R. Wang, H. Guo, L. S. Davis, and Q. Dai. Covariance discriminative learning: A natural and efficient approach to image set classification. In Computer Vision and Pattern Recognition (CVPR), 2012 IEEE Conference on, pages 2496-2503. IEEE, 2012. 