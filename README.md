# Data point cloud separability analysis based on linear Fisher discriminants

This github implements the methodology of separability analysis of data point clouds described in [Gorban A.N. et al, Correction of AI systems by linear discriminants: probabilistic foundations. 2018. Information Sciences 466:303-322.](https://www.sciencedirect.com/science/article/pii/S0020025518305607) The actual method is described in [L.Albergante, J.Bac, Zinovyev A. Estimating the effective dimension of large
biological datasets using Fisher separability analysis. Proceedings of International Joint Conference on Neural Networks-2019, Budapest, Hungary](https://arxiv.org/pdf/1901.06328.pdf).

The method implemented allows a) quantifying the effective data dimension based on comparing the separability distributions with a uniformly sampled n-dimensional sphere, and b) quantify the empirical probability distribution of a data point to be separated from the rest of the data point cloud.

## MATLAB implementation (provided by [A.Zinovyev](https://github.com/auranic))

In order to test the code, execute

	addpath tests;
	testSeparabilityAnalysis;

## Python 3 implementation (provided by [J.Bac](https://github.com/j-bac) and [L.Albergante](https://github.com/albluca))

In order to test and use the code, open the python notebook

## Acknowledgements

Supported by the University of Leicester (UK), the Ministry of Education and Science of the Russian Federation, project N 14.Y26.31.0022
