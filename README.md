# Multi-sensor labelled multi-Bernoulli filters

This repository contains Matlab implementations of various labelled multi-Bernoulli (LMB) and LMB mixture (LMBM) filters.
The repository contains both single- and multi-sensor implementations of the LMB and LMBM filters. 
These filters can be implemented using a variety of data association algorithms; however, all filters assume the linear-Gaussian dynamics.
Hopefully, the code is documented well enough that is easy to interpret. 
This code runs in Matlab R2022a and requires the **Statistics and Machine Learning Toolbox**, but only to simulate Poisson distributed random variables.

## Single-sensor LMB and LMBM filters

The script **runFilters.m** runs the single-sensor LMB and LMBM filters and plots their results and Euclidean and Hellinger optimal subpattern assignment (OSPA) errors.
The LMB filter can be run using the following three data association algorithms:

   1. Loopy belief propagation (LBP). This Williams et al.'s LBP algorithm to approximate each object's posterior existence probability and marginal association probabilities. We recommend this data association algorithm, as it computationally inepxensive and it is more accurate than the other two data association algorithms.
   2. Gibbs sampling. This uses a relatively inexpensive Gibbs sampling routine to approximate each object's posterior existence probability and marginal association probabilities.
   3. Murty's algorithm. This uses Vo and Vo's **.mex** implementation of Murty's algorithm to approximate each object's posterior existence probability and marginal association probabilities. This code may have a memory leak.

All of the single-sensor LMB filters approximate each object's spatial distribution using a Gaussian mixture (GM), and the filters' parameters can be set in the script **generateModel.m**.
The LMBM filter can be implemented using both the Gibbs sampler and Murty's algorithm.
However, it cannot be implemented using the LBP algorithm, as that algorithm cannot be used to generate data association events and can only approximate marginal distributions.

## Multi-sensor LMB and LMBM filters

The script **runMultisensorFilters.m** runs the multi-sensor LMB and LMBM filters and plots their results and Euclidean and Hellinger optimal subpattern assignment (OSPA) errors.
We have developed the following three approximate mutli-sensor LMB filters:

  1. The parallel update LBM (PU-LMB) filter. This filter results from the mathematical manipulation of the multi-sensor multi-object Bayes filter's posterior distribution. This filter assumes independent sensors, and it is the most accurate of our approximate multi-sensor LMB filters.
  2. The geometric average LMB (GA-LMB) filter. This filter is based of geometric average (GA) fusion, and it approximates the multi-sensor multi-object Bayes filter's posterior distribution using the weighted GA of each sensor's measurement-updated distribution. This filter is accurate in terms of object localisation (error in the objects' state), but provides a poor cardinality (number of objects) estimate. This filter does not assume the sensors are independent, and it provides a poor covariance estimate for independent sensors.
  3. The arithmetic average LMB (AA-LMB) filter. This filter is based of arithmetic average (AA) fusion, and it approximates the multi-sensor multi-object Bayes filter's posterior distribution using the weighted AA of each sensor's measurement-updated distribution. This filter does not assume the sensors are independent, and it provides the worst results of all our filters for independent sensors. However, its cardinality estimate is superior to the GA-LMB filter's.

The three multi-sensor LMB filters listed above all split the multi-sensor measurement update into independent single-sensor measurement updates, before combining the resulting measurement-updated distributions together.
Their measurement updates can be computed in parallel and, when combined with the LBP data association algorithm, their computational complexities are constant in the number of sensors, and linear in the number of objects and meausurements.
