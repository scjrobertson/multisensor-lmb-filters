# Multi-sensor labelled multi-Bernoulli filters

This repository contains Matlab implementations of various labelled multi-Bernoulli (LMB) and LMB mixture (LMBM) filters.
The repository contains both single- and multi-sensor implementations of the LMB and LMBM filters. 
These filters can be implemented using a variety of data association algorithms; however, all filters assume the linear-Gaussian dynamics.
The code is hopefully documented well enough that it is easy to interpret. 
You should be able to use Matlab's **help** function to access the scripts' documentation.
This code runs in Matlab R2022a and requires the **Statistics and Machine Learning Toolbox**, but only to simulate Poisson distributed random variables.

## Single-sensor LMB and LMBM filters

The script **runFilters.m** runs the single-sensor LMB and LMBM filters and plots their results and Euclidean and Hellinger optimal subpattern assignment (OSPA) metrics.
The LMB filter can be run using the following three data association algorithms:

   1. Loopy belief propagation (LBP). This is Williams et al.'s LBP algorithm that approximates each object's posterior existence probability and marginal association probabilities. We recommend this data association algorithm, as it is computationally inexpensive and it is more accurate than the other two data association algorithms.
   2. Gibbs sampling. This uses a relatively inexpensive Gibbs sampling routine to approximate each object's posterior existence probability and marginal association probabilities.
   3. Murty's algorithm. This uses Vo and Vo's **.mex** implementation of Murty's algorithm to approximate each object's posterior existence probability and marginal association probabilities.

All of the single-sensor LMB filters approximate each object's spatial distribution using a Gaussian mixture (GM), and the filters' parameters can be set in the script **common/generateModel.m**.
The LMBM filter can be implemented using both the Gibbs sampler and Murty's algorithm.
However, it cannot be implemented using the LBP algorithm, as that algorithm cannot be used to generate data association events and can only approximate marginal distributions.

## Multi-sensor LMB and LMBM filters

The script **runMultisensorFilters.m** runs the multi-sensor LMB and LMBM filters and plots their results and Euclidean and Hellinger OSPA metrics.
We have developed the following three approximate multi-sensor LMB filters:

  1. The parallel update LBM (PU-LMB) filter. This filter results from the mathematical manipulation of the multi-sensor multi-object Bayes filter's posterior distribution. This filter assumes that the sensors are independent, and it is the most accurate of our approximate multi-sensor LMB filters.
  2. The geometric average LMB (GA-LMB) filter. This filter is based on geometric average (GA) fusion, and it approximates the multi-sensor multi-object Bayes filter's posterior distribution using the weighted GA of each sensor's measurement-updated distribution. This filter is accurate in terms of object localisation (i.e. the objects' kinematic states), but provides a poor instantaneous cardinality (number of objects) estimate. However, it usually does track all the objects present. This filter does not assume the sensors are independent, and it provides a poor covariance estimate for independent sensors.
  3. The arithmetic average LMB (AA-LMB) filter. This filter is based on arithmetic average (AA) fusion, and it approximates the multi-sensor multi-object Bayes filter's posterior distribution using the weighted AA of each sensor's measurement-updated distribution. This filter does not assume the sensors are independent, and it provides the worst results of all our filters for independent sensors. However, its cardinality estimate is superior to the GA-LMB filter's.

The three multi-sensor LMB filters listed above all split the multi-sensor measurement update into independent single-sensor measurement updates, before combining the resulting measurement-updated distributions together.
Their measurement updates can be computed in parallel significantly reducing their computational cost.
However, in these implementations, the filters' measurements updates are not computed in parallel, and the **Parallel Computing Toolbox** is not required.
Only the AA-LMB filter propagates GMs, the PU- and GA-LMB filters assume an object's prior spatial distribution is Gaussian and approximates each object's posterior spatial distribution as Gaussian.
It is also possible to implement all three filters using a Gibbs sampler or Murty's algorithm.

We have also implemented an iterated-corrector LMB (IC-LMB) filter that computes each sensor's measurement update in sucession.
The IC-LMB filter propagates a GM for each object's spatial distribution, and it represents a typical implementation of a multi-sensor LMB filter.
We have also implemented a multi-sensor LMBM filter using a variant of our Gibbs sampler.
This filter represents an exact closed-form solution to the multi-sensor multi-object Bayes filter's recursion.
It accounts for a variable number of sensors; however, it is prohibitively expensive.
If you track a large number of objects using many sensors, then you will exceed Matlab's memory limit.

All of the multi-sensor filters' parameters are set in the script **common/generateMultisensorModel.m**.

## Other things

We also have some additional scripts that allow us to compare and contrast both our data association algorithms and filters.
The scripts are organised into the following two folders:

   1. **marginalEvaluations/:** The scripts in this folder compare the LBP data association's approximate marginal distrubtions to those produced by Murty's algorithm and our Gibbs sampler. 
The Gibbs sampler is based on the same underlying model as the LBP algorithm.
   2. **trials/:** The scripts in these folders compare the various single- and multi-sensor filters' OSPA metrics and runtimes in various scenarios. 
