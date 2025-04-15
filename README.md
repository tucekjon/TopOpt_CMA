# Topology optimization with Characteristic Modes Analysis
Implementation of the density-based topology optimization in method of moments modelling in conjuction with Characteristic mode analysis. The first characteristic number of a rectangular plate is manipulated

## Implementation notes
The density topology optimization is fully implemented in MATLAB according to [1]. A modified version of Method of Moving Asymptotes is utilized to update design variable and Adjoint sensitivity analysis provides the sensitivities of the objectives. The standard filtering technique (density and projection filters) is performed to regularize the solution space and accelerate the optimization. The objective is to make the first significant characteristic mode resonant, i.e., the goal is to minimize the magnitude of the first characteristic mode towards zero, see [2] for more details.

## Example
The example utilizes pre-calculated data from method-of-moments simulation to reduce both the code complexity and the computational cost of the demonstration. The code is fully compatible with the outputs of AToM package [3]. Pre-calculated data are provided for perfectly conducting rectangular region. The plate is discretized into 800 triangles and covered with 1170 basis functions. 

The objective is to minimize fitness function $f=\mathrm{Re}(\lambda_1)^2 + \nu \mathrm{Im}(\lambda_1)^2$, where $\lambda_1$ is first significant mode obtained from characteristic mode analysis and $\nu$ is a heuristic parameter which punishes imaginary part of $\lambda_1$, i.e., it damps losses in the design and has a great impact on the optimization, see [2].

## Initiation and start
Optimization parameters can be set at the beginning of START.m script, which serves as a starting script and runs automatically after pressing "F5". No extra code is required. In particular parameter $\nu$ has a great impact on the optimization and performance after postprocessing, see the convergence plot and obtained designs below.

<p align="center">
  <img src="https://github.com/tucekjon/TopOpt_CMA/blob/main/TopOpt-CMA-Convergence.png?raw=true" width="450" />
</p>
<em>An example of the convergence plot of fitness fuction, where dashed vertical lines represent iterations in which the sharpness of the projection filter is doubled.</em>

<p align="center">
  <img src="https://github.com/tucekjon/TopOpt_CMA/blob/main/TopOpt-CMA-OptimizedDesign.png?raw=true" width="400" />
  <img src="https://github.com/tucekjon/TopOpt_CMA/blob/main/TopOpt-CMA-BinaryDesign.png?raw=true" width="400" /> 
</p>
<em>An example of the  optimized structure with residual gray elements (left). An example of the thresholded binary structure (right).</em>


## References
[1] Tucek, J.,Capek, M., Jelinek, L., Sigmund, O.: Q-factor Minimization via Topology Optimization, 
     pp. 1-13, 2023.

[2] Tucek, J.,Capek, M., Jelinek, L.: Density-Based Topology Optimization for Characteristic Modes Manipulation, 
     pp. 1-13, 2023.
     
[3] Antenna Toolbox for MATLAB (AToM), [on-line]: www.antennatoolbox.com, (2022)
