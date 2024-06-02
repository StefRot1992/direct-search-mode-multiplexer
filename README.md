# Direct Search Optimisation Algorithm for Mode Multiplexers
This repository accompanies the paper entitled"Unlocking Mode Programming with Multi-Plane Light Conversion Using Computer-Generated Hologram Optimisation" (S. Rothe, F. Barbosa, J. W. Czarske, F. M. Ferreira)

With this repository, we want to make our method for the creation of specialised phase masks for mode multiplexers that use multiple spatial modes publicly available and thus reproducible. 

We have created and uploaded an example for use in Matlab. 

# File List
simulation_direct_search_MAIN.m --- this is the main of our code example, from which a calculation process can be started. For successful execution, the other files must also be within the same folder path. Details can be found in the paper and we have tried as far as possible to make intermediate steps logically comprehensible with comments.

propagate.m --- this is one fucntion we have separated. It executes free-space propagation based on a kernel. This function was taken from the repository "Laguerre-Gaussian mode sorter Supplemental Information" (Fontaine, Nicolas K., Ryf, Roland, Chen, Haoshuo, Neilson, David T., Kim, Kwangwoong, and Carpenter, Joel) (2019). Data Collection. https://doi.org/10.14264/uql.2019.81. 

WMA_5modes_5passages.mat --- this is a file that is loaded automatically to your workspace when running our code. It provides all necessary parameters that must be set to reproduce our method. 

workspace_DS_5passages_5targetModes_100k_iterations.mat --- this is one result on a 5-mode mode multiplexer that we provide here. It is the result of exactly the code in the repository, when running 100k iterations. In our main file, you can display these results when setting flag_produce_new_results=0 and switch flag_plot=1.
