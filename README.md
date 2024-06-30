# Stuff
stuff stuff stuff

# Notes for installation into conda environment
In order to develop stuff in R on servers without root permissions and without R-studio, I employ **conda**. It's easy to create a R environments and to install R packages in it. As easy as  
`conda install r-<package-name> -c conda-forge`  
It manages the environment as it does with python. Nonetheless, sometimes you have to install stuff that's not present in CRAN.  
Maybe you want to install such packages into your environment. You can see the library of the environment by activating it and then running in R simply  
`.libPaths()`  
When you activate conda environments, R's libpaths change. At this point, for installing something in it, I would proceed from source.  
First clone the github repo in local or download its zip file. Then, activate the conda environment and from R run the following:  
`install.packages(<local_path_of_source>, repo=NONE, type='source', lib=.libPaths())`  