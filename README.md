# symplastic-guard-cell-connections

MorophDynamX and MorphoMechanX code used in "Symplastic guard cell connections buffer pressure fluctuations to promote stomatal function in grasses" Wilson et. al. 2024. DOI:

The simulations were run on an Ubuntu 22.04 machine, with an NVIDIA graphics card. CUDA version: 12.3.107, and NVIDIA Driver Version: 560.35.03 

To install MorphoDynamX/MorphoMechanX for Ubuntu 22.04.
Install MorphoDynamX using the provided .deb package

```$ sudo dpkg -i MDX-2.0.2-1582-Ubuntu22.04-Cuda12-CellMaker-DivisionAnalysis-FemLib-Gmsh.deb```

If MorphoDynamX fails to install, the terminal should give you an indication of why. This is most likely due to a dependency issue. Some known dependencies that could be missing at this stage are:

g++
This can be installed using apt-get. Example,

```$ sudo apt-get install g++```

## Barley simulations

With MorphoDynamX installed, we now need to compile the source code for the models to run using MorphoMechanX.

- In the terminal, navigate to the directory 'code'. Then run the following,
- ```$ make clean```
- ```$ make```
- ```$ make run```

If an error occurs, most likely there are missing dependencies. Once again, the terminal should give an indication of what is missing. When setting up the Virtual Machine the following dependencies were missing and were installed:

- CUDA/Thrust (see above for installation)
- libqt5opengl5-dev
- libboost-all-dev
- libgsl-dev
- libtbb-dev
- libqtmultimedia5
- libqtmultimediawidgets5
- libqtmultimedia5-plugins
- qtmultimedia5-dev

For a full list of potential dependencies, please see the file dependencies.txt included in this repository.

The following Python packages need to be installed:
1. numpy
2. (optional) PIL (for making snapshots and gifs of the pressurized stoma meshes)
These can be installed using pip.

Example: ```$ pip install numpy```

To run the automated grass models:
Update globals.py to set your target output directory and change any model parameters. It is necessary to change the output directory to follow your own file structure.
- ```$ make clean```
- ```$ make```
- ```$ make run and MorphoDynamX should open.```
In MorphoDynamX, under /Process/Tools/Python/Python Script input inflate-outputs.py.
Press 'Step' (the play button) to run the simulation. This should now run a simulation.
To create pictures/gif of the simulation.

In MorphoDynamX, under /Process/Tools/Python/Python Script input take-pics.py.
Press 'Step' (the play button) to create the pictures/gif. The output input/output directory is that specified in globals.py. (Also see limitations below)

## Onion models

To run the automated onion models:
Run the run_mdx_onion.sh script - ```$ ./run_mdx_onion.sh```. There is no need to compile the code first, as this script will do this automatically. Select "Symplastic Connections/Figure simulations" from the drop down menu on the right. Figures 3b, 3c, S2a and S2b can be run by selecting the correct respective label in the "Figure" field.

## Current limitations

- Only one Python script can run with each launch of MorphoDynamX. Once one Python script has run and completed, MorphoDynamX needs to be closed and re-launched using the above steps.
- To create pictures/gifs:
  - The Python script take-pics.py works by taking screenshots within MorphoDynamX. Therefore, the overlayed text will not be in the same location for each setup. The user will need to edit take-pics.py to overlay text for their images.
  - Similarly, the heatmap is overlayed as well. If the heatmap for stress is changed from 0-75 MPa, heatmap.png will need to be replaced with a picture of the new legend.
 
## Microscopy data

For the images used in this study, please see: [microscopy data](https://doi.org/10.5281/zenodo.14793472).
 
## Help
For additional resources on MorphoGraphX, MorphoDynamx and MorphoMechanX please visit https://morphographx.org and Strauss, Sören, et al. "Using positional information to provide context for biological image analysis with MorphoGraphX 2.0." Elife 11 (2022): e72601.

## Citation

If you use this repository, please cite the following work:

Wilson, M.J, McGregor, S., Durney, C. H., Tomkins, M., Armand, J., Smith, R. S., Gray, J. E., Morris, R. J., Fleming, A. J. Symplastic guard cell connections buffer pressure fluctuations to promote stomatal function in grasses.

For any questions please contact Melissa at melissa.tomkins@jic.ac.uk.
