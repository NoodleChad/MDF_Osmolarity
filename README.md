# Understanding systems level metabolic adaptation resulting from osmotic stress
In this repository, you will every code and ressources used to reproduced the data presented in our paper.\
To proceed, copy the repository and install the conda environment.yml file.\
To reproduce the figure, simply select the python file in MDF_Osmolarity repository with the name of the figure you want to recreate. These will perform the analysis used to generate the simulated data.\
The rest of the files founds on the main repository contain all of the required functions.

## Folders description
**Experimental** Data generated in the wet lab\
**Results** Results of each of the analysis\
**To_make_model** These codes should allow you to recreate your own thermodynamic models which are going to be the input of all further analysis. We used iML1515 directly retrieved from the BiGG database\
**Tools_Figure_3** Contains the code used to directly take the data of what is feasible or not (Figure 3) and transforms it in way that makes it amenable for Escher maps.
**Uncertainty** The results of our Delta G uncertainty analysis\
**resources** This file contains everything that was required as an input for our analysis\
Question? alexandre.tremblay@mail.utoronto.ca

## Reproduce the results
### Dependencies
1. Annaconda ([here](https://www.anaconda.com/)) to set-up the conda environement
2. IBM CPLEX solver version >= 12.10
### Steps
1. Clone repository
2. Download and install the Python environment (environment.yml) in your desired folder using the following line:
```sh
conda env create -n osmdf -f environment.yml
```
3. Change the location of saved data.
4. Run the name of the figure .py








