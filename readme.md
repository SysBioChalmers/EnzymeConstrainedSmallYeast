# Enzyme Constrained Small Yeast (Metabolic Trade-offs in Yeast are Caused by F1FO-ATP synthase)

- Models and scripts for simulations of ATP and biomass synthesis in yeast, the repo contains a small scale enzyme constrained stoichiometric model of intermediary metabolism and scripts to reproduce the figures.

- Abstract:
Intermediary metabolism provides living cells with free energy and precursor metabolites required for
synthesizing proteins, lipids, RNA and other cellular constituents, and it is highly conserved among
living species. Only a fraction of cellular protein can, however, be allocated to enzymes of intermediary
metabolism and consequently metabolic trade-offs may take place. One such trade-off, aerobic
fermentation, occurs in both yeast (the Crabtree effect) and cancer cells (the Warburg effect) and has
been a scientific challenge for decades. Here we show, using flux balance analysis combined with in vitro
measured enzyme specific activities, that fermentation is more catalytically efficient than respiration,
i.e. it produces more ATP per protein mass. And that the switch to fermentation at high growth rates
therefore is a consequence of a high ATP production rate, provided by a limited pool of enzymes. The
catalytic efficiency is also higher for cells grown on glucose compared to galactose and ethanol, which
may explain the observed differences in their growth rates. The enzyme F1F0-ATP synthase (Complex V)
was found to have flux control over respiration in the model, and since it is evolutionary conserved, we
expect the trade-off to occur in organisms from all kingdoms of life.


- KeyWords: Protein allocation, FBA, Crabtree Effect

**Utilisation:** maximising ATP, maximising growth, predictive simulation; **Model Source:** iFF708; **Taxonomy:** Saccharomyces Cerevisiae **Condition:** Chemostat 

- Reference: Submitted

- Pubmed ID: 26928598

- Last update: 2019-03-20



This repository is administered by @avlant, Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology


## Installation
The repo contains scripts that run under Matlab, no installation is required. Open files in matlab and press the green play button, more details [here](https://se.mathworks.com/help/matlab/matlab_prog/create-scripts.html).

### Required Software:
Code has been tested under windows 7 using matlab and RAVEN.

* *_PROGRAMMING LANGUAGE/Version_*  (e.g.):
  * You need a functional Matlab installation of **Matlab_R_2018_b**
  * (optional, required for some functionality) The [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox for MATLAB. 
  
### Run
First run runMeFirst.m to set up the path to the embeded RAVEN functions  
To reproduce the figures of the paper run corresponding files (e.g. figure1c_XYZ.m for figure 1C)



### Addaption of the software to your data
The code is not intended to be used as a tool, but it is written in modular form, so minimal effort should in principle be required to replace parameter values, experimental data files, and stoichiometric models with custom versions. Benevolent behavior of the program cannot be guaranteed.