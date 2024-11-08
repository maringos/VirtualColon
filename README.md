# Virtual Colon
Spatiotemporal modelling of metabolic interactions in a computational colonic environment
------------------------------------------------------------------------------------------------
This repository contains the source and input files needed to run the Virtual Colon. 

- The Mono-colonalisation.R file simulates mono-colonisation experiments with context-specific colonic cells and bacterial ones.
  
- The SIHUMIx.R file simulates context-specific colonic cells with the SIHUMIx minimal model microbiome.

It requires the following software / R packages:
1. R
2. sybil
3. BacArena
4. cplexAPI

The gplkAPI would theoritically suffice instead of cplexAPI. However, it is not tested and may lead to significant increase of computational time.

# Example of usage
You can try our code with the example datasets. To this end, you have to choose the required models from the file `simulations_sequence.RDS` by selecting the respective row.

For instance, please use the command `Rscript SIHUMIx.R $nodes` in the command line, where $nodes is the number of the selected row in `simulations_sequence.RDS`.

## References and acknowlegments
This research was supported in part through high-performance computing resources available at the Kiel University Computing Centre. See our manuscript for a full list of other work cited.

1. Eugen Bauer, Johannes Zimmermann, et al.,  BacArena: Individual-based metabolic modeling of heterogeneous microbes in complex communities, doi:10.1371/journal.pcbi.1005544
2. R Core Team. R: A Language and Environment for Statistical Computing [Internet]. Vienna, Austria: R Foundation for Statistical Computing; Available from: https://www.R-project.org/
3.	Gelius-Dietrich G, Desouki AA, et al., sybil – Efficient constraint-based modelling in R. BMC Syst Biol. 2013 Dec;7(1):125.
4.	Becker N., Kunath J., et al., Human intestinal microbiota: Characterization of a simplified and stable gnotobiotic rat model. Gut Microbes 2, 25–33 (2011).
5.	Geva-Zatorsky N., Sefik E., et al.,Mining the Human Gut Microbiota for Immunomodulatory Organism,  Cell. 2017 Feb 23;168(5):928-943.e11


## Citation
Georgios Marinos*, Johannes Zimmermann*, Jan Taubenheim, and Christoph Kaleta, **Virtual Colon: Spatiotemporal modelling of metabolic interactions in a computational colonic environment**, 2024, www.doi.org/10.1101/2024.06.11.598488 

*shared first authorship with an interchangeable order

GNU General Public License version 3.0 (GPLv3) is applied to all copyrightable parts of this software. 

Contact details: https://www.iem.uni-kiel.de/de/medizinische-systembiologie/medizinische-systembiologie
