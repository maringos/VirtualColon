# Virtual Colon
Spatiotemporal modelling of metabolic interactions in a computational colonic environment
------------------------------------------------------------------------------------------------
This repository contains the source and input files needed to run the Virtual Colon. 

- The Mono-colonalisation.R file simulates mono-colonisation experiments with context-specific colonic cell.
  
- The SIHUMIx.R file context-specific colonic cells with the SIHUMIx minimal model microbiome.

It requires the following software / R packages:
1. R
2. sybil
3. BacArena
4. cplexAPI

The gplkAPI would theoritically suffice instead of cplexAPI. However, it is not tested and may lead to significant increase of computational time.

## References and acknowlegments
This research was supported in part through high-performance computing resources available at the Kiel University Computing Centre.

1. Eugen Bauer, Johannes Zimmermann, et al..,  BacArena: Individual-based metabolic modeling of heterogeneous microbes in complex communities, doi:10.1371/journal.pcbi.1005544
2. R Core Team. R: A Language and Environment for Statistical Computing [Internet]. Vienna, Austria: R Foundation for Statistical Computing; Available from: https://www.R-project.org/
3.	Gelius-Dietrich G, Desouki AA, Fritzemeier CJ, Lercher MJ. sybil – Efficient constraint-based modelling in R. BMC Syst Biol. 2013 Dec;7(1):125.

## Citation
Georgios Marinos*, Johannes Zimmermann*, Jan Taubenheim, and Christoph Kaleta, **Virtual Colon: Spatiotemporal modelling of metabolic interactions in a computational colonic environment**, 2024, in preparation

*shared first authorship with an interchangeable order
