# R implementation for Johnson & Berry photosynthesis model

The script [johnson_model.R](johnson_model.R) contains code for the model function. The function will calculate gross rate of CO2 assimilation, potential rate of net CO2 assimilation under Rubisco limitation,  potential rate of net CO2 assimilation under Cyt b6f limitation, and fluorescence parameters.

The script [test_johnson_model.R](test_johnson_model.R) will test the function. It is suggested that users run this script first to ensure the function will run properly.

The folder [example_figures](example_figures) contains figures showing the PAR response for both the static and dynamic models that are created in the script [test_johnson_model.R](test_johnson_model.R).

All function descriptions, including parameter descriptions, can be found in the script files.

## Inputs
- PAR = PAR, umol PPFD m-2 s-1
- Temp = Leaf temperature, C
- CO2 = Mesophyll CO2, ubar
- O2 = Atmospheric O2, mbar
- Abs = Total leaf absorptance to PAR, mol PPFD absorbed mol-1 PPFD incident
- beta = PSII fraction of total leaf absorptance, mol PPFD absorbed by PSII mol-1 PPFD absorbed
- CB6F = Cyt b6f density, mol sites m-2
- RUB = Rubisco density, mol sites m-2
- Rds = Scalar for dark respiration, dimensionless
- Ku2 = Rate constant for exciton sharing at PSII, s-1
- theta1 = Curvature parameter for Aj/Ac transition, dimensionless
- eps1 = PSI transfer function, mol PSI F to detector mol-1 PSI F emitted
- eps2 = PSII transfer function, mol PSII F to detector mol-1 PSII F emitted
- alpha_opt: option for static or dymanic absorption cross-sections of PSI and PSII

## Citation

Johnson, J. E. and J. A. Berry. 2021. The role of Cytochrome b<sub>6</sub>f in the 
control of steady-state photosynthesis: a conceptual and quantitative model.
 *Photosynthesis Research*, DOI: 10.1007/s11120-021-00840-4

## Contact

Any questions or issues can be submitted via GitHub or directed to Risa McNellis (risa.mcnellis@ttu.edu).