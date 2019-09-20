# nias-et-al-2019-grl

Code used to create the ensemble of BISICLES simulations.

*gen-ase-axby-ensembleB.R* creates input data files for the perturbed-parameter ensemble using Latin hypercube sampling. It also produces BISICLES configuration files based on *inputs.template*.
The resulting files were used to run the ensemble presented in Nias *et al*. (2016), which informed the prior sea level distribution in this study.

For the 200-year extended simulations, the melt rate forcing is applied using *forcing.py* (the relevant lines for the basal flux are commented out in *inputs.template*).

#### References

Nias, IJ, Cornford, SL, & Payne, AJ (2016) Contrasting the modelled sensitivity of the Amundsen Sea Embayment ice streams. *Journal of Glaciology*, 62(233), 552–562. doi: 10.1017/jog.2016.40
