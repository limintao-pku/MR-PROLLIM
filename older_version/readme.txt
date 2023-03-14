MR-PROLLIM (version 1.0.1) was used to derive the results in our paper.

The 1.0.2 version has the following main changes:
1. Version 1.0.2 has changed the behavior of the random number generator (RNG). The previous version produces
    reproducible results only with the "L'Ecuyer-CMRG" kind, while this version is compatible with other RNG kinds.
2. Version 1.0.2, by default, uses both the nlminb and nlm optimizers for log-linear regressions, while the previous
    version mainly uses the nlminb optimizer. According to our simulations, nlminb is less likely to get stuck in regions
    that produce > 1 probabilities than nlm. But it does not have a gradient control as good as nlm, thus leading to
    less precise estimates. MR-PROLLIM (version 1.0.2) first uses nlminb to derive initial estimates and then refine them
    with nlm (the initial results are put into nlm as the starting point).
3. Version 1.0.2 runs the genoud optimizer in the global environment (the user's workspace) to speed up parallel
    computations, while the previous version does this in the local environment of MR-PROLLIM. The default pop.size
    for genoud is also increased to a > 1000 value if > 3 logical cores are available.

We expect that these changes will not have substantial impacts on the results in our paper.    

R scripts used to derive the results in our paper are also provided.