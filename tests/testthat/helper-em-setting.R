## Most tests check properties that hold at any EM tolerance: equivalence
## with the ported engine, invariances, the E-step and the samplers. They run
## at the pre-0.2.0 setting to keep the suite fast (about 7 minutes; about 75
## at the converged default). test-defaults.R clears these options and tests
## the current defaults (tol = 1e-8, max_iter = 1000).
options(gemcox.tol = 1e-4, gemcox.max_iter = 100L)
