Model description:
Stan code that is based on "Comparison of chromatographic stationary phases using a Bayesian-based multilevel model". 
Simplified version of stan code that uses pKa as a predictor
one sigma
pKas and alphas without etas
between analyte variabilty for dlogkw and dS1
dlogkw ans dS1 the same for acids and bases

noncentered parametrication - seems ok
gaussain process matrix normal distribution
The max correlation set up as 0.8. The one calculated based on similarity is 1 for large group of compounds.
