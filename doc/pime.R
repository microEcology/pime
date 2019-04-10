## ------------------------------------------------------------------------
library(pime)
data("restroom")
pime.oob.error(restroom, "Environment")

## ------------------------------------------------------------------------
per_variable_obj= pime.split.by.variable(restroom, "Environment")
per_variable_obj

## ------------------------------------------------------------------------
prevalences=pime.prevalence(per_variable_obj)
prevalences

## ------------------------------------------------------------------------
pime.best.prevalence(prevalences, "Environment")

## ------------------------------------------------------------------------
prevalence.60 = prevalences$`60`
prevalence.60

## ------------------------------------------------------------------------
#randomized=pime.error.prediction(restroom, "Environment", bootstrap = 10, parallel = TRUE, max.prev = 95)
#randomized$Plot
#randomized$`Table results'

## ------------------------------------------------------------------------
replicated.oob.error= pime.oob.replicate(prevalences, "Environment", bootstrap = 10, parallel = TRUE)

