## ------------------------------------------------------------------------
library(pime)
data("restroom")
pime.oob.error(restroom, "Environment")

## ------------------------------------------------------------------------
per_variable_obj <- pime.split.by.variable(restroom, "Environment")
per_variable_obj

## ------------------------------------------------------------------------
prevalences <- pime.prevalence(per_variable_obj)
prevalences

## ------------------------------------------------------------------------
set.seed(42)
best.prev=pime.best.prevalence(prevalences, "Environment")

## ------------------------------------------------------------------------
imp50=best.prev$`Importance`$`Prevalence 50`
knitr::kable(imp50) %>% kableExtra::kable_styling(full_width = F)
#To get the table with OOB error results.
#best.prev$`OOB error`

## ------------------------------------------------------------------------
#randomized=pime.error.prediction(restroom, "Environment", bootstrap = 10, parallel = TRUE, max.prev = 95)
#randomized$Plot
#randomized$`Table results'

## ------------------------------------------------------------------------
#replicated.oob.error= pime.oob.replicate(prevalences, "Environment", bootstrap = 10, parallel = TRUE)

