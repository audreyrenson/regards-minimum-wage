sl_wrapper_aipw = function(formula, data, family=binomial(), library=q_library, ...) {
  sl_wrapper(formula    = formula,
             family = family,
             data       = select(data, t, female:sweight, ...),
             SL.library = library,  weights = sweight,
             cvControl  = list(V=4), parallel = FALSE)
}
predict_aipw = function(object, data, ...) {
  predict(object, newdata=data %>% intervene(),
          type='response')
}
