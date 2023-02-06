require(glmnet)
glmnet_f = function(formula, data, weights, subset, ...) {

  if(!missing(subset)) {
    subset = rlang::enexpr(subset)
    subset = eval(subset, envir = data)
    data = data[subset, ]
  }

  if(missing(weights)) {
    weights = rep(1, nrow(data))
  } else {
    weights = rlang::enexpr(weights)
    weights = eval(weights, envir = data)
  }
  keep = apply(model.frame(formula, data, na.action = na.pass), 1, function(x) all(!is.na(x)))
  weights = weights[keep]


  x = model.matrix(formula, data)
  y = model.frame(formula, data)[,1]
  mod = glmnet(x,y, weights=weights, ...)
  #mod = structure(mod, class = c('glmnet_f', class(mod)))
  mod$formula = formula
  mod$data = data
  mod
}

predict.glmnet_f = function(object, newdata, type='response', ...) {
  predict.glmnet(object, newx = model.matrix(object$formula, newdata), type=type, ...)
}
