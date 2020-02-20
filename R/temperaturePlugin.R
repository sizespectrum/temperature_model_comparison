#' tempFun is a function that takes temperature parameters (temperature, t_ref, t_d) and physiological parameters (w, Ea, c_a)
#' and returns a scalar based on the Padfield (2016) temperature performance equation. This will be upgraded in the future with
#' a deactivation portion with Anna Gardmark's equation once it's published
tempFun <- function(w, temperature, t_ref, Ea, c_a)
{
  k = 8.617332e-5 # Boltzmann constant 
  # equation
  # (w^(c_a*(temperature-t_ref)))  *exp((-Ea/8.617332e-5)*((1/temperature) - (1/t_ref)))
  
  # converting to Kelvin from Celcius
  temperature <- temperature + 273 
  t_ref <- t_ref + 273
  
  temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature-(t_ref)))}) *exp((-Ea/k)*((1/temperature) - (1/(t_ref)))))
  return(temperatureScalar)
}


#' Temperature Mizer plugin
#' 
#' uses the latest version of Mizer (02/2020) and changes 4 parameters that are temperature sensitive
#' (search volume / max intake / metabolism / background mortality) according to the scalars from tempFun.
#' At the moment this function only allow time independent changes (i.e temperature is constant during the simulation and the affected parameters
#' do not vary through time)
#' this function takes and return a Mizer params object
#' 
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param scalar Default is set to false. If true, will return the value of the different scalars in a list instead of the modified param object
#' @param t_ref Temperature where scalar equals 1. Default is 10 degrees C.
#' @param temperature Ambient temperature influencing the scalar. Only an integer at the moment meaning that temperature do not vary
#' during projections. Default is set to t_ref (i.e. temperature has no effect)
#' @param ea_int Activation energy relative to food intake. This parameter affects the search volume and intake max
#' Default is set to 0 (i.e. no effect)
#' @param ca_int Size dependence activation relative to food intake. Negative values flatten the performance curve slope with increasing size.
#' Default is set to 0 (i.e. no size dependence)
#' @param ea_met Activation energy relative to metabolism
#' @param ca_met Size dependence activation relative to metabolsim
#' @param ea_mor Activation energy relative to background mortality
#' @param ca_mor Size dependence activation relative to mortality
#' 
#' @return A \linkS4class{MizerParams} object
#' @export

temperatureMizer <- function(params, scalar = F, t_ref = 10, temperature = t_ref, 
                             ea_int = 0,ca_int = 0, ea_met = 0, ca_met = 0, ea_mor = 0, ca_mor = 0)
{
  intakeScalar <- tempFun(w = params@w, t_ref = t_ref, temperature = temperature, Ea = ea_int, c_a = ca_int)
  metScalar <- tempFun(w = params@w, t_ref = t_ref, temperature = temperature, Ea = ea_met, c_a = ca_met)
  mortScalar <- tempFun(w = params@w, t_ref = t_ref, temperature = temperature, Ea = ea_mor, c_a = ca_mor)
  
  params@search_vol <- sweep(params@search_vol,2, intakeScalar, "*", check.margin = FALSE)
  params@intake_max <- sweep(params@intake_max,2, intakeScalar, "*", check.margin = FALSE)  
  params@metab <- sweep(params@metab,2, metScalar, "*", check.margin = FALSE)
  params@mu_b <- sweep(params@mu_b,2, mortScalar, "*", check.margin = FALSE)
  
  if(scalar) return(list("intake" = intakeScalar, "metabolsim" = metScalar, "mortality" = mortScalar)) else return(params)
}