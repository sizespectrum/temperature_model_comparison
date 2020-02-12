#' tempFun is a function that takes temperature parameters (temperature, t_ref, t_d) and physiological parameters (w, Ea, c_a)
#' and returns a scalar based on the Padfield (2016) temperature performance equation. This will be upgraded in the future with
#' a deactivation portion with Anna Gardmark's equation once it's published
tempFun <- function(w, temperature, t_d = 25, t_ref, Ea, c_a) # default are 0 for now as deactivation is buggy
{
  k = 8.617332e-5 # Boltzmann constant 
  # equation
  # (w^(c_a*(temperature-t_ref)))  *exp((-Ea/8.617332e-5)*((1/temperature) - (1/t_ref)))
  
  # converting to Kelvin from Celcius
  temperature <- temperature + 273 
  t_ref <- t_ref + 273
  t_d <- t_d + 273
  
  temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature-(t_ref)))}) *exp((-Ea/k)*((1/temperature) - (1/(t_ref)))))
  return(temperatureScalar)
}

#' temperatureMizer is a function that uses the latest version of Mizer (02/2020) and changes 4 parameters that are temperature sensitive
#' (search volume / max intake / metabolism / background mortality) according to the scalars from tempFun.
#' At the moment this function only allow time independent changes (i.e temperature is constant during the simulation and the affected parameters
#' do not vary through time)
#' this function takes and return a Mizer params object
#' 
temperatureMizer <- function(params, t_ref = 10, temperature = t_ref, 
                             ea_int = 0,ca_int = 0, ea_met = 0, ca_met = 0, ea_mor = 0, ca_mor = 0)
{
  intakeScalar <- tempFun(w = params@w, t_ref = t_ref, temperature = temperature, Ea = ea_int, c_a = ca_int)
  metScalar <- tempFun(w = params@w, t_ref = t_ref, temperature = temperature, Ea = ea_met, c_a = ca_met)
  mortScalar <- tempFun(w = params@w, t_ref = t_ref, temperature = temperature, Ea = ea_mor, c_a = ca_mor)
  
  params@search_vol <- sweep(params@search_vol,2, intakeScalar, "*", check.margin = FALSE)
  params@intake_max <- sweep(params@intake_max,2, intakeScalar, "*", check.margin = FALSE)  
  params@metab <- sweep(params@metab,2, metScalar, "*", check.margin = FALSE)
  params@mu_b <- sweep(params@mu_b,2, mortScalar, "*", check.margin = FALSE)
  return(params)
}