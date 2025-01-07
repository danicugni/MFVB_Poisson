r.delta.invga <- function(media, varianza) {
  r <- media^2/varianza + 2
  delta <- media*(media^2/varianza + 1)
  return(list(r = r, delta = delta))
}