# Create a random genetic dataset of size nsamples x nvars
random.genetic.dataset <- function(nsamples,nvars,characterSize){
  # loop through, adding a random vector as a column each time
  D <- list()
  for (i in 1:(nvars-1)){
    D <- c(D,list(sample(0:(characterSize-1), nsamples, replace = TRUE)))
  }
  D <- as.data.frame(D)
  colnames(D) <- 1:(nvars-1)
  # add class
  D$Class <- sample(0:1,nsamples,replace=TRUE)
  return(D)
}
