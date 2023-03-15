##BF function based on Dienes & McLatchie, 2018
##modified so that H1 is represented by normal distribution (rather than t), hence there is no 'dftheory' argument

Bf<-function(sd, obtained, dfdata, meanoftheory, sdtheory, tail = 2)
{
  area <- 0
  normarea <- 0
  theta <- meanoftheory - 5 * sdtheory
  incr <- sdtheory / 200
  for (A in -1000:1000){
    theta <- theta + incr
    dist_theta <- dnorm((theta-meanoftheory)/sdtheory)
    if(identical(tail, 1)){
      if (theta <= 0){
        dist_theta <- 0
      } else {
        dist_theta <- dist_theta * 2
      }
    }
    height <- dist_theta * dt((obtained-theta)/sd, df = dfdata)
    area <- area + height * incr
    normarea <- normarea + dist_theta*incr
  }
  LikelihoodTheory <- area/normarea
  Likelihoodnull <- dt(obtained/sd, df = dfdata)
  BayesFactor <- LikelihoodTheory / Likelihoodnull
  BayesFactor
}
