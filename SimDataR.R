##### Applied Spatial Statistics --- generate postive and negative spatial autocorrelation simulation data
### Erjia Ge,   Winter, 2019 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## 1. Generate the observations(y) ##
# xy-axis:
u <- seq(0, 5, by=0.25)
v <- seq(0, 5, by=0.25)
xy_pts <- cbind(rep(u, each=21), rep(v, 21))

# 1) positive spatiotemporal autocorrelation signals + white noise (mean = 0, deviation = 0.5)
Sig.Fn_Po <- function(k) {
  # Positive signal function at time k:
  sig.fn_po <- numeric(0)
  for (i in u) {
    for (j in v) {
      s <- (1/(1+0.5*k^2))*(exp(((i-1)^2+(j-1)^2-2.5*(1+k^2))/(-(1+0.75*k^2)))
                            +exp(((i-1)^2+(j-4)^2-2.5*(1+k^2))/(-(1+0.75*k^2)))
                            +exp(((i-4)^2+(j-1)^2-2.5*(1+k^2))/(-(1+0.75*k^2)))
                            +exp(((i-4)^2+(j-4)^2-2.5*(1+k^2))/(-(1+0.75*k^2))))
      sig.fn_po <- c(sig.fn_po, s)
    }
  }
  po <- sig.fn_po 
  
  # white noise
  noise <- rnorm(length(po), 0, 0.2^2)
  po <- po + noise
  
  m <- cbind(xy_pts, po)
  return(m)
}

# The observations generated below are in the format of (x,y,obs).
po0 <- cbind(Sig.Fn_Po(0),0)
po1 <- cbind(Sig.Fn_Po(1),1)
po2 <- cbind(Sig.Fn_Po(2),2)
po3 <- cbind(Sig.Fn_Po(3),3)
po <- rbind(po0,po1,po2,po3)
colnames(po) <- c('u', 'v', 'signoi', 't')

# Visualization:
po0.df <- matrix(po0[,3], nrow=21, ncol=21)
po1.df <- matrix(po1[,3], nrow=21, ncol=21)
po2.df <- matrix(po2[,3], nrow=21, ncol=21)
po3.df <- matrix(po3[,3], nrow=21, ncol=21)

library(plotly)

po0.plot <- add_surface(plot_ly(x=u, y=v, z=~po0.df))
po1.plot <- add_surface(plot_ly(x=u, y=v, z=~po1.df))
po2.plot <- add_surface(plot_ly(x=u, y=v, z=~po2.df))
po3.plot <- add_surface(plot_ly(x=u, y=v, z=~po3.df))

# 2) negative spatiotemporal autocorrelation:
Sig.Fn_Ne <- function(k) {
  
  # Negative signal function at time k:
  
  # notice** the 1-10 inter-value matrix is only true when the count of column is odd.
  sig.fn_ne <- matrix(rep(c(1,10), len=441), ncol = 21) 
  
  for (i in u) {
    for (j in v) {
      if ((abs(i-2) <= 0.5+0.25*k && abs(j-3.25) <= 0.5+0.25*k) ||
          (abs(i-3) <= 0.5+0.25*k && abs(j-1.75) <= 0.5+0.25*k)) {}
        else {
          sig.fn_ne[which(i==u), which(j==v)] <- 5.5
        }
    }
  }
  
  # noise 
  noise <- matrix(rnorm(length(sig.fn_ne), 0, 0.2^2), ncol = 21)
  
  sig.fn_ne <- as.numeric(sig.fn_ne+noise)
  
  m <- cbind(xy_pts, sig.fn_ne)
  return(m)
}

ne0 <- cbind(Sig.Fn_Ne(0),0)
ne1 <- cbind(Sig.Fn_Ne(1),1)
ne2 <- cbind(Sig.Fn_Ne(2),2)
ne3 <- cbind(Sig.Fn_Ne(3),3)
ne <- rbind(ne0,ne1,ne2,ne3)
colnames(ne) <- c('u', 'v', 'signoi', 't')


# Visualization:
ne0.df <- matrix(ne0[,3], nrow=21, ncol=21)
ne1.df <- matrix(ne1[,3], nrow=21, ncol=21)
ne2.df <- matrix(ne2[,3], nrow=21, ncol=21)
ne3.df <- matrix(ne3[,3], nrow=21, ncol=21)

ne0.plot <- add_surface(plot_ly(x=u, y=v, z=~ne0.df))
ne1.plot <- add_surface(plot_ly(x=u, y=v, z=~ne1.df))
ne2.plot <- add_surface(plot_ly(x=u, y=v, z=~ne2.df))
ne3.plot <- add_surface(plot_ly(x=u, y=v, z=~ne3.df))

save(ne, po, file = "Signals.Rda")

