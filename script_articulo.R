# Title: Script of GAMLSS models

# Authors:
# - Yeison Yovany Ocampo Naranjo
# - Laura Nataly Montoya Colorado
# - Freddy Hernández  Barajas

# Date: 05/16/2018


# Libraries neccesariess --------------------------------------------------

library(gamlss)
library(moments)

# Lectura de datos --------------------------------------------------------

tiempo_de_exposicion<- c(50,70,50,70,16.36,50,50,50,83.64,70,50,30,30,70,30,50,50,30)
densidad_de_corriente<- c(41.67,53.03,41.67,30.30,41.67,41.67,60.78,41.67,41.67,30.30,22.55,53.03,53.03,53.03,30.30,41.67,41.67,30.30)
temperatura<- c(50,45,50,55,50,50,50,41.6,50,45,50,45,55,55,45,50,58.41,55)
espesor_promedio<- c(8.12,21.31,8.12,16.73,5.58,10.52,12.36,14.20,22.41,13.23,10.91,11.62,12.80,17.70,8.44,11.49,11.63,7.67)

# variables names (easy to use)

y<-espesor_promedio
x1<-tiempo_de_exposicion
x2<-densidad_de_corriente
x3<-temperatura

# Dataframe with the dataset

nombres<- c("Espesor","Tiempo","Dens","Temp")
datos<- data.frame(espesor_promedio,
                   tiempo_de_exposicion,
                   densidad_de_corriente,
                   temperatura)

names(datos)<- nombres



# Moments -----------------------------------------------------------------

skewness(y)
kurtosis(y)

# The first fitted model by Ruiz et al.(2003) ------------------------------------

modelo.0 <- lm(y~x1+I(x1^2)+x2+x3)
coef(modelo.0)
summary(modelo.0)

plot(modelo.0)
shapiro.test(residuals(modelo.0))
# GAMLSS section (selecting distribution) ---------------------------------

library(gamlss)
m1<-fitDist(y,type = "realplus")
m1$fit

# formula to the full model (possible mistake)

form<- as.formula("~x1*x2+x1*x3+x2*x3+I(x1^2)*I(x2^2)+I(x1^2)*I(x3^2)+I(x2^2)*I(x3^2)")
form

form2 <- as.formula("~x1*x2+x1*x3+x2*x3+I(x1^2)+I(x2^2)+I(x3^2)")
form2
# Fitted models -----------------------------------------------------------

## Only intercept models

### Inverse Gaussian model
modelog1<- gamlss(y~1,family =IG, trace=F,n.cyc = 50)
summary(modelog1)

### Log-normal model
modelog2<- gamlss(y~1,family =LOGNO, trace=F,n.cyc = 50)
summary(modelog2)

### Inverse Gamma model
modelog3<- gamlss(y~1,family =IGAMMA, trace=F,n.cyc = 50)
summary(modelog3)

### Gamma model
modelog4<- gamlss(y~1,family =GA, trace=F,n.cyc = 50)
summary(modelog4)


# Finals model using GAIC.ALL() Forward step ------------------------------

### Inverse Gaussian model
modelog1.final<-stepGAICAll.A(modelog1,
                              scope=list(lower=~1,
                                         upper=form,
                                         k=log(length(y))),trace=F)

### Log-normal model
modelog2.final<-stepGAICAll.A(modelog2,
                              scope=list(lower=~1,
                                         upper=form,
                                         k=log(length(y))),trace=F)

### Inverse Gamma model
modelog3.final<-stepGAICAll.A(modelog3,
                              scope=list(lower=~1,
                                         upper=form,
                                         k=log(length(y))), trace=F)

### Gamma model
modelog4.final<-stepGAICAll.A(modelog4,
                              scope=list(lower=~1,
                                         upper=form,
                                         k=log(length(y))),trace=F)

# Generalized Akaike Information Criterio ---------------------------------

AIC(modelo.0,k = 2)
GAIC(modelog1.final, k = 2) # Inverse Gaussiana
GAIC(modelog2.final, k = 2) # Log-normal
GAIC(modelog3.final, k = 2) # Inverse gamma
GAIC(modelog4.final, k = 2) # Gamma


# Correlation between real y and fitted y
esp.hat0 <- fitted(modelo.0)
cor(y,esp.hat0)
esp.hat1 <- fitted(modelog1.final, what='mu') 
cor(y,esp.hat1)
esp.hat2<-exp(fitted(modelog2.final,what="mu"))*((exp(fitted(modelog2.final,what="sigma")))^1/2)
cor(y,esp.hat2)
esp.hat4<-fitted(modelog4.final,what="mu")
cor(y,esp.hat4)
esp.hat3<-fitted(modelog3.final,what = "mu")
cor(y,esp.hat3)

# Validation models (worm plot)

par(mfrow=c(2,2))
wp(modelog3.final, main="Gráfico de gusano usando la distribución Gamma")
title("Gamma")
wp(modelog2.final, main="Gráfico de gusano usando la distribución Lognormal")
title("Log-Normal")
wp(modelog1.final, main="Gráfico de gusano usando la distribución ")
title("Inversa Gaussiana")
wp(modelog4.final, main="Gráfico de gusano usando la distribución ")
title("Inversa Gamma")

par(mfrow=c(1,1))


# Final model (using the Gamma distribution) -------------------------------------------------------------

modelo.final<-gamlss(y~x1+x2+I(x1^2)+I(x2^2),
                     sigma.formula = ~x1+I(x1^2),
                     family = GA, trace=F,n.cyc = 50)
summary(modelo.final)


# Surface response to mean of model ---------------------------------------

modcom<-gamlss(y~x1+x2+I(x1^2)+I(x2^2),
               sigma.formula = ~x1+I(x1^2),
               family = GA,
               trace=F,n.cyc = 50)

# grid to values of x1,x2
X1  <- seq(from=16.36, to=83.64, length.out=15)
X2 <- seq(from=5.58, to=22.41, length.out=15)

f.espesor <- function(x1, x2) {
  x <- c(1,x1,x2, x1^2, x2^2)
  exp(sum(coef(modcom) * x))
}

f.espesor <- Vectorize(f.espesor)

Z <- outer(X1, X2, f.espesor)

pdf(file = "superficie1.pdf",width = 7,height = 7)
persp(x=X1, y=X2, z=Z,
      theta=310, phi=0, ticktype = "detailed",
      col='dodgerblue',
      main = expression(widehat(E)(Y)),
      xlab = "Tiempo",ylab = "Densidad",
      zlab = "Espesor promedio estimado \n ")
dev.off()


# Surface response to the variance of the model ---------------------------

# grid to x1,x2
X1  <- seq(from=16.36, to=83.64, length.out=15)
X2 <- seq(from=5.58, to=22.41, length.out=15)

f_var.espesor <- function(x1, x2) {
  x <- c(1,x1,x2, x1^2, x2^2)
  beta <- c(-2.0092,0.01782,0.00202,-0.0512,0.0008)
  (sum(beta * x))
}

f_var.espesor <- Vectorize(f_var.espesor)

Z1 <- outer(X1, X2, f_var.espesor)


# initial angle (310)
pdf(file = "superficie_varianza_sin_exp.pdf",
    width = 7,height = 7)
persp(x=X1, y=X2, z=Z1,
      theta=310, phi=-1, ticktype = "detailed",
      col='dodgerblue',
      main = expression(log(widehat(Var)(Y))),
      xlab = "Tiempo",ylab = "Densidad",
      zlab = "Espesor promedio estimado\n",
      cex.axis=0.8)

dev.off()

