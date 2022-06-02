# Eliminación de todos los objetos anteriores
rm(list=ls(all=TRUE))

### ----------------------------------------------------------------------
### Paquetes y funciones necesarias
### ----------------------------------------------------------------------

# Paquetes

library(forecast)
library(TSA)
library(car)
library(lmtest)
library(FitAR)

# Funciones de usuario

source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loess.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexpo.ErrorARMA.R")
### ----------------------------------------------------------------------
### 1. Lectura de datos
### ----------------------------------------------------------------------

# Lectura de la base de datos
enlace <- "C:/Users/Sofia/Documents/GitHub/E3_T2/Sofía/anexos-emmet-noviembre-2021-1-total industria-modif.csv"
datos <- read.table(file = enlace,
                    header=T,
                    sep=";",
                    skip=15,
                    dec=",")
# Selección de la columna necesaria
datos <- as.data.frame(datos[, 7])
colnames(datos) <- c("ventas.nominal")

# Creación del objeto ts
datos <- ts(datos,
            freq=12,
            start=c(2001,1))

##
##
##

### ----------------------------------------------------------------------
### 2. Análisis descriptivo
### ----------------------------------------------------------------------

# Gráfica de la serie de tiempo lineal

par(adj = 0.5, col = 'black')
plot(datos,
     xlab = "Año\n",
     ylab = "Índice de ventas nominal",
     cex.main = 1,
     lwd = 1)
grid(col = 'gray', lwd = 1)
par(adj = 1,
    col = 'black')

# Gráfica de la serie de tiempo en escala logarítmica

par(adj = 0.5, col = 'black')
plot(log(datos),
     xlab = "Año\n",
     ylab = "Logaritmo del índice de ventas nominal",
     cex.main = 1,
     lwd = 1)
grid(col = 'gray', lwd = 1)
par(adj = 1,
    col = 'black')

# Barplot para la estacionalidad

par(adj = 0.5, col = 'black')
boxplot(log(datos)~cycle(log(datos)),
        names = month.abb,
        cex.axis = 0.8,
        main = "Valor del índice mensual de ventas nominal
        en Colombia por mes",
        xlab = "Mes",
        ylab = "Índice mensual de ventas nominal")
par(adj = 1)
title(sub = "Fuente. DANE")

# Periodograma
par(adj = 0.5)
log.desest <- diff(log(datos))
periodogram(log.desest,
            lwd = 2,
            main = "Periodograma para el índice de ventas nominal
            en Colombia por mes",
            xlab = "Frecuencias",
            ylab = "Asociación")
abline(v = c(1:6)/12,
       col=2, lty=2)

# Gráfico de la componente de tendencia de la gráfica

par(adj = 0.5)
plot(decompose(log(datos))$trend,
     xlab = "Año",
     ylab = "Logaritmo del índice mensual de ventas nominales",
     main = "Tendencia del logaritmo índice mensual de ventas
     nominales en Colombia entre 2001 y 2021",
     lwd = 2)
par(adj = 1,
    col = 'black')
title(sub = "Fuente. DANE. Estimación de la tendencia hecha con R.")
grid(col = 'gray', lwd = 1)

# Gráfico de la función de autocorrelación muestral (ACF) del logaritmo

par(adj = 0.5)
acf(as.numeric(log(datos)),
    lag.max = round((length(datos)-12)/4),
    ci.type = "ma",
    col = 1,
    ci.col = 2,
    main = '')

# Gráfico de la componente estacional de la serie

descomposicion <- decompose(datos, type = "multiplicative")

par(adj = 0.5)
plot(descomposicion$seasonal,
     main = "Componente estacional del índice mensual de
     ventas nominales en Colombia entre 2001 y 2021",
     xlab = "Año",
     ylab = "Índice mensual de ventas nominales")
par(adj = 1,
    col = 'black')
title(sub = "Fuente. DANE. Estimación de la estacionalidad hecha con R.")

par(adj = 1,
    col = 'black')

##
##
##

### ----------------------------------------------------------------------
### 3. Preparación para el modelo exponencial polinomial estacional
### ----------------------------------------------------------------------

# Definición de variables necesarias

m <- 12                # Número de periodos a pronosticar dentro de la muestra
n <- length(datos) - m # Tamaño de la muestra para el ajuste
t <- 1:n               # Índice de tiempo en los periodos de ajuste
t2 <- t^2
t3 <- t^3
t4 <- t^4
t5 <- t^5
t6 <- t^6

# Funciones trigonométricas para la componente estacional

sen1 <- sin(pi*t/6)
cos1 <- cos(pi*t/6)
sen2 <- sin(pi*t/3)
cos2 <- cos(pi*t/3)
sen3 <- sin(pi*t/2)
cos3 <- cos(pi*t/2)
sen4 <- sin(2*pi*t/3)
cos4 <- cos(2*pi*t/3)
sen5 <- sin(5*pi*t/6)
cos5 <- cos(5*pi*t/6)

# Serie de tiempo con los n datos para construir el modelo

yt <- ts(datos[t],
         frequency = 12,
         start=c(2001,1))

# Matriz de diseño

X1 <- data.frame(t, t2, t3, t4, t5, t6,
                 sen1, cos1, sen2, cos2,
                 sen3, cos3, sen4, cos4,
                 sen5, cos5)

# Valores de las variables en la validación cruzada

tnuevo <- (n+1):length(datos) # Índice de tiempo en los pronósticos
t2nuevo <- tnuevo^2
t3nuevo <- tnuevo^3
t4nuevo <-  tnuevo^4
t5nuevo <-  tnuevo^5
t6nuevo <-  tnuevo^6

#Funciones trigonométricas en los periodos de pronóstico

sen1n <- sin(pi*tnuevo/6)
cos1n <- cos(pi*tnuevo/6)
sen2n <- sin(pi*tnuevo/3)
cos2n <- cos(pi*tnuevo/3)
sen3n <- sin(pi*tnuevo/2)
cos3n <- cos(pi*tnuevo/2)
sen4n <- sin(2*pi*tnuevo/3)
cos4n <- cos(2*pi*tnuevo/3)
sen5n <- sin(5*pi*tnuevo/6)
cos5n <- cos(5*pi*tnuevo/6)

# Serie de tiempo para los valores que van a ser usados para el ajuste

ytf <- ts(datos[tnuevo],
          freq = 12,
          start = c(2020, 12))

# Matriz de diseño con las potencias de t en polinomio de grado seis
# y trigonométricas para la estacionalidad

X1nuevo <- data.frame(t = tnuevo, t2 = t2nuevo, t3 = t3nuevo,
                      t4 = t4nuevo, t5 = t5nuevo, t6 = t6nuevo,
                      sen1 = sen1n, cos1 = cos1n, sen2 = sen2n,
                      cos2 = cos2n, sen3 = sen3n, cos3 = cos3n,
                      sen4 = sen4n, cos4 = cos4n, sen5 = sen5n,
                      cos5=cos5n)

# ------------------------------------------
# 4. Modelo global
#    Exponencial polinomial de grado seis con trigonométricas
# ------------------------------------------

# Crear vector con nombres de los parámetros a usar en la fórmula R 
# del modelo exponencial

param2 <- c(paste0("beta",0:6),
            "alfa1", "gamma1", "alfa2",
            "gamma2", "alfa3", "gamma3",
            "alfa4", "gamma4", "alfa5",
            "gamma5")
param2

modelo_global <- regexponencial(respuesta = yt,
                                data = X1,
                                names.param = param2)
summary(modelo_global)

# Cálculo valores ajustados del modelo 2

ythat2 <- ts(fitted(modelo_global),
             frequency = 12,
             start = start(yt))

# Pronósticos. Solo son posibles los pronósticos puntuales

predicciones2 <- predict(modelo_global,
                         newdata = X1nuevo,
                         interval = "prediction")
predicciones2

# Convirtiendo en serie de tiempo las predicciones
ytpron2 <- ts(predicciones2,
              frequency = 12,
              start = start(ytf))
ytpron2

# Exactitud de los pronósticos puntuales
accuracy(ytpron2, ytf)


# Cálculo AIC y BIC
npar2 <- length(coef(modelo_global)[coef(modelo_global)!=0])
npar2                                 # Número par?metros modelo 2b
Criterios2 <- exp.crit.inf.resid(residuales = residuals(modelo_global),
                                 n.par = npar2)
Criterios2

#Gráficos de residuos

par(adj = 0.5)
plot.ts(residuals(modelo_global),
        ylim = c(min(residuals(modelo_global),
                     -2*summary(modelo_global)$sigma,
                     2*summary(modelo_global)$sigma),
                 max(residuals(modelo_global),
                     -2*summary(modelo_global)$sigma,
                     2*summary(modelo_global)$sigma)),
        lwd = 1,
        xlab = "Periodo",
        ylab = "Residuales")
abline(h = c(-2*summary(modelo_global)$sigma,
             0,
             2*summary(modelo_global)$sigma),
       col = 2)
legend("topleft", legend = c("Modelo global"),
       lty = 1,
       col = 1,
       lwd = 2)


plot(fitted(modelo_global),
     residuals(modelo_global),
     ylim = c(min(residuals(modelo_global),
                  -2*summary(modelo_global)$sigma,
                  2*summary(modelo_global)$sigma),
              max(residuals(modelo_global),
                  -2*summary(modelo_global)$sigma,
                  2*summary(modelo_global)$sigma)),
     xlab = 'Índice de ventas nominales ajustado',
     ylab = 'Residuales')
abline(h = c(-2*summary(modelo_global)$sigma,
             0,
             2*summary(modelo_global)$sigma),
       col = 2)
legend("topleft",
       legend = c("Modelo global"),
       lty = 1,
       col = 1,
       lwd = 2)

# Gráfico del ajuste modelo

plot(datos,
     xlab = "Año",
     ylab = "Índice de ventas nominales")
lines(ythat2, col=2)
grid(col = 'gray', lwd = 1)
legend("topleft",
       legend = c("Original",
                  "Ajuste del modelo exponencial\npolinomial de grado seis estacional"),
       col = c(1, 2),
       lty = 1)



### ----------------------------------------------------------------------
### 5. Pronósticos del modelo global
### ----------------------------------------------------------------------

# Comparación gráfica de los pronósticos

vector.auxiliar <- c(ytf, ytpron2)
par(adj = 0.5)
plot(datos, lwd = 2,
     xlab = "Año",
     ylab = "Índice de ventas nominales")
lines(ythat2, col = 2, lwd = 2)
lines(ytpron2, col = 4, lwd = 2)
grid(col = 'gray', lwd = 1)
legend('topleft',
       legend = c('Original', 'Ajustada', 'Pronósticos'),
       lty = 1,
       lwd = c(1, 2, 2),
       col = c(1, 2, 4))



plot(ytf, type = "b", pch = 19, lty = 1, col = 1, lwd = 2,
     ylab = "Índice de ventas nominales",
     xlab = "Periodo [mmm - yy]",
     ylim = c(min(vector.auxiliar), max(vector.auxiliar)),
     xaxt = "n")
lines(ytpron2, col = 4, pch = 3, lty = 3, type = "b", lwd = 2)
grid(col = 'gray', lwd = 1)
legend("bottomright",
       legend = c("Real", "Modelo global exponencial\ncúbico con trigonométricas"),
       col = c(1, 4),
       pch = c(1, 3),
       lty = c(1, 3),
       lwd = 2)
axis(1,at = time(ytf),
     labels = c("dic-20", "ene-21", "feb-21", "mar-21",
                "abr-21", "may-21", "jun-21", "jul-21",
                "ago-21", "sep-21", "oct-21", "nov-21"))

### ----------------------------------------------------------------------
### 6. Verificación de ruido blanco de los errores estructurales
###    del modelo global
### ----------------------------------------------------------------------

# Gráficas de la ACF y la PACF de los errores estructurales

acf(as.numeric(residuals(modelo_global)),
    ci.type = "ma",
    main = "ACF del modelo global",
    lag.max = 36,
    xlab = "Rezago")
pacf(as.numeric(residuals(modelo_global)),
     main = "PACF del modelo global",
     lag.max = 36,
     xlab = 'Rezago',
     ylab = 'ACF parcial')

# Test de Ljung-Box

BP.LB.test(residuals(modelo_global),
           maxlag=36,
           type="Ljung")

# Test de Dubin-Watson

#pruebaDW1(modelo_global) ## XXXX. ¡PROBLEMAS! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### ----------------------------------------------------------------------
### 7. Identificación de posibles modelos ARMA
### ----------------------------------------------------------------------

# EACF

eacf <- eacf(residuals(modelo_global),ar.max = 36, ma.max = 36)

# Con selectModel() de la librería FitAR para modelos AR(p)

SelectModel(residuals(modelo_global),lag.max=36,Criterion="AIC",ARModel="AR")
SelectModel(residuals(modelo_global),lag.max=36,Criterion="BIC",ARModel="AR")

# Usando la función auto.arima() de la librería forecast.

serie_et=ts(residuals(modelo_global), freq = 12,start = c(2001, 1))              # Serie de tiempo de los residuales
auto.arima(serie_et,ic="aic")
auto.arima(serie_et,ic="bic")

# Usando armasubsets() de la librería TSA

plot(armasubsets(residuals(modelo_global),
                 nar=12,nma=12,
                 y.name='AR',
                 ar.method='ml'))


##-----------------------------------------------------------------
##                             PUNTO 6
##-----------------------------------------------------------------


#Modelo 1 AR(19) no es posible obtener pronosticos por ip


mod1=regexpo.ErrorARMA(respuesta=yt,names.param=param2,data=X1,newdata=X1nuevo,order=c(19,0,0),method="ML") 
ytpron1 = mod1$forecast

#Medidas precisión pronósticos puntuales 
accuracy(ytpron1,ytf) 

#Modelo 2 ARMA(7,11) no es posible obtener pronosticos por ip
mod2=regexpo.ErrorARMA(respuesta=yt,names.param=param2,data=X1,newdata=X1nuevo,order=c(7,0,11),method="ML") 
ytpron2 = mod2$forecast

#Medidas precisión pronósticos puntuales 
accuracy(ytpron2,ytf) 

#modelo 3 ARMA(4,8)(1,1)[12]
mod3=regexpo.ErrorARMA(respuesta=yt,names.param=param2,data=X1,newdata=X1nuevo,order=c(3,0,9),seasonal=list(order=c(1,0,0)), 
                         method="ML") 
ytpron3 = mod3$forecast
para<-mod3$coefficients

#Medidas precisión pronósticos puntuales 
accuracy(ytpron3,ytf) 

#modelo 4 ARMA(3, 12)
mod4 = regexpo.ErrorARMA(respuesta=yt,names.param=param2,data=X1,
                            newdata=X1nuevo,order=c(12,0,10),
                            fixed= c(NA,NA,NA,rep(0,3),NA,rep(0,4),NA,rep(0,3),NA,rep(0,4),NA,NA),
                            method="ML")
ytpron4 = mod4$forecast

#Medidas precisión pronósticos puntuales 
accuracy(ytpron4,ytf)


# ------------------------------------------
# Comparacion grafica de los pronosticos
# ------------------------------------------

win.graph()

vector.auxiliar <- c(ytf, ytpron1, ytpron2,
                    ytpron3 , ytpron4)
par(adj = 0.5)
plot(ytf, type = "b", pch = 19, lty = 1, col = 1, lwd = 2,
     ylab = "Índice de ventas nominales",
     xlab = "Periodo [mmm - yy]",
     ylim = c(min(vector.auxiliar), max(vector.auxiliar)),
     xaxt = "n",
     main = "Comparación de los pronósticos para el periodo ex post")
lines(ytpron1, col = 2, pch = 2, lty = 2, type = "b", lwd = 2)
lines(ytpron2, col = 3, pch = 3, lty = 3, type = "b", lwd = 2)
lines(ytpron3, col = 4, pch = 4, lty = 4, type = "b", lwd = 2)
lines(ytpron4, col = 5, pch = 5, lty = 5, type = "b", lwd = 2)
legend("bottomright",
       legend = c("Real", "Modelo1", "Modelo 2",
                  "Modelo 3", "Modelo 4" ), 
       col = 1:5,
       pch = c(19, 2:5),
       lty = 1:5,
       lwd = 1,
       cex= 0.7)
axis(1,at = time(ytf), 
     labels = c("dic-20", "ene-21", "feb-21", "mar-21",
                "abr-21", "may-21", "jun-21", "jul-21",
                "ago-21", "sep-21", "oct-21", "nov-21"))
###---------------------------------------------------------------
###                             PUNTO 7
###---------------------------------------------------------------
#MEJOR MODELO LOCAL DEL TRABAJO 1 FUE EL MODELO 4 Filtro descomposicion 
#combinado con loess lineal, usando criterio GCV para escoger parametro
#suavizamiento loess


modloc=Descomp.Loess(serie.ajuste=yt,h=m,tipo.descomp="multiplicative",grado=1,criterio="gcv")

#Graficos de la serie y su ajuste final

#LOCAL
plot(datos,lwd=2)
lines(fitted(modloc),col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste por DLL(GCV)"),col=c(1,2),lty=1,lwd=2,cex = 0.6)
Criteriosmodelo=exp.crit.inf.resid(residuales=residuals(modloc),n.par=modloc$p)
#MOD1AR19
plot(datos,lwd=2)
lines(fitted(mod1),col=2,lwd=2)
legend("topleft",legend=c("Original","AR(19)"),col=c(1,2),lty=1,lwd=2,cex = 1)
Criteriosmodelo1=exp.crit.inf.resid(residuales=residuals(mod1),n.par=mod1$p); Criteriosmodelo1


#MOD2ARMA(7,11)

plot(datos,lwd=2)
lines(fitted(mod2),col=2,lwd=2)
legend("topleft",legend=c("Original","ARMA(7,11)"),col=c(1,2),lty=1,lwd=2,cex = 1)
Criteriosmodelo2=exp.crit.inf.resid(residuales=residuals(mod2),n.par=mod2$p); Criteriosmodelo2


#MOD3 ARMA(4,8)(1,1)[12]
plot(datos,lwd=2)
lines(fitted(mod3),col=2,lwd=2)
legend("topleft",legend=c("Original","ARMA(4,8)(1,1)[12]"),col=c(1,2),lty=1,lwd=2,cex = 0.7)
Criteriosmodelo3=exp.crit.inf.resid(residuales=residuals(mod3),n.par=mod3$p); Criteriosmodelo3


#modelo 4 
plot(datos,lwd=2)
lines(fitted(mod4),col=2,lwd=2)
legend("topleft",legend=c("Original","modelo 4"),col=c(1,2),lty=1,lwd=2,cex = 0.7)
Criteriosmodelo4=exp.crit.inf.resid(residuales=residuals(mod4),n.par=mod4$p); Criteriosmodelo4


###----------------------------------------------------------------
###               GRAFICO RESIDUALES
###-----------------------------------------------------------------
#Gráficos de residuales 
plot(residuals(modloc),lwd=2,ylim=c(min(-2*sqrt(modloc$MSE),residuals(modloc)),max(2*sqrt(modloc$MSE),residuals(modloc))))
abline(h=c(-2*sqrt(modloc$MSE),0,2*sqrt(modloc$MSE)),col=2)
legend("topleft",legend="Modelo local: DLL(GCV)",lwd=2,cex=0.5)


plot(as.numeric(fitted(modloc)),residuals(modloc),cex=1.5,ylim=c(min(-2*sqrt(modloc$MSE),residuals(modloc)),max(2*sqrt(modloc$MSE),residuals(modloc))))
abline(h=c(-2*sqrt(modloc$MSE),0,2*sqrt(modloc$MSE)),col=2)
legend("topleft",legend="Modelo local: DLL(GCV)",lwd=2, cex = 0.5)


#ACF sobre residuales de ajuste. Use valor para m el que se indica en la guía del trabajo 

acf(as.numeric(residuals(modloc)),ci.type="ma",lag.max=m,main="ACF modelo local",ci.col=2) 

#PACF sobre residuales de ajuste. Use valor para m el que se indica en la guía del trabajo 

win.graph(width=4.875,height=3.5,pointsize=8) 
pacf(as.numeric(residuals(modloc)),lag.max=m,main="PACF modelo local",ci.col=2) 

BP.LB.test(residuals(modelo),maxlag=m,type="Ljung") #test Ljung-Box use máximo m igual al de ACF y PACF 

#Normalidad sobre residuales de ajuste en el modelo. Sólo si no se rechaza supuesto de ruido blanco 
shapiro.test(residuals(modelo)) 



###pronost
ytpronloc=modloc$ytpron
accuracy(ytpronloc,ytf)
