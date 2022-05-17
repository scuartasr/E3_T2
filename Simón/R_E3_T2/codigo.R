# Eliminación de todos los objetos anteriores
rm(list=ls(all=TRUE))

### ----------------------------------------------------------------------
### Paquetes y funciones necesarias
### ----------------------------------------------------------------------

# Paquetes

library(forecast)
library(TSA)
library(fANCOVA)

# Funciones de usuario

source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-
Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loess.R")

### ----------------------------------------------------------------------
### Lectura de datos
### ----------------------------------------------------------------------

# Lectura de la base de datos
enlace <- "./anexos-emmet-noviembre-2021-1-total industria-modif.csv"
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
### 3. Planteamiento de modelo exponencial polinomial estacional
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
# Modelo dos
# Exponencial polinomial de grado seis con trigonométricas
# ------------------------------------------

# Crear vector con nombres de los parámetros a usar en la fórmula R 
# del modelo exponencial

param2 <- c(paste0("beta",0:6),
            "alfa1", "gamma1", "alfa2",
            "gamma2", "alfa3", "gamma3",
            "alfa4", "gamma4", "alfa5",
            "gamma5")
param2

modelo2 <- regexponencial(respuesta = yt,
                          data = X1,
                          names.param = param2)
summary(modelo2)

# Cálculo valores ajustados del modelo 2

ythat2 <- ts(fitted(modelo2),
             frequency = 12,
             start = start(yt))

# Pronósticos. Solo son posibles los pronósticos puntuales

predicciones2 <- predict(modelo2,
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
npar2 <- length(coef(modelo2)[coef(modelo2)!=0])
npar2                                 # Número par?metros modelo 2b
Criterios2 <- exp.crit.inf.resid(residuales = residuals(modelo2),
                                 n.par = npar2)
Criterios2

#Gráficos de residuos

plot.ts(residuals(modelo2),
        ylim = c(min(residuals(modelo2),
                     -2*summary(modelo2)$sigma,
                     2*summary(modelo2)$sigma),
                 max(residuals(modelo2),
                     -2*summary(modelo2)$sigma,
                     2*summary(modelo2)$sigma)))
abline(h = c(-2*summary(modelo2)$sigma,
             0,
             2*summary(modelo2)$sigma),
       col = 2)
legend("topleft", legend = c("Modelo 2"),
       lty = 1,
       col = 1,
       lwd = 2)


plot(fitted(modelo2),
     residuals(modelo2),
     ylim = c(min(residuals(modelo2),
                  -2*summary(modelo2)$sigma,
                  2*summary(modelo2)$sigma),
              max(residuals(modelo2),
                  -2*summary(modelo2)$sigma,
                  2*summary(modelo2)$sigma)))
abline(h = c(-2*summary(modelo2)$sigma,
             0,
             2*summary(modelo2)$sigma),
       col = 2)
legend("topleft",
       legend = c("Modelo 2"),
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
### 5. Pronóstico
### ----------------------------------------------------------------------




# ------------------------------------------
# Comparación gráfica de los pronósticos
# ------------------------------------------

vector.auxiliar <- c(ytf, ytpron1, ytpron2,
                     ytpron3, ytpron4)
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
legend("topleft",
       legend = c("Real", "Modelo1", "Modelo 2",
                  "Modelo 3", "Modelo 4"),
       col = 1:5,
       pch = c(19, 2:5),
       lty = 1:5,
       lwd = 2)
axis(1,at = time(ytf), 
     labels = c("dic-20", "ene-21", "feb-21", "mar-21",
                "abr-21", "may-21", "jun-21", "jul-21",
                "ago-21", "sep-21", "oct-21", "nov-21"))