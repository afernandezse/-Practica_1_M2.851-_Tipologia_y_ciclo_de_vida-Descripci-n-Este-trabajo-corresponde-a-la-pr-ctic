rm(list=ls())

setwd("/home/flatline/Documentos/Master_Data_Science/Ciclo_de_Vida_Dato/Practica_2/gdac.broadinstitute.org_BRCA.mRNAseq_Preprocess.Level_3.2016012800.0.0")

# Librerías
library(gdata)
library(tidyr)
library(Boruta)
library(caret)
library(lubridate)
library(car)
library(missForest)
library(VIM)
library(tidyr)
library(purrr)
library(tibble)
library(PerformanceAnalytics)
library(knitr)


# Cargamos las dos fuentes de datos

## Datos transcriptómicos obtenidos de http://firebrowse.org/
tcga_data <- read.table("BRCA.mRNAseq_raw_counts.txt", header = T, sep = "\t",
                        check.names=FALSE)

rownames(tcga_data) <- tcga_data[,1]
tcga_data <- tcga_data[,-1]

tcga_data <- t(tcga_data)


setwd("/home/flatline/Documentos/Master_Data_Science/Ciclo_de_Vida_Dato/Practica_2")

## Datos clinicopatológicos y expresión de 4 proteínas obtenidas de kaggle

clin_data <- read.csv("BRCA.csv",header=T, fill=T)
rownames(clin_data) <- clin_data[,1]
clin_data <- clin_data[,-1]

## Como tenemos más casos en el dataset de expresión génica (tcga_data), aunque no
## son todos comunes, los usamos para la integración

selection <- rownames(clin_data)

################# 2. Integración y selección de los datos de interés a analizar. ##################################################

## Nos quedamos con 239 muestras comunes
## entre el dataset clínico y de proteínas y el dataset transcriptómico

tcga <- tcga_data[rownames(tcga_data) %in% selection,]


## Creación de la variable respuesta categórica dicotómica
## Cálculo de la supervivencia en días

clin_data$OS <- as.numeric(dmy(clin_data$Date_of_Last_Visit) -
  dmy(clin_data$Date_of_Surgery),units="days")


# Codifico la Supervivencia Global en días según sea mayor o igual o menor a la mediana
# Esta será la variable respuesta con etiquetas "High OS" y "Low_OS"

clin_data$c_OS <- ifelse(clin_data$OS < median(clin_data$OS, na.rm=T), "Low OS",
                         "High OS")
 
# Contamos el número de NAs

prop.table(table(is.na(clin_data$OS)))


# Como se trata de un pequeño porcentaje (5.08%) opto por eliminar
## las instancias sin fecha de última revisión del dataset clin_data
## antes de la integración, luego de la integración volverá a limpiarse
## el dataset

clin_data <- clin_data %>% drop_na()

## Quedan 317 muestras


# Integración de los datos de expresión génica con los datos clínicos

data <- merge(clin_data, tcga_data, by=0, all=TRUE)
data <- data %>% drop_na()

## Quedan 227 muestras comunes a ambos datasets

# Selección de genes más informativos

## Creamos dataset temp_data con la expresión de los 20532 genes
## Como es una cantidad muy alta, habrá que hacer una limpieza muy
## drástica


temp_data <- data[,18:length(data)]
rownames(temp_data) <- data$Row.names


### LIMPIEZA DEL DATASET TRANSCRIPTÓMICO #########################

# En primer lugar eliminamos todos los genes que tengan en total menos de 227 lecturas,
# lo que indicará que la información de expresión de ese gen es poco informativa. Este
# cut-off corresponde a una lectura promedio por caso

temp_clean <- temp_data[2:length(temp_data)][,-which(colSums(temp_data[sapply(temp_data, is.numeric)]) < 227)]

# Eliminamos 1365 genes que tienen menos de 227 lecturas

#### Feature selection near to zero variance (by rows)

# nearZeroVar() del paquete caret es una función que identifica aquellas variables
# con una varianza cercana a 0 y que por tanto tampoco tienen poder discriminante
# entre los 2 grupos de la variable respuesta

nzv_metrics <- nearZeroVar(temp_clean,freqCut = 95/5, uniqueCut = 10, saveMetrics = TRUE,
                           names = FALSE, foreach = FALSE, allowParallel = TRUE)

f <- names(temp_clean)[nearZeroVar(temp_clean)]

temp_clean <- temp_clean[,!names(temp_clean) %in% f]
rownames(temp_clean) <- data$Row.names


## 22 genes con Varianza cercana a 0



# ## Convert to integers
# table(sapply(raw, class))
# chars <- sapply(raw, is.numeric)
# raw[ , chars] <- as.data.frame(apply(raw[ , chars], 2, as.integer))


# Ahora eliminamos las variables muy correlacionadas que no aportan información al modelo,
# crean ruido y aumentan el coste computacional, especialmente importante en dataframes
# de dimensionalidad tan grande como este

# En el primer filtro, vamos a eliminar todas las variables con un coeficiente de correlación
# mayor de 0.85

cor_temp_clean <- cor(temp_clean)
highlyCorrelated <- findCorrelation(cor_temp_clean, cutoff=0.85)

# Tenemos 3724 atributos con R > 0.85 que vamos a limpiar del dataset

temp_clean <- temp_clean[, -highlyCorrelated]

# El siguiente criterio es eliminar todos los atributos que tengan una baja desviación std
## Eliminación de la desviación estandár hasta el primer cuartil (227.4476)

sd_vector <- apply(temp_clean, 2, sd)
table(sd_vector)

## Son 3856 variables a eliminar

rm_sd <- sd_vector > quantile(sd_vector)[2]
tmp_clean <- temp_clean[,rm_sd]  # 3866 genes eliminados
temp_clean <- tmp_clean

c_OS <- temp_data$c_OS
temp_clean$c_Os <- c_OS
temp_clean$c_Os <-as.factor(temp_clean$c_Os)

# Selección de atributos informativos con el algoritmo BORUTA
# basado en Random Forest

set.seed(111)
boruta.train <- Boruta(c_Os~., data = temp_clean, doTrace = 2)
print(boruta.train)



## Como quedan 25 atributos en los que no está clara
## Su mayor importancia que su shadow feature, repito
## el proceso con más iteraciones

set.seed(123)
boruta.train2 <- Boruta(c_Os~., data = temp_clean,
                        doTrace = 1, maxRuns=8000)

print(boruta.train2)

select_boruta <- getSelectedAttributes(boruta.train2, withTentative = F)
select_boruta <- gsub("`","",select_boruta)



set.seed(123)
boruta.train3 <- Boruta(c_Os~., data = temp_clean[,c("c_Os",select_boruta)],
                        doTrace = 1, maxRuns=2000)

print(boruta.train3)

gene_exp <- temp_clean[,select_boruta]

## Con 6000 iteraciones retenemos 30 variables, como es un número manejable, seguimos
## con ellas el análisis.


# Ahora creamos el dataframe defiinitivo de trabajo, que incluirá las variables clínicas,la expresión de las 4 proteínas,
# y la expresión de los 30 genes que más poder discriminante respecto a la variable c_Os (supervivencia global), este dataset
# ya es manejable.


def_data <- cbind(data[,1:18],gene_exp)
def_data <- def_data[,-1]

# Antes de seguir limpiamos el workspace de todos los dataframes y variables auxiliares que hemos ido creando

rm(selection, f, highlyCorrelated, nzv_metrics, sd_vector,
   cor_temp_clean, temp_data, tmp_clean, select_boruta,
   temp_clean, gene_exp, tcga, tcga_data, boruta.train,
   boruta.train2)


## Normalización de los datos de expresión génica

num_var <- unlist(lapply(def_data, is.numeric))
def_data_num <- def_data[,num_var]

# Estandarización
process <- preProcess(def_data_num,method=c("center","scale"))
def_data_norm <- predict(process, def_data_num)

def_data_temp <- cbind(def_data[,-which(names(def_data) %in% colnames(def_data_norm))], def_data_norm)
def_data <- def_data_temp

# Eliminamos los campos de fechas que hemos usado para calcular la supervivencia global (OS), porque ya no
# son informativos.
def_data <- def_data[,-c(8,9)]


rm(def_data_norm, def_data_num, def_data_temp)



################## 3. Limpieza de los datos #######################################

# 3.1 Comprobación de 0 y NAs

sum(is.na(def_data)) / length(def_data)
sum(def_data==0) / length(def_data)


# No tenemos valores NA, porque ya hemos limpiado el 5.04% que teníamos antes de la integración
# de los 2 datasets usados en esta práctica. Tampoco tenemos ningún valor 0, aunque es de resaltar
# que en este caso no trataríamos los 0, porque en expresión génica y proteíca, un valor de 0 es informativo
# indicando que las células que componen esa muestra tumoral, no expresan ese gen o proteína. Vamos a 
# introducir otro 2.1% de Nas con la función prodNA

def_data_temp <- prodNA(def_data, noNA=0.005)
sum(is.na(def_data_temp)) / length(def_data_temp)

# # Tenemos un 4.75% de datos perdidos
# A continuación crearemos 2 dataframes, uno con los NA eliminados (def_data1) y el 2 con los valores imputados
# con un algoritmo k-nearest neighbors usando la librería VIM y la función KNN

def_data1 <- def_data_temp %>% drop_na()
sum(is.na(def_data1)) / length(def_data1)


## Quedan 180 casos, aún reduciendo hasta un 0.5% el número de casos que quedan son 180, con lo que perdemos 47 que en 
 # en una variable u otra tengan NAs

df <- data.frame(apply(def_data_temp, 2, as.factor))
def_data2 <- kNN(def_data_temp, k=3)
sum(is.na(def_data2)) / length(def_data2)


# 3.2 Identificación y tratamiento de valores extremos.

# Identificación outliers


outliers <- function(dataframe){
  dataframe %>%
    dplyr::select_if(is.numeric) %>% 
    map(~ boxplot.stats(.x)$out) 
  
  }

outliers(def_data_temp)

## DataFrame con los outliers eliminados (def_data3)


outlierreplacement <- function(dataframe){
  dataframe %>%          
    map_if(is.numeric, ~ replace(.x, .x %in% boxplot.stats(.x)$out, NA)) %>%
    dplyr::bind_cols()
}


def_data3 <- outlierreplacement(def_data) %>% drop_na()

# Si eliminamos todos los outliers, nos quedamos con sólo 53 casos, no es viable




################## 4. Análisis de los datos #######################################
  
## 4.1 en el word Practica2_draw.doc
  
################## 4.2. Comprobación de la normalidad y homogeneidad de la varianza. #######################################
  
# Test de Shapiro-Wilks para medir normalidad
  
res_shapiro <- def_data%>%
  dplyr::select_if(is.numeric)%>%
  sapply(shapiro.test)%>%
  t()%>%
  data.frame()%>%
  dplyr::select(p.value)%>%
  mutate(Is_normally_distributed=p.value>=.05)
  

# Test de Levene
  
responses <- as.matrix(def_data[,10:length(def_data)])
res_levene <- data.frame(var = colnames(responses), p_levene = rep(NA,dim(responses)[2]))
res_levene$p_levene <- apply(responses,2,function(x) {leveneTest(x ~ as.factor(def_data$c_OS))[1,3]})
res_levene$homocedasticity <- ifelse(results$p < 0.05, "heterocedastico",
                                  "homocedastico")

df_shapiro_levene <- as.data.frame(cbind(res_shapiro[,c(1,2)], res_levene[,c(2,3)]))
colnames(df_shapiro_levene)[1] <- "p_shapiro"


######### 4.3 Aplicación de pruebas estadísticas para comparar los grupos de datos. En ######### 
######### función de los datos y el objetivo del estudio, aplicar pruebas de contraste de ######### 
######### hipótesis, correlaciones, regresiones, etc. Aplicar al menos tres métodos de ######### 
######### análisis diferentes. ######### set.seed(66666)

## Contraste de hipótesis, aplicamos el test U de Mann-Whitney (no paramétrico), porque tenemos datos no distribuidos
## normalmente

res_kw <- def_data%>%
  dplyr::select_if(is.numeric)%>%
  sapply(wilcox.test)%>%
  t()%>%
  data.frame()%>%
  dplyr::select(p.value)%>%
  mutate(Different_OS=p.value<.05)


res_chi <- def_data[,-c(1,4,5)]%>%
  dplyr::select(where(is.character)) %>%
  dplyr::summarise_all(dplyr::funs(chisq.test(.,def_data$c_OS)$p.value))


res_cor <- def_data[,-c(1,4,5)]%>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::summarise_all(dplyr::funs(cor.test(.,def_data$OS)$p.value))





# división en series de entrenamiento y validación

# Eliminamos el género que al dividirlo queda como 1 solo factor y nos da error
def_data1 <- def_data1[,-1]
def_data2 <- def_data2[,-1]

def_data1 <- def_data1[,-c(3,4)]
def_data2 <- def_data2[,-c(3,4)]

# Convertimos chr a factores
def_data1[sapply(def_data1, is.character)] <- lapply(def_data1[sapply(def_data1, is.character)], 
                                       as.factor)


n <- nrow(def_data1)
train_rows <- sample(seq(n), size = .8 * n)
def_train1 <- def_data1[ train_rows, ]
def_test1  <- def_data1[-train_rows, ]


n <- nrow(def_data2)
train_rows <- sample(seq(n), size = .8 * n)
def_train2 <- def_data2[ train_rows, ]
def_test2  <- def_data2[-train_rows, ]


## Random Forest

hiperparam <- expand.grid(mtry = c(2, 5, 10),
                          min.node.size = c(2, 3, 4, 5, 10),
                          splitrule = "gini")


set.seed(123)
seeds <- vector(mode = "list", length = 21)
for (i in 1:20) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparam))
}
seeds[[21]] <- sample.int(1000, 1)


# Definición del entrenamiento
control_train <- trainControl(method = "cv", number = 5,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = FALSE, allowParallel = TRUE)





set.seed(123)
rf_1 <- train(
  form = c_OS ~ .,
  data = def_train1,
  method = "ranger",
  tuneGrid = hiperparam,
  metric = "Accuracy",
  trControl = control_train,
  # Número de árboles ajustados
  num.trees = 500,
  importance = "impurity")

confusionMatrix(predict(rf_1, def_test1),
                as.factor(def_test1$c_OS))


set.seed(123)
rf_2 <- train(
  form = c_OS ~ .,
  data = def_train2,
  method = "ranger",
  tuneGrid = hiperparam,
  metric = "Accuracy",
  trControl = control_train,
  # Número de árboles ajustados
  num.trees = 500,
  importance = "impurity")

confusionMatrix(predict(rf_2, def_test2),
                as.factor(def_test2$c_OS))

## Extreme gradient boosting

hiperparam <- expand.grid(nrounds = c(100,200),
                       max_depth = c(10, 15, 20, 25),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       eta = 0.1,
                       gamma=0,
                       min_child_weight = 1,
                       subsample = 1)

set.seed(123)
seeds <- vector(mode = "list", length = 21)
for (i in 1:20) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparam))
}
seeds[[21]] <- sample.int(1000, 1)


set.seed(453) 
xgb_1 = train(
  form = c_OS ~ .,
  data = def_train1,  
  trControl = control_train,
  tuneGrid = hiperparam,
  method = "xgbTree"
)

confusionMatrix(predict(xgb_1, def_test1),
                as.factor(def_test1$c_OS))


set.seed(453) 
xgb_2 = train(
  form = c_OS ~ .,
  data = def_train2,  
  trControl = control_train,
  tuneGrid = hiperparam,
  method = "xgbTree"
)

confusionMatrix(predict(xgb_2, def_test2),
                as.factor(def_test2$c_OS))


########################### 5. Representación de los resultados a partir de tablas y gráficas. #######################

## Histograma de la variable respuesta

hist(table(as.factor(def_data$c_OS)), freq=T,
     xlab = levels(as.factor(def_data$c_OS)),
     ylab = "Frequencies")

## Gráfica de los 30 atributos seleccionados por el algoritmo Boruta

plot(boruta.train3, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta.train3$ImpHistory),function(i)
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
names(lz) <- colnames(boruta.train3$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta.train3$ImpHistory), cex.axis = 0.7)


## Boxplot de los datos normalizados

boxplot(def_data[,10:length(def_data)], pars  =  list(xaxt = "n"), col="darkgreen")
axis(1, at=(10:length(def_data))-9, labels = FALSE)
text((10:length(def_data))-9, par("usr")[3] - 1,
     labels = colnames(def_data[,10:length(def_data)]), srt = 45, pos = 1, xpd = TRUE, cex=0.8)

# Resumen de los datos limpios

summary(def_data_temp)

glimpse(def_data_temp)


## Gráfico de los datos perdidos creados con prodNA

aggr(def_data_temp, numbers=T, sortVars=T, labels=names(def_data_temp),
     cex.axis=.7, gap=3, ylab=c("Missing data", "Pattern"))


## Número de outliers en cada una de las variables
num_outliers <- as.data.frame(lengths(outliers(def_data_temp), use.names = T))
num_outliers <- tibble::rownames_to_column(num_outliers, "variable")
names(num_outliers) <- c("variable", "num outliers")

print(num_outliers)

## Dataframe con los resultados de los tests de normalidad y homocedasticidad

kable(df_shapiro_levene, "rst")

## Correlaciones entre la expresión de las proteínas

chart.Correlation(def_data[,11:14])

## Correlaciones entre la expresión génica

chart.Correlation(def_data[,15:20])
chart.Correlation(def_data[,21:25])
chart.Correlation(def_data[,26:30])
chart.Correlation(def_data[,31:35])
chart.Correlation(def_data[,36:40])
chart.Correlation(def_data[,41:length(def_data)])


# Se ve claramente, que sólo han quedado variables poco correlaciondas, la mayor R es de menos de 0.75


# Resultados de los test de inferencia (U-Mann Whitney)

kable(res_kw, "rst")

# Resultados del test de Chi-cuadrado sobre las variables no numéricas

kable(res_chi, "rst")

# Resultados de la correlación entre la supervivencia global y las variables numéricas de edad y expresión de proteínas y genes

kable(res_cor, "simple")

# Una conclusión que sacamos, es que la limpieza del dataframe ha sido efectiva, casi todos los atributos que quedan,
# correlacionan bien con la variable respuesta.

write.table(def_data, "output.csv",
            append = FALSE, sep="\t", dec = ".", row.names = FALSE, col.names = TRUE)




save.image(file = "PRACTICA2_CVD_afernandezse.RData")
