#TRABAJO R OLGA MIRANDA LOZANO 
#1.Cargamos los datos y examinamos. ¿Cuántas variables hay? ¿Cuántos tratamientos? Hay 2 variables y 5 tratamientos. 

getwd()
list.files()
datos <- read.table("datos-trabajoR.txt", header=TRUE)
datos
head(datos)
summary(datos)
dim(datos)
str(datos)
 
#2. Hacemos un boxplot para cada variable (variable 1 en rosa, variable 2 en azul).
boxplot(Variable1 ~ Tratamiento, data=datos, col = c ("pink"), main = "Boxplot de Variable 1")
boxplot(Variable2 ~ Tratamiento, data=datos, col = c ("blue"), main = "Boxplot de Variable 2")

#3. Hacemos un gráfico de dispersión con las dos variables (cada tratamiento de un color distinto)
Variable1 <- datos$Variable1
Variable2 <- datos$Variable2
Tratamiento <- datos$Tratamiento
plot(Variable1 ~ Variable2)
plot (datos$Variable1, datos$Variable2, xlab = "V1", ylab= "V2", pch= 19, col = (datos$Tratamiento))

#4. Ponemos leyenda al gráfico del apartado anterior. En el margen inferior derecho
legend(x = "bottomright", legend = c("T1", "T2", "T3", "T4", "T5"), fill = c("black","red", "green", "lightblue", "blue"), 
       title = "Tratamientos")

#5. Hacemos un histograma para cada variable. Mantenemos los colores
hist(datos$Variable1, col="pink", main="Histograma Variable1", xlab="Variable1")
hist(datos$Variable2, col="blue", main="Histograma Variable2", xlab="Variable2")

#6.  Hacemos un factor en la columna tratamiento y lo guardamos en una variable
datos$Tratamiento_Factor <- factor(datos$Tratamiento)
summary(datos) #para verificar que se ha creado 

#7. Calculamos la media y la desviación estándar para cada tratamiento
V1.mean = tapply(datos$Variable1, datos$Tratamiento, mean)
V2.mean = tapply(datos$Variable2, datos$Tratamiento, mean)
V1.sd = tapply(datos$Variable1, datos$Tratamiento, sd)
V2.sd = tapply(datos$Variable2, datos$Tratamiento, sd)
head(V1.mean)
head(V2.mean)
head(V1.sd)
head(V2.sd)

#8. Averiguamos cuántos elementos tiene cada tratamiento.
elementos_tratamiento <- table(datos$Tratamiento)
print(elementos_tratamiento) #Cada tratamiento tiene 10 elementos

#9. Extraemos los datos para el tratamiento 1 y el tratamiento 4 y los guardamos cada uno en una variable diferente.
tratamiento1_datos <- datos[datos$Tratamiento == "Tratamiento1", ]
tratamiento4_datos <- datos[datos$Tratamiento == "Tratamiento4", ]
tratamiento1_datos
tratamiento4_datos #Para ver si se ha creado correctamente

#10. Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la Variable 1 son iguales. 
#¿Puedes comprobarlo? Para ello, necesitarás comprobar primero si los datos se distribuyen de forma normal. 
#En función del resultado de la prueba de normalidad, ¿qué test usarías? 
#En general, asumimos que las muestras son independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo a tus resultados.
#Para comprobarlo utilizamos shapiro.tets. 
#La distribución NO es normal ya que nuestros p values son mayores de 0.05 por lo que no podemos recharzar nuestra hipotesis nula. 
shapiro.test(datos$Variable1[datos$Tratamiento_Factor == 1])
shapiro.test(datos$Variable1[datos$Tratamiento_Factor == 4])

#Comprobamos las varianzas 
var.test(tratamiento1_datos$Variable1, tratamiento4_datos$Variable1)

#Realizamos un t-student para comprobar las medias 
#Dado que pvalue es menor de 0.05, decimos que las medias son diferentes y rechazamos la hipótesis nula. 
t.test(Variable1, Tratamiento, var.equal = TRUE)
