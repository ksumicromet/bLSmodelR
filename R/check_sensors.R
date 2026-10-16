library(bLSmodelR)
library(data.table)

Sensors <- genSensors(
  LineSensor1 = list(x=c(56.554,5.523), y=c(33.988,46.169), z=c(1.50,2.36), d=1.81),
  LineSensor2 = list(x=c(54.936,0.0),   y=c(4.784,0.0),     z=c(1.33,2.06), d=1.90)
)

# Call internal procSensors
pSens <- bLSmodelR:::procSensors(Sensors)

cat('Unique heights (z values):\n')
print(sort(as.numeric(names(pSens$heights))))
cat('Number of unique heights:', length(pSens$heights), '\n')

cat('\nAll expanded point sensors:\n')
print(pSens$Calc.Sensors[, c('Sensor Name','Point Sensor Name','x-Coord (m)','y-Coord (m)','Sensor Height (m)','Node')])

cat('\nPS_list (sub-sensors per line sensor):\n')
print(pSens$PS_list)

