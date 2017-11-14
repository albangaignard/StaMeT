
CD C:\Users\f-benazzouz\Desktop\biostats modules

"C:\Users\f-benazzouz\Documents\R\R-3.4.2\bin\R.exe" R --vanilla< Create_directories.R


"C:\Users\f-benazzouz\Documents\R\R-3.4.2\bin\R.exe" R --vanilla< Directories.R


"C:\Users\f-benazzouz\Documents\R\R-3.4.2\bin\R.exe" R --vanilla --args 10000 135 15 0.9 0.5 "FC.txt" < counts_data_simulator.R


"C:\Users\f-benazzouz\Documents\R\R-3.4.2\bin\R.exe"  R --vanilla --args 10000 135 15 0.9 0.5  1.4 0.8 < Array_simulator.R



"C:\Users\f-benazzouz\Documents\R\R-3.4.2\bin\R.exe" R --vanilla --args "simulation" "DESeq" "Zscor" 135 15 135 15 < Cross_plateformes.R