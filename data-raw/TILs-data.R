library(usethis)

<<<<<<< HEAD
#Regulon_Activity_Matrix_TILs = readRDS(file="./Regulon_Activity_Matrix_TILs.Rds")
#Annotation_TILs = readRDS(file="./Annotation_TILs.Rds")

#usethis::use_data(Regulon_Activity_Matrix_TILs, overwrite = TRUE)
#usethis::use_data(Annotation_TILs, overwrite = TRUE)
=======
Regulon_Activity_Matrix_TILs = readRDS(file="./Regulon_Activity_Matrix_TILs.Rds")
Annotation_TILs = readRDS(file="./Annotation_TILs.Rds")

usethis::use_data(Regulon_Activity_Matrix_TILs, overwrite = TRUE)
usethis::use_data(Annotation_TILs, overwrite = TRUE)
>>>>>>> c4160b7bd8cc8fffa9885973eb1f6c861b2dc56b
