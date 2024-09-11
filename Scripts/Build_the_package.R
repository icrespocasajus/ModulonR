# creates description and namespace files
usethis::use_description()
usethis::use_namespace()

# this will setup the folders needed for the data and raw-data
usethis::use_data_raw()
#usethis::use_data(data, overwrite = TRUE)# This line should be included in an RScript to generate the data; a second script (data.R) has to be created, 
# executed and kept in data-raw; a new directory named data with the required files is created. Then you have to document your data

# Create R directory
base::dir.create("R")

# creates Package-level documentation so you can run ?nameofpackage
usethis::use_package_doc()

# created README.Rmd for Github landing page
# an .Rbuildignore file gets created
usethis::use_readme_rmd()
# Creating an R package
#https://sahirbhatnagar.com/rpkg/

# creates license file
#usethis::use_mit_license("Sahir Bhatnagar")

# creates news file
usethis::use_news_md()


# sets up testing infrastructure
usethis::use_testthat()


# this will generate documentation for your functions
library(pacman)
pacman::p_load(sinew)
sinew::makeOxyFile("R/ModulonR_functions_wo_comments.R")# This is going to create ModulonR_functions.R, which has to be renamed to ModulonR_functions.R


# Get out and in of the project
# Click install within the Build tab

# See available functions
pacman::p_functions("ModulonR")


#devtools::check()# WARNING! Be aware that using check from RStudio build panel is going to wipe out the doc folder; vignettes may stop working properly
usethis::use_vignette(name = "README")

tools::buildVignettes(dir = ".", tangle=TRUE)
dir.create("inst")
dir.create("inst/doc")
file.copy(dir("vignettes", full.names=TRUE), "inst/doc", overwrite=TRUE)

# Commit and re-build

browseVignettes('ModulonR')
vignette('README')


#Update Readme.Rmd
devtools::build_readme()