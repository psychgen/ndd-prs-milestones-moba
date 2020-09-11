#grateful
library(grateful)

#1 Scan the project for packages

pkgs <- scan_packages()
#2 Get citations for each package

cites <- get_citations(pkgs, out.dir =paste0(getwd(),"/reports"))
#3 Create an rmarkdown document citing all these packages

rmd <- create_rmd(cites, out.dir =paste0(getwd(),"/reports"))
#4 Rendering the rmarkdown document to the desired output format

render_citations(rmd, output = "html", out.dir =paste0(getwd(),"/reports"))
