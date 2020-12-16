library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())

## data location
d <- "~/rds/rds-cew54-wallace-share/Data/elisa"
NFOLDS <- 9
list.files(d)

## pregnant and blood donors
file_xlsx_v1 <- file.path(d, "200606 Chris Donors and Pregnant Cohorts.xlsx")
file_xlsx_v2  <- file.path(d, "200608 BD + Preg All Values.xlsx")
file_xlsx_v3  <- file.path(d, "200609 BD + Preg All Values.xlsx")
file_xlsx_v4 <- file.path(d,"200705 BD Preg Wks 14-25.xlsx")
file_xlsx_v5 <- file.path(d,"200708 S Re-run .xlsx")
file_xlsx_v6 <- file.path(d,"201216_BDPW_Wk45-50.xlsx")

## patient samples
file2_xlsx_v1 <- file.path(d,"200619 Patient Data for Chris - regression.xlsx")
file2_xlsx_v2 <- file.path(d,"200626 Chris, Stas Updated Patient File.xlsx")
file2_xlsx_v3 <- file.path(d,"200627 Updated data table, Chris, Stas.xlsx")
file2_xlsx_v5 <- file.path(d,"200712 Clin Data + Neuts - DUPS REMOVED + PCR.xlsx")

file_rdata_v1 <- file.path(d,"data-200607.RData")
file_rdata_v2 <- file.path(d,"data-200609.RData")
file_rdata_v4 <- file.path(d,"data-200705.RData")
file_rdata_v5 <- file.path(d,"data-200708.RData")
file_rdata_v6 <- file.path(d,"data-201216.RData")
