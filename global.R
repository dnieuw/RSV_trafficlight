library(data.table)
library(stringr)
library(Biostrings)
library(pwalign)
library(DECIPHER)
library(jsonlite)

library(shiny)
library(shinydashboard)
#library(bslib)
library(shinyWidgets)
library(shinycssloaders)
library(DT)
#Source data
#source("R/placeholder.R")

rsv_nirsevimab_resistance_table_path <- "resistance_study_data.json"

subtype_mapping <- c(
"KY883567.1 Human respiratory syncytial virus isolate HRSV/A/BuenosAires/ARG/001sanger/2015, complete genome" = "RSV_A",
"KY883569.1 Human respiratory syncytial virus isolate HRSV/B/BuenosAires/ARG/002sanger/2015, complete genome" = "RSV_B"
)

# Shared helpers
source("R/resistance_reference_helpers.R")

#Source modules
source("R/modules/welcome_module.R")
source("R/modules/sequence_quality_module.R")
source("R/modules/traffic_light_module.R")
source("R/modules/additional_info_module.R")

options(shiny.maxRequestSize = 20 * 1024^2)
options(shiny.reactlog = TRUE)

# Add shinyvalidate to replace the fasta validation field.
# Add shinywidgets to have better modal windows
# Add shiny conductor to create "tour" like tutorial of the app.