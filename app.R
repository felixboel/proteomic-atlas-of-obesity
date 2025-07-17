# List of required packages
packages <- c("shiny", "shinyTree", "ggplot2", "plotly", "dplyr", "tidyr", "readr")

# Install any missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}
invisible(lapply(packages, install_if_missing))

# Source app file
source("sub_folder/proteomics_explorer.R")

# Run the shiny app
app <- shinyApp(ui = ui, server = server)
runApp(app, launch.browser = TRUE)
