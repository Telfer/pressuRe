% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pressuRe_functions.R
\name{load_pliance}
\alias{load_pliance}
\title{Load pliance data}
\usage{
load_pliance(pressure_filepath)
}
\arguments{
\item{pressure_filepath}{String. Filepath pointing to pliance pressure file}
}
\value{
A list with information about the pressure data.
\itemize{
  \item pressure_array. 3D array covering each timepoint of the measurement.
           z dimension represents time
  \item pressure_system. String defining pressure system
  \item sens_size. String with sensor type
  \item time. Numeric value for time between measurements
  \item masks. List
  \item events. List
  \item sensor_polygons. Data frame with corners of sensors
  \item max_matrix. Matrix
 }
}
\description{
Imports and formats .asc files collected on pliance system and
   exported from Novel software
}
\examples{
pliance_data <- system.file("extdata", "pliance_test.asc", package = "pressuRe")
pressure_data <- load_pliance(pliance_data)
}
