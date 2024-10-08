% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pressuRe_functions.R
\name{create_mask_auto}
\alias{create_mask_auto}
\title{Automatically mask pressure footprint}
\usage{
create_mask_auto(
  pressure_data,
  masking_scheme,
  foot_side = "auto",
  res_value = 10000,
  plot = TRUE,
  template_mask = NULL
)
}
\arguments{
\item{pressure_data}{List. First item is a 3D array covering each timepoint
of the measurement. z dimension represents time}

\item{masking_scheme}{String. "automask_simple", "automask_novel",
"pedar_mask1", "pedar_mask2", "pedar_mask3", "template".
"simple_automask" applies a simple 3 part mask (hindfoot, midfoot, forefoot)
"automask_novel" attempts to apply a 9-part mask (hindfoot, midfoot, mets,
hallux, lesser toes), similar to the standard novel automask
"pedar_mask1" splits the insole into 4 regions using sensel boundaries:
hindfoot, midfoot, forefoot, and toes- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9470545/
"pedar_mask2" splits the insole into 4 regions using percentages:
hindfoot, forefoot, hallux, and lesser toes- https://jfootankleres.biomedcentral.com/articles/10.1186/1757-1146-7-18
"pedar_mask3" splits the foot into 9 regions using sensel boundaries:
 medial hindfoot, lateral hindfoot, medial midfoot, lateral midfoot, MTPJ1,
 MTPJ2-3, MTPJ4-5, hallux, and lesser toes- https://jfootankleres.biomedcentral.com/articles/10.1186/1757-1146-7-20}

\item{foot_side}{String. "RIGHT", "LEFT", or "auto". Auto uses
auto_detect_side function}

\item{res_value}{Numeric. Adjusting this can help if the line between the forefoot and toes
isn't correct. Default is 100000. This line is calculated using a least cost function and this
parameter basically adjusts the resistance of the pressure value for that algorithm}

\item{plot}{Logical. Whether to play the animation}

\item{template_mask}{Data frame. Mask to be used if "template_mask" is selected as the masking scheme}
}
\value{
List. Masks are added to pressure data variable
\itemize{
  \item pressure_array. 3D array covering each timepoint of the measurement.
           z dimension represents time
  \item pressure_system. String defining pressure system
  \item sens_size. Numeric vector with the dimensions of the sensors
  \item time. Numeric value for time between measurements
  \item masks. List
  \item events. List
 }
}
\description{
Automatically creates mask for footprint data
}
\examples{
emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
pressure_data <- load_emed(emed_data)
pressure_data <- create_mask_auto(pressure_data, "automask_novel",
res_value = 100000, foot_side = "auto", plot = FALSE)
}
