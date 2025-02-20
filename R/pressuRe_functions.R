# to do (current version)
## low cost paths to define heel and forefoot automasks
## low cost path to define hallux automask
## pliance if sensors matrix then add max_matrix
## fscan processing needs to be checked (work with NA?)
## in edit_mask, make edit_list a vector that works with numbers or names?
## change pedar_polygon to sensel_polygon
## create mask manual work with internal hole
## make capture frequency a vector to allow for uneven sampling

# to do (future)
## global pressure_import function (leave for V2)
## create masks for iscan during startup
## CPEI manual edit to be built into function
## add more input tests to throw errors
## cop for pedar
## UNITS!!!

# data list:
## Array. pressure data
## String. data type (usually collection system, e.g. emed)
## Numeric. sens_size. sensor size
## Numeric. Single number time between samples
## List. Mask list
## Data frame. Events (for example, to define start/end of individual steps for insole data)
## Data frame. Sensor polygon coordinates


# =============================================================================

#' @title Load emed data
#' @description Imports and formats .lst files collected on emed system and
#'   exported from Novel software
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @param trim_active Logical. Restricts frames to only the first continuous
#'   foot contact
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the areas of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix Matrix with maximum image
#'  }
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' @importFrom stringr str_extract_all str_detect str_which
#' @importFrom utils read.fwf read.table
#' @export

load_emed <- function(pressure_filepath, trim_active = FALSE) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "lst")
    stop("incorrect file extension, expected .lst")

  # Read unformatted emed data
  pressure_raw <- readLines(pressure_filepath, warn = FALSE)

  # get sensor size
  sens_size_ln <- which(grepl("Sensor size", pressure_raw, useBytes = TRUE))[1]
  sens_size <- str_extract_all(pressure_raw[sens_size_ln], "\\d+\\.\\d+")
  sens_size <- as.numeric(unlist(sens_size))
  if (str_detect(pressure_raw[sens_size_ln], "cm") == TRUE) {
    sens_size <- sens_size * 0.01
  }

  # get capture frequency
  time_line <- which(grepl("Time", pressure_raw, useBytes = TRUE))[1]
  time_ln <- str_split(pressure_raw[time_line], "picture:")[[1]][2]
  time <- as.numeric(unlist(str_extract_all(time_ln, "\\d+\\.\\d+")))
  if (str_detect(time_ln, "ms") == TRUE) {time <- time / 1000}

  # determine position breaks
  breaks <- str_which(pressure_raw, "Page")

  # identify and remove any summary frames
  frame_type <- pressure_raw[breaks + 8]
  MVP <- which(grepl("MVP", frame_type, fixed = TRUE))
  MPP <- which(grepl("MPP", frame_type, fixed = TRUE))
  if (length(MVP) > 0 | length(MPP) > 0) {
    breaks <- breaks[1:(min(c(MVP, MPP)) - 1)]
  }
  breaks <- breaks[str_detect(frame_type, "Pict")]
  breaks <- breaks[!is.na(breaks)]

  # get blank lines
  ends <- which(pressure_raw == "\x0C")

  # how many frames in each measurement
  nfs <- sum(str_detect(pressure_raw[breaks + 8], "Pict\\-No\\.\\: 1 "),
             na.rm = TRUE)

  # how many measurements
  n_meas <- floor(length(breaks)/nfs)

  # empty array
  pressure_array <- array(0, dim = c(95, 64, n_meas))
  for (i in 1:n_meas) {
    for (j in 1:nfs) {
      # start line
      str <- (nfs * i) - nfs + j

      # load as table
      y <- pressure_raw[(breaks[str] + 10):(ends[which(ends > breaks[str])[1]] - 1)]
      z <- read.table(textConnection(y), sep = "", header = TRUE)

      # remove zeros
      z[is.na(z)] <- 0

      # remove force column
      if (colnames(z)[ncol(z)] == "Force") {z <- z[, 1:(ncol(z) - 1)]}

      # remove force row
      if (rownames(z)[nrow(z)] == "Force") {z <- z[1:(nrow(z) - 1),]}

      # column and row numbers
      cn <- as.numeric(gsub("\\D", "", colnames(z)))
      rn <- as.numeric(gsub("\\D", "", rownames(z)))

      # add to array
      pressure_array[rn, cn, i] <- as.matrix(z)
    }
  }

  # remove zero columns and rows
  fp <- apply(simplify2array(pressure_array), 1:2, max)

  ## rows
  rsums <- rowSums(fp)
  minr <- min(which(rsums > 0))
  maxr <- max(which(rsums > 0))

  ## columns
  csums <- colSums(fp)
  minc <- min(which(csums > 0))
  maxc <- max(which(csums > 0))

  ## update pressure array
  pressure_array <- pressure_array[minr:maxr, minc:maxc,]

  # make max mat
  max_mat <- apply(simplify2array(pressure_array), 1:2, max)

  # active sensor polygons
  active_sensors <- which(max_mat > 0, arr.ind = TRUE)
  sens_array <- sensor_2_polygon3(max_mat, sens_size[1], sens_size[2])

  # active sensor areas
  sens_areas <- sensor_area(sens_array)

  # array to 2d matrix
  full_mat <- matrix(NA, nrow = dim(pressure_array)[3], ncol = nrow(active_sensors))
  for (i in 1:dim(pressure_array)[3]) {
    full_mat[i, ] <- pressure_array[, , i][active_sensors]
  }

  ## reduce to active continuous contact
  if (trim_active == TRUE) {
    # force vector
    fv_ <- full_mat * sens_areas * 1000
    fv <- rowSums(fv_)

    # detect start and end
    str_f <- which(fv > 0)[1]
    end_f_ <- which(fv < 10)
    end_f <- end_f_[which(end_f_ > (str_f + 20))][1]

    # trim array
    full_mat <- full_mat[c(str_f:end_f),]
  }

  # return formatted emed data
  return(list(pressure_array = full_mat, pressure_system = "emed",
              sens_size = sens_areas, time = time, masks = NULL, events = NULL,
              sens_polygons = sens_array, max_matrix = max_mat))
}


# =============================================================================

#' @title Load pedar data
#' @description Imports and formats .asc files collected on pedar system and
#'    exported from Novel software
#' @param pressure_filepath String. Filepath pointing to pedar pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. String with sensor type
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix Matrix with maximum image
#'  }
#' @examples
#' pedar_data <- system.file("extdata", "pedar_example.asc", package = "pressuRe")
#' pressure_data <- load_pedar(pedar_data)
#' @importFrom stringr str_split str_trim
#' @export

load_pedar <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "asc")
    stop("incorrect file extension, expected .asc")

  # measurement data
  ## header
  header <- readLines(pressure_filepath, n = 3)

  # Identify insole
  insole_type_ln <- str_split(header[2], "\\:")[[1]][2]
  insole_type <- str_split(insole_type_ln, "\\-")[[1]][1]
  insole_type <- str_trim(insole_type)[[1]][1]
  if (nchar(insole_type) != 1) {insole_type <- substr(insole_type, 1, 1)}
  insole_type <- tolower(insole_type)

  # Get time between measurements
  time_ln <- str_split(header[3], "\\[Hz\\]\\:")[[1]][2]
  time <- 1 / as.numeric(time_ln)

  # import data
  pressure_array <- as.data.frame(read.table(pressure_filepath,
                                             sep = "", skip = 10,
                                             header = FALSE))
  pressure_array <- as.matrix(pressure_array[, 2:ncol(pressure_array)])

  # pedar sensor polygons
  sens_polygons <- sensor_2_polygon4()

  # pedar sensor areas
  pedar_insole_areas <- pedar_insole_area()
  pedarSensorAreas <- as.vector(pedar_insole_areas[[insole_type]] /
                                  1e06)
  pedarSensorAreas <- c(pedarSensorAreas, pedarSensorAreas)

  # return
  return(list(pressure_array = pressure_array, pressure_system = "pedar",
              sens_size = pedarSensorAreas, time = time, masks = NULL,
              events = NULL, sensor_polygons = sens_polygons, max_matrix = NA))
}


# =============================================================================

#' @title Load pliance data
#' @description Imports and formats .asc files collected on pliance system and
#'    exported from Novel software
#' @param pressure_filepath String. Filepath pointing to pliance pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. String with sensor type
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix. Matrix
#'  }
#' @examples
#' pliance_data <- system.file("extdata", "pliance_test.asc", package = "pressuRe")
#' pressure_data <- load_pliance(pliance_data)
#' @importFrom stringr str_split str_trim
#' @importFrom sf st_centroid
#' @export

load_pliance <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "asc")
    stop("incorrect file extension, expected .asc")

  # global variables
  id <- NULL

  # Read unformatted emed data
  pressure_raw <- readLines(pressure_filepath, warn = FALSE)

  # capture frequency
  freq_line <- which(grepl("scanning rate", pressure_raw, useBytes = TRUE))[1]
  freq <- as.numeric(str_split(pressure_raw[freq_line], "\\[Hz\\]\\: ")[[1]][2])
  time <- 1 / freq

  # pressure array
  breaks <- which(grepl("time\\[", pressure_raw, useBytes = TRUE))
  y <- pressure_raw[(breaks[1] + 1):(breaks[2] - 2)]
  ## remove MVP and MPP lines if present
  MVP <- which(str_detect(y, "MVP"))
  MPP <- which(str_detect(y, "MPP"))
  if (length(MVP) > 0) {y <- y[-(MVP)]}
  if (length(MPP) > 0) {y <- y[-(MPP)]}
  pressure_array_full <- as.matrix(read.table(textConnection(y), sep = ""))
  pressure_array_full <- pressure_array_full[, 2:ncol(pressure_array_full)]
  active_sensors <- which(colSums(pressure_array_full) > 0)
  pressure_array <- pressure_array_full[, active_sensors]

  # sensor coords
  sensor_array_line <- which(grepl("corner", pressure_raw, useBytes = TRUE))[1]
  sensor_array <- pressure_raw[(sensor_array_line + 1):(sensor_array_line + 8)]
  sensor_array <- read.table(textConnection(sensor_array), sep = "")
  sensor_array <- t(sensor_array[, 2:ncol(sensor_array)]) / 100
  sens_array_full <- sensor_2_polygon2(sensor_array)
  sens_array <- sens_array_full %>% filter(id %in% active_sensors)

  # sensor areas
  sens_areas <- sensor_area(sens_array)

  # if sensors are all same size and rectilinear, make max matrix
  if (length(unique(signif(sens_areas, 6))) == 1) {
    # make grid based on sensor sizes
    x_max <- max(sens_array_full$x)
    x_min <- min(sens_array_full$x)
    y_max <- max(sens_array_full$y)
    y_min <- min(sens_array_full$y)
    sens_1 <- sens_array_full %>% filter(id == 1)
    base_x <- max(sens_1$x) - min(sens_1$x)
    base_y <- max(sens_1$y) - min(sens_1$y)
    n_rows <- round(y_max / base_y)
    n_cols <- round(x_max / base_x)

    # sensor centroids
    sens_poly <- sens_df_2_polygon(sens_array_full)
    centroids <- data.frame(x = double(), y = double())
    for (sens in 1:length(sens_poly)) {
      centroids[sens, ] <- st_coordinates(st_centroid(sens_poly[[sens]]))
    }

    # find which grid square each sensor falls in and map to matrix
    pressure_max <- apply(pressure_array_full, 2, max)
    max_mat <- matrix(NA, nrow = n_rows, ncol = n_cols)
    for (i in 1:nrow(centroids)) {
      # col
      mat_col <- ceiling((centroids[i, 1] / x_max) * n_cols)

      # row
      mat_row <- (n_rows + 1) - ceiling((centroids[i, 2] / y_max) * n_rows)

      # map to matrix
      max_mat[mat_row, mat_col] <- pressure_max[i]
    }
  } else {max_mat <- NA}

  # return
  return(list(pressure_array = pressure_array, pressure_system = "pliance",
              sens_size = sens_areas, time = time, masks = NULL,
              events = NULL, sensor_polygons = sens_array,
              max_matrix = max_mat))
}


# =============================================================================

#' @title Load Tekscan data
#' @description Imports and formats files collected on tekscan systems and
#'    exported from Tekscan software
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix. Matrix
#'  }
#' @examples
#' tekscan_data <- system.file("extdata", "fscan_testL.asf", package = "pressuRe")
#' pressure_data <- load_tekscan(tekscan_data)
#' @export

load_tekscan <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  file_ext <- str_split(basename(pressure_filepath), "\\.")[[1]][2]
  if (file_ext != "asf" & file_ext != "csv")
    stop("incorrect file extension, expected .asf or .csv")

  # Read unformatted data
  pressure_raw <- readLines(pressure_filepath, warn = FALSE)

  # get sensor size
  sens_size <- c(NA, NA)
  sens_width <- which(grepl("ROW_SPACING", pressure_raw))
  sens_width_ <- str_extract_all(pressure_raw[sens_width], "\\d+\\.\\d+")
  sens_size[1] <- as.numeric(unlist(sens_width_))
  sens_height <- which(grepl("COL_SPACING", pressure_raw))
  sens_height_ <- str_extract_all(pressure_raw[sens_height], "\\d+\\.\\d+")
  sens_size[2] <- as.numeric(unlist(sens_height_))
  if (str_detect(pressure_raw[sens_width], "mm") == TRUE) {
    sens_size <- sens_size * 0.001
  }
  if (str_detect(pressure_raw[sens_width], "centimeters") == TRUE) {
    sens_size <- sens_size * 0.01
  }

  # get capture frequency
  time_line <- which(grepl("SECONDS_PER_FRAME", pressure_raw))
  time_ln <- str_split(pressure_raw[time_line], "FRAME ")[[1]][2]
  time <- as.numeric(unlist(str_extract_all(time_ln, "\\d+\\.\\d+")))

  # determine position breaks
  breaks <- grep("Frame", pressure_raw)

  # determine matrix dimensions
  col_n <- length(strsplit(pressure_raw[breaks[1] + 1], ",")[[1]])
  row_n <- breaks[2] - breaks[1] - 2

  # empty array
  pressure_array <- array(0, dim = c(row_n, col_n, length(breaks)))

  # fill array
  for (i in 1:length(breaks)) {
    # frame
    y <- pressure_raw[(breaks[i] + 1):(breaks[i] + row_n)]
    z <- read.table(textConnection(y), sep = ",")
    z[z == "B"] <- -1
    z <- as.data.frame(lapply(z, as.numeric))
    pressure_array[, , i] <- as.matrix(z)
  }

  # make max mat
  max_mat <- apply(simplify2array(pressure_array), 1:2, max)

  # active sensor polygons
  active_sensors <- which(max_mat > 0, arr.ind = TRUE)
  sens_polygons <- sensor_2_polygon3(max_mat, sens_size[1], sens_size[2])

  # active sensor areas
  sens_areas <- sensor_area(sens_polygons)

  # array to 2d matrix
  full_mat <- matrix(NA, nrow = dim(pressure_array)[3], ncol = length(active_sensors))
  for (i in 1:dim(pressure_array)[3]) {
    full_mat[i, ] <- pressure_array[, , i][active_sensors]
  }

  # return formatted data
  return(list(pressure_array = pressure_array, pressure_system = "tekscan",
              sens_size = sens_size, time = time, masks = NULL, events = NULL,
              sensor_polygons = sens_polygons, max_matrix = max_mat))
}


# =============================================================================

#' @title Load footscan data
#' @description Imports and formats files collected on footscan systems
#' (formerly RSScan)
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix. Matrix
#'  }
#' @examples
#' footscan_data <- system.file("extdata", "footscan_test.xls", package = "pressuRe")
#' pressure_data <- load_footscan(footscan_data)
#' @importFrom readxl read_excel
#' @export

load_footscan <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "xls")
    stop("incorrect file extension, expected .xls")

  # Read unformatted data
  pressure_raw <- readxl::read_excel(pressure_filepath,
                                     .name_repair = "minimal")

  # sensor size
  sens_size <- c(0.00508, 0.00762)
  sens_area <- sens_size[1] * sens_size[2]

  # determine position breaks
  c1 <- unname(unlist(as.vector(pressure_raw[, 1])))
  breaks <- grep("Frame", c1)

  # get capture frequency
  time <- regmatches(c1[breaks[2]], gregexpr("(?<=\\().*?(?=\\))",
                                             c1[breaks[2]], perl = T))[[1]]
  time <- as.numeric(unlist(str_split(time, " "))[1]) * 0.001

  # determine matrix dimensions
  col_n <- ncol(pressure_raw)
  row_n <- breaks[2] - breaks[1] - 2

  # empty array
  pressure_array <- array(0, dim = c(row_n, col_n, length(breaks)))

  # frame
  for (i in 1:length(breaks)) {
    z <- pressure_raw[(breaks[i] + 1):(breaks[i] + row_n), ]
    z <- as.matrix(as.data.frame(lapply(z, as.numeric)))
    z <- z[nrow(z):1, ]
    pressure_array[, , i] <- z / (sens_area * 1000)
  }

  # make max mat
  max_mat <- apply(simplify2array(pressure_array), 1:2, max)

  # active sensor polygons
  active_sensors <- which(max_mat > 0)
  sens_polygons <- sensor_2_polygon3(max_mat, sens_size[1], sens_size[2])

  # active sensor areas
  sens_areas <- sensor_area(sens_polygons)

  # array to 2d matrix
  full_mat <- matrix(NA, nrow = dim(pressure_array)[3], ncol = length(active_sensors))
  for (i in 1:dim(pressure_array)[3]) {
    full_mat[i, ] <- pressure_array[, , i][active_sensors]
  }

  # return formatted data
  return(list(pressure_array = pressure_array, pressure_system = "footscan",
              sens_size = sens_areas, time = time, masks = NULL, events = NULL,
              sensor_polygon = sens_polygons, max_matrix = max_mat))
}


# =============================================================================

#' @title Load xsensor data
#' @description Imports and formats files collected on xsensor insole systems
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 2D array covering each timepoint of the measurement.
#'            row dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix. Matrix
#'  }
#' @examples
#' xsensor_data <- system.file("extdata", "xsensor_data.csv", package = "pressuRe")
#' pressure_data <- load_xsensor(xsensor_data)
#' @importFrom abind abind
#' @export
load_xsensor <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "csv")
    stop("incorrect file extension, expected .csv")

  # Read unformatted emed data
  pressure_raw <- readLines(pressure_filepath, warn = FALSE)

  # get sensor size
  sens_w_ln <- which(grepl("Sensel Width", pressure_raw, useBytes = TRUE))[1]
  sens_h_ln <- which(grepl("Sensel Height", pressure_raw, useBytes = TRUE))[1]
  sens_w <- str_extract_all(pressure_raw[sens_w_ln], "\\d+\\.\\d+")
  sens_w <- as.numeric(unlist(sens_w))
  sens_h <- str_extract_all(pressure_raw[sens_h_ln], "\\d+\\.\\d+")
  sens_h <- as.numeric(unlist(sens_h))
  sens_size <- sens_w * sens_h
  if (str_detect(pressure_raw[sens_w_ln], "cm") == TRUE) {
    sens_size <- sens_size * 0.0001
  }

  # get capture frequency
  time_ln <- which(grepl("^Time,", pressure_raw, useBytes = TRUE))[c(10, 11)]
  time_1 <- str_split(pressure_raw[time_ln[1]], ",")[[1]][2]
  time_1 <- as.POSIXct(time_1, format = "%H:%M:%OS")
  time_2 <- str_split(pressure_raw[time_ln[2]], ",")[[1]][2]
  time_2 <- as.POSIXct(time_2, format = "%H:%M:%OS")
  time <- as.numeric(difftime(time_2, time_1))

  # right foot
  ## determine position breaks
  breaks <- str_which(pressure_raw, "FRAME")
  rf <- str_which(pressure_raw, "RF")
  lf <- str_which(pressure_raw, "LF")
  sensels <- str_which(pressure_raw, "SENSELS")

  ## matrix size
  n_row <- which(grepl("Rows", pressure_raw))[1]
  n_row <- as.numeric(unlist(str_extract_all(pressure_raw[n_row], "\\d+")))
  n_col <- which(grepl("Columns", pressure_raw))[1]
  n_col <- as.numeric(unlist(str_extract_all(pressure_raw[n_col], "\\d+")))

  ## get frames
  pressure_array_r <- array(0, dim = c(n_row, n_col, length(breaks)))
  pressure_array_l <- pressure_array_r

  for (frm in seq_along(breaks)) {
    # right
    rf_ln <- rf[which(rf > breaks[frm])][1]
    sensel_ln <- sensels[which(sensels > rf_ln)][1]
    y <- pressure_raw[(sensel_ln + 2):(sensel_ln + 2 + (n_row * 2))]
    y <- y[seq(1, (n_row * 2), 2)]
    z <- read.table(textConnection(y), sep = ",")[, 1:n_col]
    pressure_array_r[,, frm] <- as.matrix(z)

    # left
    lf_ln <- lf[which(lf > breaks[frm])][1]
    sensel_ln <- sensels[which(sensels > lf_ln)][1]
    y <- pressure_raw[(sensel_ln + 2):(sensel_ln + 2 + (n_row * 2))]
    y <- y[seq(1, (n_row * 2), 2)]
    z <- read.table(textConnection(y), sep = ",")[, 1:n_col]
    pressure_array_l[,, frm] <- as.matrix(z)
  }

  # join left and right (separate by 2 empty cols)
  empty <- array(0, c(n_row, 2, length(breaks)))
  pressure_array <- abind::abind(pressure_array_l, empty, pressure_array_r, along = 2)

  # make max mat
  max_mat <- apply(simplify2array(pressure_array), 1:2, max)

  # active sensor polygons
  active_sensors <- which(max_mat > 0, arr.ind = TRUE)
  sens_array <- sensor_2_polygon3(max_mat, sens_w, sens_h)

  # active sensor areas
  sens_areas <- sensor_area(sens_array)

  # array to 2d matrix
  full_mat <- matrix(NA, nrow = dim(pressure_array)[3], ncol = nrow(active_sensors))
  for (i in 1:dim(pressure_array)[3]) {
    full_mat[i, ] <- pressure_array[, , i][active_sensors]
  }

  # return formatted xsensor data
  return(list(pressure_array = full_mat, pressure_system = "xsensor",
              sens_size = sens_areas, time = time, masks = NULL, events = NULL,
              sens_polygons = sens_array, max_matrix = max_mat))
}


# =============================================================================

#' @title Interpolate pressure data
#' @description Resamples pressure data over time. Useful for normalizing to
#' stance phase, for example
#' @param pressure_data List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time.
#' @param interp_to Integer. Number of frames to interpolate to
#' @return
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix. Matrix
#'  }
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- pressure_interp(pressure_data, interp_to = 101)
#' @export

pressure_interp <- function(pressure_data, interp_to) {
  # check inputs
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")
  if (is.numeric(interp_to) & interp_to < 2)
    stop("interp_to must be a numeric value greater than 2")
  interp_to <- round(interp_to)

  # make new empty array
  interp_array <- array(NA, dim = c(interp_to, ncol(pressure_data[[1]])))

  # interpolate array
  for (i in 1:ncol(pressure_data[[1]])) {
    interp_array[, i] <- approxP(pressure_data[[1]][, i], interp_to)
  }

  # timing
  time_seq <- seq(0, by = pressure_data[[4]],
                  length.out = nrow(pressure_data[[1]]))
  time_seq_int <- approxP(time_seq, interp_to)
  time_sample_int <- time_seq_int[2] - time_seq_int[1]

  # update pressure data
  pressure_data[[1]] <- interp_array
  pressure_data[[4]] <- time_sample_int

  # return interpolated pressure data
  return(pressure_data)
}


# =============================================================================

#' @title Select steps
#' @description Select steps, usually from insole data, and format for analysis
#' @param pressure_data List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time.
#' @param threshold Numeric. Threshold force to define start and end of step.
#' If "auto", function will set threshold at minimum force in trial + 10N
#' @param min_frames Numeric. Minimum number of frames that need to be in step
#' @param n_steps Numeric. Target number of steps/cycles. User will be
#' asked to keep selected steps until this target is reached or they run out of
#' candidate steps
#' @param skip Numeric. Usually the first few steps of a trial are accelerating
#' and not representative of steady state walking so this removes them
#' @return
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. Data frame
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix. Matrix
#'   }
#' @examplesIf interactive()
#' pedar_data <- system.file("extdata", "pedar_example.asc", package = "pressuRe")
#' pressure_data <- load_pedar(pedar_data)
#' pressure_data <- select_steps(pressure_data)
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab ggtitle
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom utils menu
#' @export

select_steps <- function (pressure_data, threshold = "auto", min_frames = 10,
                          n_steps = 5, skip = 2) {
  # set up global variables
  frame <- NULL

  # check session is interactive
  if (interactive() == FALSE)
    stop("user needs to select suitable steps")

  # check this is pedar (or other suitable) data
  if (!(pressure_data[[2]] == "pedar" || pressure_data[[2]] == "tekscan"))
    stop("data should be from pedar or f-scan")

  # get force df
  if (pressure_data[[2]] == "pedar") {
    te_R <- threshold_event(pressure_data, threshold[1], min_frames, "RIGHT")
    te_L <- threshold_event(pressure_data, threshold[length(threshold)],
                            min_frames, "LEFT")
    df_R <- te_R[[1]]
    df_L <- te_L[[1]]
    FS_events_R <- te_R[[2]]
    FO_events_R <- te_R[[3]]
    FS_events_L <- te_L[[2]]
    FO_events_L <- te_L[[3]]

    # approve steps
    events_R <- approve_step(df_R, FS_events_R, FO_events_R, n_steps, "RIGHT")
    events_L <- approve_step(df_L, FS_events_L, FO_events_L, n_steps, "LEFT")

    # event df
    event_df <- rbind(events_R, events_L)
  } else {
    thr_ev <- threshold_event(pressure_data, threshold[1], min_frames)
    df <- thr_ev[[1]]
    FS_events <- thr_ev[[2]]
    FO_events <- thr_ev[[3]]

    # approve steps
    event_df <- approve_step(df, FS_events, FO_events)
  }

  # add events to pressure data
  pressure_data[[6]] <- event_df

  # Return
  return(pressure_data)
}


# =============================================================================

#' @title Detect foot side
#' @description Detects which foot plantar pressure data is from (left or
#' right), usually would only be needed for barefoot pressure plate data.
#' Generally reliable but may be thrown off by severe deformities or abnormal
#' walking patterns
#' @param pressure_data List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time
#' @return String. "LEFT" or "RIGHT"
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' auto_detect_side(pressure_data)
#' @importFrom sf st_polygon st_as_sf st_area st_filter st_convex_hull
#' st_intersection
#' @importFrom stats dist
#' @export

auto_detect_side <- function(pressure_data) {
  # throw error if pedar data
  if (!(pressure_data[[2]] == "emed" | pressure_data[[2]] == "footscan" |
        pressure_data[[2]] == "pliance"))
    stop("This function currently only works for pressure plate data")

  # Bounding box
  mbb <- getMinBBox(as.matrix(pressure_data[[7]][, c(1, 2)]))
  side1 <- mbb[c(1, 2), ]
  side2 <- mbb[c(3, 4), ]

  # get points
  side1_30 <- side1[1, ] + ((side1[2, ] - side1[1, ]) * 0.3)
  side2_30 <- side2[1, ] + ((side2[2, ] - side2[1, ]) * 0.3)
  side1_60 <- side1[1, ] + ((side1[2, ] - side1[1, ]) * 0.7)
  side2_60 <- side2[1, ] + ((side2[2, ] - side2[1, ]) * 0.7)
  midpt_30 <- (side1_30 + side2_30) / 2
  midpt_60 <- (side1_60 + side2_60) / 2
  midpt_top <- (side1[2, ] + side2[2, ]) / 2

  # make lat and med forefoot boxes
  ff1_pts <- rbind(side1[2, ], midpt_top, midpt_60, side1_60, side1[2, ])
  ff2_pts <- rbind(side2[2, ], midpt_top, midpt_60, side2_60, side2[2, ])
  ffbox1 <- st_polygon(list(ff1_pts))
  ffbox2 <- st_polygon(list(ff2_pts))

  # make lat and med midfoot boxes
  mf1_pts <- rbind(side1_60, midpt_60, midpt_30, side1_30, side1_60)
  mf2_pts <- rbind(side2_60, midpt_60, midpt_30, side2_30, side2_60)
  mfbox1 <- st_polygon(list(mf1_pts))
  mfbox2 <- st_polygon(list(mf2_pts))

  # make chull
  df.sf <- pressure_data[[7]][, c(1, 2)] %>%
    st_as_sf(coords = c( "x", "y" ))
  fp_chull <- st_convex_hull(st_combine(df.sf))

  # count
  #ffbox1_count <- length(st_filter(df.sf, ffbox1)[[1]])
  #ffbox2_count <- length(st_filter(df.sf, ffbox2)[[1]])
  mfbox1_count <- length(st_filter(df.sf, mfbox1)[[1]])
  mfbox2_count <- length(st_filter(df.sf, mfbox2)[[1]])

  # area in each half box
  ffbox1_count <- st_area(st_intersection(fp_chull, ffbox1))
  ffbox2_count <- st_area(st_intersection(fp_chull, ffbox2))

  # ratio
  ff_ratio <- ffbox1_count / ffbox2_count
  mf_ratio <- mfbox2_count / mfbox1_count
  weighted_ratios <- mean(c(ff_ratio, mf_ratio))
  #count_1 <- ffbox1_count
  #count_2 <- ffbox2_count

  # side
  if (weighted_ratios < 1) {side <- "RIGHT"}
  if (weighted_ratios > 1) {side <- "LEFT"}
  #if (count_1 < count_2) {
  #  side <- "RIGHT"
  ##} else {
  #  side <- "LEFT"
  #}

  # Return side
  return(side)
}


# =============================================================================

#' @title Whole pressure curve
#' @description Generates vectors with option to plot for force, peak/mean
#' pressure and area for complete measurement. Useful for checking data
#' @param pressure_data List. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param variable String. "peak_pressure", "force", or "area"
#' @param side For insole data only
#' @param threshold Numeric. Threshold value for sensor to be considered active.
#' Currently only applies to insole data
#' @param plot Logical. If TRUE also plots data as line curve
#' @return Numeric vector containing variable values
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' whole_pressure_curve(pressure_data, variable = "peak_pressure", plot = FALSE)
#' whole_pressure_curve(pressure_data, variable = "area", plot = FALSE)
#' whole_pressure_curve(pressure_data, variable = "force", plot = FALSE)
#' @importFrom ggplot2 aes ggplot geom_line theme_bw xlab ylab
#' @export

whole_pressure_curve <- function(pressure_data, variable, side, threshold = 10,
                                 plot = FALSE) {
  # set global variables
  time <- value <- NULL

  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")
  if(pressure_data[[2]] == "pedar" & missing(side) == TRUE)
    stop("pedar data needs to have side defined")

  # force
  if (variable == "force") {
    if (pressure_data[[2]] == "pedar") {
      if (side == "RIGHT") {
        values <- force_pedar(pressure_data, variable = "force")[, 1]
      }
      if (side == "LEFT") {
        values <- force_pedar(pressure_data, variable = "force")[, 2]
      }
    } else {
      force_array <- pressure_data[[1]] * pressure_data[[3]] * 1000

      # find total force for each frame and store in vector
      values <- rowSums(force_array)
    }
    variable_units <- "force (N)"
  }

  # peak pressure
  if (variable == "peak_pressure") {
    if (pressure_data[[2]] == "pedar") {
      if (side == "RIGHT") {
        mat_r <- pressure_data[[1]][, 100:198]
        values <- apply(mat_r, 1, max, na.rm = TRUE)
      }
      if (side == "LEFT") {
        mat_r <- pressure_data[[1]][, 1:99]
        values <- apply(mat_r, 1, max, na.rm = TRUE)
      }
    } else {
      values <- apply(pressure_data[[1]], 1, max, na.rm = TRUE)
    }
    variable_units <- "peak pressure (kPa)"
  }

  # area
  if (variable == "area") {
    if (pressure_data[[2]] == "pedar") {
      if (side == "RIGHT") {
        values <- force_pedar(pressure_data, variable = "area", threshold)[, 1]
      }
      if (side == "LEFT") {
        values <- force_pedar(pressure_data, variable = "area", threshold)[, 2]
      }
    } else {
      # find active area for each frame and store in vector
      active_cells <- pressure_data[[1]]
      active_cells[active_cells > 0] <- 1
      active_cells <- t(t(active_cells) * pressure_data[[3]])
      values <- rowSums(active_cells) * 10000
    }
    variable_units <- "contact area (cm2)"
  }

  # plot, if required
  if (plot == TRUE) {
    # make df
    variable_df <- data.frame(value = values,
                              time = seq(0, by = pressure_data[[4]],
                                         length.out = length(values)))

    # plot
    g <- ggplot(variable_df, aes(x = time, y = value))
    g <- g + geom_line()
    g <- g + theme_bw()
    g <- g + xlab("time (s)") + ylab(variable_units)
    print(g)
  }

  # return
  return(values)
}


# =============================================================================

#' @title Center of pressure
#' @description Generates xy coordinates for center of pressure during each
#' frame of measurement
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement. z dimension represents time
#' @return Data frame with x and y coordinates of COP throughout trial
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' cop(pressure_data)
#' @export

cop <- function(pressure_data) {
  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")

  # check not pedar
  if (pressure_data[[2]] == "pedar")
    stop("pedar data currently not supported for this function")

  # centroids
  centroids <- sensor_centroid(pressure_data)

  # COP coordinates in x direction
  x_coord <- c()
  for (i in 1:nrow(pressure_data[[1]])) {
    p_total <- sum(pressure_data[[1]][i, ])
    x_coord[i] <- (sum(centroids$x * pressure_data[[1]][i, ])) / p_total
  }

  # COP coordinates in y direction
  y_coord <- c()
  for (i in 1:nrow(pressure_data[[1]])) {
    p_total <- sum(pressure_data[[1]][i, ])
    y_coord[i] <- (sum(centroids$y * pressure_data[[1]][i, ])) / p_total
  }

  # combine coordinates into dataframe
  COP_df <- data.frame(x_coord, y_coord)

  # return COP coordinates
  return(COP_df)
}


# =============================================================================

#' @title Footprint
#' @description Determines footprint of pressure data
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param variable String. "max" = maximum value of each sensor across full
#' dataset. "mean" = average value of sensors over full dataset."frame" = an
#' individual pressure frame. "meanmax" average max values across cycles (
#' currently just for pedar)
#' @param frame Integer. Only used if variable = "frame".
#' @param plot Logical. Display pressure image
#' @return Matrix. Maximum or mean values for all sensors
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' footprint(pressure_data, plot = FALSE)
#' @export

footprint <- function(pressure_data, variable = "max", frame = NULL,
                      plot = FALSE) {
  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")
  stopifnot("Unknown Variable" = variable %in% c("max", "mean", "frame",
                                                 "meanmax"))
  if (variable == "frame" & missing(frame))
    stop("For a single frame, a value is needed for frame variable")

  # calculate footprint for different variables
  if (variable == "max") {
    mat <- apply(pressure_data[[1]], 2, max)
  }
  if (variable == "mean") {
    mat <- apply(pressure_data[[1]], 2, mean)
  }
  if (variable == "frame") {
    if (frame <= nrow(pressure_data[[1]])) {
      mat <- pressure_data[[1]][frame, ]
    } else {
      stop("The frame selected is greater that the number of frames in trial")
    }
  }
  if (variable == "meanmax") {
    if (is.data.frame(pressure_data[[6]]) == FALSE) {
      stop("data must have events defined")
    } else {
      if (pressure_data[[2]] == "pedar") {
        evs <- pressure_data[[6]]
        ev_rows_r <- which(evs[, 1] == "RIGHT")
        ev_rows_l <- which(evs[, 1] == "LEFT")
        mat_peak_r <- matrix(NA, nrow = length(ev_rows_r), ncol = 99)
        mat_peak_l <- matrix(NA, nrow = length(ev_rows_l), ncol = 99)

        # right
        for (i in 1:length(ev_rows_r)) {
          start_frame <- evs[ev_rows_r[i], 2]
          end_frame <- evs[ev_rows_r[i], 3]
          pd <- pressure_data[[1]][start_frame:end_frame, 100:198]
          mat_peak_r[i, ] <- apply(pd, 2, "max")
        }
        mean_peak_r <- apply(mat_peak_r, 2, "mean")
        # left
        for (i in 1:length(ev_rows_l)) {
          start_frame <- evs[ev_rows_l[i], 2]
          end_frame <- evs[ev_rows_l[i], 3]
          pd <- pressure_data[[1]][start_frame:end_frame, 1:99]
          mat_peak_l[i, ] <- apply(pd, 2, "max")
        }
        mean_peak_l <- apply(mat_peak_l, 2, "mean")
        mat <- c(mean_peak_l, mean_peak_r)
      }
    }
  }

  # plot if requested
  if (plot == TRUE) {g <- plot_pressure(pressure_data, variable, frame)}

  # return footprint
  return(mat)
}


# =============================================================================

#' @title Plot pressure
#' @description Produces visualization of pressure data
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param variable String. "max" = footprint of maximum sensors. "mean" =
#' average value of sensors over time (usually for static analyses). "frame" =
#' an individual frame
#' @param frame Integer.
#' @param step_n If numeric, the step number to plot (only for insole data). If
#' "max", the max across complete trial, if "meanmax", the max on a per step
#' basis
#' @param smooth Logical. Not implemented. If TRUE, plot will interpolate
#'   between sensors to increase data density
#' @param plot_COP Logical. If TRUE, overlay COP data on plot. Default = FALSE
#' @param plot_outline Logical. If TRUE, overlay convex hull outline on plot
#' @param plot_masks Logical. If TRUE, overlay mask outline on plot
#' @param plot_colors String. "default": novel color scheme; "custom": user
#' supplied
#' @param break_values Vector. If plot_colors is "custom", values to split
#' colors at
#' @param break_colors Vector. If plot_colors is "custom", colors to use.
#' Should be one shorter than break_values
#' @param sensor_outline Logical. Sensor outline to be shown
#' @param plot Logical. If TRUE, plot will be displayed
#' @param legend Logical. If TRUE, legend will be added to plot
#' @return ggplot plot object
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' plot_pressure(pressure_data, variable = "max", plot_COP = FALSE)
#' plot_pressure(pressure_data, variable = "frame", frame = 20,
#'               plot_colors = "custom", break_values = c(100, 200, 300),
#'               break_colors = c("light blue", "light green", "yellow", "pink"))
#' @importFrom ggplot2 ggplot aes geom_raster geom_polygon scale_fill_manual
#' theme geom_point element_rect binned_scale unit
#' @importFrom scales manual_pal
#' @export

plot_pressure <- function(pressure_data, variable = "max", smooth = FALSE,
                          frame, step_n = "max", plot_COP = FALSE,
                          plot_outline = FALSE, plot_masks = FALSE,
                          plot_colors = "default", break_values, break_colors,
                          sensor_outline = TRUE, plot = TRUE, legend = TRUE) {
  # set global variables
  x <- y <- X <- Y <- id <- cols <- x_coord <- y_coord <- value <- NULL

  # check inputs
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure data must contain array")

  # if pedar
  if (pressure_data[[2]] == "pedar") {
    if (step_n == "max") {pedar_var <- "max"; step_no <- NA}
    if (step_n == "meanmax") {pedar_var <- "meanmax"; step_no <- NA}
    if (is.numeric(step_n)) {pedar_var <- "step_max"; step_no <- step_n}
    cor <- plot_pedar(pressure_data, pedar_var, step_no, foot_side = "both")
    x_lims <- c(0, max(cor$x))
    y_lims <- c(0, max(cor$y))
    legend_spacing <- (cor$x[2] - cor$x[1]) * 10
  } else {
    # max size of df
    x_lims <- c(-0.005, max(pressure_data[[7]]$x) + 0.005)
    y_lims <- c(-0.005, max(pressure_data[[7]]$y) + 0.005)

    # footprint
    fp <- footprint(pressure_data, variable = variable, frame)

    # combine with pressure values
    ids <- unique(pressure_data[[7]]$id)#c(1:length(fp))
    vals <- data.frame(id = ids, value = fp)

    # merge value and coordinate frames
    cor <- merge(pressure_data[[7]], vals, by = c("id"))

    # add colors
    cor <- generate_colors(cor, col_type = plot_colors, break_values,
                           break_colors)

    # plot
    legend_spacing <- 0.005 * 100
  }

  if (plot_colors == "default") {
    break_colors <- c("grey","lightblue", "darkblue","green","yellow",
                      "red", "pink")
    break_values <- c(0, 40, 60, 100, 150, 220, 300)
  }

  # legend range
  range_max <- max(pressure_data[[1]])

  # plot
  g <- ggplot()
  if (sensor_outline == TRUE) {
    g <- g + geom_polygon(data = cor,
                          aes(x = x, y = y, group = id, fill = value),
                          color = "black")
  } else {
    g <- g + geom_polygon(data = cor,
                          aes(x = x, y = y, group = id, fill = value))
  }

  g <- g + binned_scale("fill", "foo",
                        binned_pal(manual_pal(break_colors)),
                        guide = "coloursteps", breaks = break_values,
                        #limits = c(0, max(cor$value)), show.limits = FALSE,
                        limits = c(0, range_max), show.limits = FALSE,
                        name = "Pressure (kPa)")
  g <- g + scale_x_continuous(expand = c(0, 0), limits = x_lims)
  g <- g + scale_y_continuous(expand = c(0, 0), limits = y_lims)
  g <- g + coord_fixed()

  # add COP?
  if (plot_COP == TRUE) {
    cop_df <- cop(pressure_data)
    g <- g + geom_point(data = cop_df, aes(x = x_coord, y = y_coord))
  }

  # add outline
  if (plot_outline == TRUE) {
    ch_out <- st_coordinates(pressure_outline(pressure_data))
    ch_out <- as.data.frame(ch_out)[, c(1, 2)]
    g <- g + geom_path(data = ch_out, aes(x = X, y = Y), colour = "black")
    g <- g + geom_point(data = ch_out, aes(x = X, y = Y), colour = "purple")
  }

  # add masks
  if (plot_masks == TRUE) {
    masks <- pressure_data[[5]]
    if (length(masks != 0)) {
      for (n_mask in 1:length(masks)) {
        mask_plt_df <- data.frame(st_coordinates(masks[[n_mask]])[,1:2])
        g <- g + geom_path(data = mask_plt_df, aes(X, Y), color = "red",
                           linewidth = 1)
      }
    }
  }

  # formatting
  g <- g + theme_void()
  if (legend == FALSE) {
    g <- g + theme(legend.position = "none")
  } else {
    g <- g + theme(panel.background = element_rect(fill = "white",
                                                   colour = "white"),
                   legend.box.spacing = unit(legend_spacing, "cm"))
  }

  # display plot if requested
  if (plot == TRUE) {print(g)}

  # return ggplot object
  return(g)
}


# =============================================================================

#' @title Animate pressure
#' @description Produces animation (gif) of pressure data
#' @param pressure_data Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param plot_colors String
#' @param fps Numeric. Number of frames per second in animation
#' @param dpi Numeric. Resolution of gif
#' @param file_name Name (inlcuding path) of export file
#' @return Animation in gif format
#' @examplesIf interactive()
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' animate_pressure(pressure_data, fps = 10, file_name = "pli_gif.gif")
#' @importFrom stringr str_ends str_pad
#' @importFrom magick image_graph image_animate image_write image_info
#' image_read
#' @importFrom ggplot2 ggplot aes geom_polygon scale_x_continuous
#' scale_y_continuous coord_fixed theme_void xlim ylim ggsave
#' @importFrom grDevices dev.off
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

animate_pressure <- function(pressure_data, plot_colors = "default", fps,
                             dpi = 96, file_name) {
  # parameter check
  if (str_ends(file_name, ".gif") == FALSE)
    stop("filename must end in .gif")

  # max size of df
  x_lim <- max(pressure_data[[7]]$x)
  y_lim <- max(pressure_data[[7]]$y)

  # set up temp directory
  temp_dir <- tempdir()

  # number of frames
  n_frames <- nrow(pressure_data[[1]])

  # plot
  img_fns <- rep(NA, length.out = n_frames)
  pb <- txtProgressBar(min = 1, max = n_frames, style = 3)
  message("processing images")
  for (i in 1:n_frames) {
    # make plot
    g <- plot_pressure(pressure_data, variable = "frame", frame = i,
                       plot_colors = plot_colors, plot = FALSE)
    img_fn <- paste0(temp_dir, "/img", str_pad(i, nchar(n_frames),
                                               pad = "0"), ".png")
    ggsave(img_fn, g, width = 6.45, height = 5.77, dpi = dpi)
    img_fns[i] <- img_fn
    setTxtProgressBar(pb, i)
  }

  # update progress
  close(pb)
  message("generating and saving animation")

  # load images back in
  allInfo <- image_info(image_read(img_fns))
  images <- image_read(img_fns)

  # create animation
  animation <- magick::image_animate(images, fps = fps)#, optimize = TRUE)

  # save animation
  magick::image_write(animation, file_name)
  file.remove(img_fns)
  invisible(gc())
}


# =============================================================================

#' @title Create masking
#' @description Allows user to manually define mask regions
#' @param pressure_data List. First item is a matrix covering each timepoint
#' of the measurement.
#' @param mask_definition String. "by_vertices" or "by_sensors". The first
#' option let's you draw a shape around the area you want to select, the
#' second allows you to define this area by clicking on specific sensors
#' @param n_masks Numeric. Number of masks to add
#' @param n_verts Numeric. Number of vertices in mask
#' @param n_sens Numeric. Number of sensors mask will contain
#' @param threshold Numeric. Distance between adjacent mask vertices before
#' sharing vertex coordinates
#' @param plot_existing_masks Logical. Show existing masks
#' @param mask_names List. Mask names. Default is "custom_mask#"
#' @param plot Logical. Show new maks on pressure image
#' @return List. Mask(s) are added to pressure data variable
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix Matrix with maximum image
#'  }
#' @examplesIf interactive()
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- create_mask_manual(pressure_data, mask_definition = "by_vertices",
#' n_masks = 1, n_verts = 4)
#' pressure_data <- create_mask_manual(pressure_data, mask_definition = "by_sensors",
#' n_masks = 1, n_sens = 8)
#' @importFrom grDevices x11
#' @importFrom ggmap gglocator
#' @importFrom ggplot2 aes geom_path
#' @importFrom sf st_polygon
#' @export

create_mask_manual <- function(pressure_data, mask_definition = "by_vertices", n_masks = 1,
                               n_verts = 4, n_sens = 4, threshold = 0.005,
                               plot_existing_masks = TRUE,
                               mask_names = "default", plot = TRUE) {
  # global variables
  x <- y <- X <- Y <- sensor_polygon <- polygon_list <- NULL

  # check session is interactive
  if (interactive() == FALSE) {stop("user needs to select mask vertices")}

  # get existing masks (if any)
  mask_list <- pressure_data[[5]]
  n_exist_mask <- length(mask_list)

  # empty data frame to store vertex coordinates
  mask_vertices <- data.frame(x = double(), y = double())

  # plot existing masks or just footprint
  if (plot_existing_masks == TRUE & n_exist_mask > 0) {
    g <- plot_masks(pressure_data)
    for (mask_idx in 1:n_exist_mask){
      mask_coord <- st_coordinates(pressure_data[[5]][[mask_idx]])[, 1:2]
      mask_vertices[(nrow(mask_vertices) + 1):(nrow(mask_vertices) + nrow(mask_coord)), ] <-
        mask_coord
    }
  } else {
    grDevices::x11()
    g <- plot_pressure(pressure_data, plot = FALSE)
    #g <- g + scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 0.15))
    #g <- g + scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 0.30))
    print(g)
  }

  # define mask by vertices
  if (mask_definition == "by_vertices") {
    for (mask_n in 1:n_masks){
      if (n_verts < 2){stop("Masks must contain at least 3 vertices")}

      # interactively select area
      message("Select mask corners")
      mask <- gglocator(n_verts)

      # if requested, allow closely selected mask vertices to share a vertex
      if (threshold > 0) {
        mask_vertices[(nrow(mask_vertices) + 1):
                        (nrow(mask_vertices) + n_verts), ] <- mask}
      if (nrow(mask_vertices) > n_verts) {
        v_distance <- as.matrix(dist(mask_vertices))
        for (change_v in 0:(n_verts - 1)) {
          close_v <- which((v_distance[nrow(mask_vertices) - change_v, ] < threshold) &
                             (v_distance[nrow(mask_vertices) - change_v, ] > 0), arr.ind = TRUE)
          if (length(close_v) > 0) {
            mask[n_verts - change_v, ] <- mask_vertices[close_v[1], ]
          }
        }
      }

      # format mask points
      mask <- data.frame(x = mask[, 1], y = mask[, 2])
      mask <- mask[c(1:nrow(mask), 1), ]

      # preview mask
      if (plot == TRUE) {
        g <- g + geom_path(data = mask, aes(x, y), color = "blue", linewidth = 1)
        g <- g + geom_point(data = mask, aes(x, y), color = "blue", size = 2)
        print(g)
      }

      # mask to polygon
      mask_pol <- st_polygon(list(as.matrix(mask)))

      # update mask list
      mask_list[[length(mask_list) + 1]] <- mask_pol
    }
  }

  # define mask by selecting sensors
  if (mask_definition == "by_sensors") {
    # check just one mask
    if (n_masks != 1)
      stop("only one mask can be created when using by_sensors")

    # selection of sensors
    message(paste0("Select ", n_sens, " adjacent sensors to define your custom mask"))
    sensor_pts <- gglocator(n_sens)

    # find which sensors points fall in
    sensor_polygons <- sens_df_2_polygon(pressure_data[[7]])
    sensor_list <- c()
    for (pts in 1:n_sens) {
      point <- sf::st_point(c(sensor_pts[pts, 1], sensor_pts[pts, 2]))
      for (sens_idx in 1:length(sensor_polygons)) {
        if (length(st_intersects(sensor_polygons[[sens_idx]], point)[[1]]) == 1) {
          sensor_list <- c(sensor_list, sens_idx)
        }
      }
    }

    # build custom mask
    sensel_polygon <- st_union(st_sfc(sensor_polygons[sensor_list]))
    #if (length(st_geometry(sensel_polygon)[[1]]) > 1)
    #  stop("sensels must share a corner or an edge")
    if ("L3" %in% colnames(st_coordinates(sensel_polygon)))
      stop("sensels must share a corner or an edge")

    # plot
    if (plot == TRUE) {
      mask <- data.frame(st_coordinates(sensel_polygon)[, 1:2])
      g <- g + geom_path(data = mask, aes(X, Y), color = "blue", linewidth = 1)
      g <- g + geom_point(data = mask, aes(X, Y), color = "blue", size = 2)
      print(g)
    }

    # update mask list
    mask_list <- pressure_data[[5]]
    mask_list[[length(mask_list) + 1]] <- sensel_polygon
  }

  # mask list to main data
  pressure_data[[5]] <- mask_list

  # add labels
  if (mask_names[1] == "default" & length(mask_names) == 1){
    names(pressure_data[[5]])[(n_exist_mask + 1):(n_masks + n_exist_mask)] <- sprintf("custom_mask%d", seq(1, n_masks))
  } else {
    names(pressure_data[[5]])[(n_exist_mask + 1):(n_masks + n_exist_mask)] <- mask_names
  }

  # return
  return(pressure_data)
}


# =============================================================================

#' @title Automatically mask pressure footprint
#' @description Automatically creates mask for footprint data
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement. z dimension represents time
#' @param masking_scheme String. "automask_simple", "automask_novel",
#' "pedar_mask1", "pedar_mask2", "pedar_mask3", "template".
#' "simple_automask" applies a simple 3 part mask (hindfoot, midfoot, forefoot)
#' "automask_novel" attempts to apply a 9-part mask (hindfoot, midfoot, mets,
#' hallux, lesser toes), similar to the standard novel automask
#' "pedar_mask1" splits the insole into 4 regions using sensel boundaries:
#' hindfoot, midfoot, forefoot, and toes-
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9470545/
#' "pedar_mask2" splits the insole into 4 regions using percentages:
#' hindfoot, forefoot, hallux, and lesser toes-
#' https://jfootankleres.biomedcentral.com/articles/10.1186/1757-1146-7-18
#' "pedar_mask3" splits the foot into 9 regions using sensel boundaries:
#'  medial hindfoot, lateral hindfoot, medial midfoot, lateral midfoot, MTPJ1,
#'  MTPJ2-3, MTPJ4-5, hallux, and lesser toes-
#'  https://jfootankleres.biomedcentral.com/articles/10.1186/1757-1146-7-20
#' @param foot_side String. "RIGHT", "LEFT", or "auto". Auto uses
#' auto_detect_side function
#' @param res_value Numeric vector. Adjusting these values can help if the heel,
#' midfoot, toe, and hallux lines aren't correct. Default values are c(10000,
#' 10000, 100000, 10000). These lines are calculated using a least cost function
#' and the parameter essentially adjusts the resistance of the pressure value
#' for that algorithm
#' @param plot Logical. Whether to play the animation
#' @param template_mask Data frame. Mask to be used if "template_mask" is
#' selected as the masking scheme
#' @return List. Masks are added to pressure data variable
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- create_mask_auto(pressure_data, "automask_novel",
#' res_value = c(20000, 20000, 100000, 20000), foot_side = "auto", plot = FALSE)
#' @importFrom zoo rollapply
#' @importFrom sf st_union st_difference
#' @export

create_mask_auto <- function(pressure_data, masking_scheme, foot_side = "auto",
                             res_value = c(20000, 20000, 100000, 20000),
                             plot = TRUE, template_mask = NULL) {
  # simple
  if (masking_scheme == "automask_simple") {
    if (pressure_data[[2]] == "pedar")
      stop("automask is not compatible with this type of data")
    pressure_data <- automask(pressure_data, res_scale = res_value,
                              "automask_simple", foot_side = foot_side,
                              plot = FALSE)
  }

  # full auto mask novel
  if (masking_scheme == "automask_novel") {
    if (!(pressure_data[[2]] == "emed" || pressure_data[[2]] == "pliance"))
      stop("automask is not compatible with this type of data")
    pressure_data <- automask(pressure_data, "automask_novel",
                              res_scale = res_value, foot_side = foot_side,
                              plot = FALSE)
  }

  # pedar masks
  ## pedar mask 1
  if (masking_scheme == "pedar_mask1") {
    if (pressure_data[[2]] != "pedar")
      stop("pedar_mask1 is not compatible with this type of data")
    pressure_data[[5]] <- pedar_mask1(pressure_data)
  }

  ## pedar mask 2
  if (masking_scheme == "pedar_mask2") {
    if (pressure_data[[2]] != "pedar")
      stop("pedar_mask2 is not compatible with this type of data")
    pressure_data[[5]] <- pedar_mask2(pressure_data)
  }

  ## pedar mask 3
  if (masking_scheme == "pedar_mask3") {
    if (pressure_data[[2]] != "pedar")
      stop("pedar_mask3 is not compatible with this type of data")
    pressure_data[[5]] <- pedar_mask3(pressure_data)
  }

  # apply a template mask
  if (masking_scheme == "template") {
    pressure_data <- align_mask(pressure_data, template_mask)
  }

  # plot masks
  if (plot == TRUE) {plot_masks(pressure_data)}

  # return
  return(pressure_data)
}


# =============================================================================

#' @title Edit mask
#' @description Allows user to manually adjust mask vertices
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param n_edit Numeric. Number of vertices to edit
#' @param threshold Numeric. Distance between point clicked and vertex that is
#' selected
#' @param edit_list List. Mask numbers that want to be edited. (Default is to
#' load all masks so that adjacent masks with shared coordinates are modified
#' together)
#' @param image String."max" = footprint of maximum sensors. "mean"
#' average value of sensors over time
#' @return List. Edited mask is added to the pressure data variable
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   \item sensor_polygons. Data frame with corners of sensors
#'   \item max_matrix Matrix with maximum image
#'  }
#' @examplesIf interactive()
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- create_mask_auto(pressure_data, "automask_novel",
#' foot_side = "auto", plot = FALSE)
#' pressure_data <- edit_mask(pressure_data, n_edit = 1, threshold = 0.002,
#' image = "max")
#' @importFrom grDevices rainbow
#' @export

edit_mask <- function(pressure_data, n_edit, threshold = 0.002,
                      edit_list = seq(1, length(pressure_data[[5]])),
                      image = "max") {
  # check session is interactive
  if (interactive() == FALSE) {stop("user needs to select mask vertices")}

  # global variables
  X <- Y <- NULL

  if (n_edit > 0) {
    # initialize plot
    grDevices::x11()

    # plot original mask data
    g <- plot_masks(pressure_data, image = image)

    # compile all mask vertices
    mask_org <- pressure_data[[5]]
    mask_vertices_edit <- data.frame(x = double(), y = double())

    for (n_mask in edit_list) {
      mask_vertices <- data.frame(st_coordinates(mask_org[[n_mask]])[,1:2])
      mask_vertices_edit[(nrow(mask_vertices_edit) + 1):
                           (nrow(mask_vertices) + nrow(mask_vertices_edit)), ] <- mask_vertices
    }

    # select point near edit vertex and select new value
    color_v <- rainbow(n_edit + 1)
    for (edit in seq(1, n_edit)){
      message("Select mask vertex to edit, then it's new location")
      mask <- gglocator(2)
      mask_vertices_edit[nrow(mask_vertices_edit) + 1, ] <- mask[1, ]
      edit_distance <- as.matrix(dist(mask_vertices_edit))
      close_v <-
        which((edit_distance[nrow(mask_vertices_edit),] < threshold) &
                (edit_distance[nrow(mask_vertices_edit),] > 0), arr.ind = TRUE)

      while (length(close_v) == 0) {
        message("The point selected isn't close to any existing vertex, try again")
        mask <- gglocator(2)
        mask_vertices_edit[nrow(mask_vertices_edit) + 1, ] <- mask[1,]
        edit_distance <- as.matrix(dist(mask_vertices_edit))
        close_v <-
          which((edit_distance[nrow(mask_vertices_edit),] < threshold) &
                  (edit_distance[nrow(mask_vertices_edit),] > 0), arr.ind = TRUE)
      }

      for (n_mask in edit_list) {
        mask_vertices <- data.frame(st_coordinates(mask_org[[n_mask]])[, 1:2])
        new_mask_v <- mask_vertices
        mask_vertices[nrow(mask_vertices) + 1, ] <- mask[1, ]
        v_distance <- as.matrix(dist(mask_vertices))
        close_selected <- which((v_distance[nrow(mask_vertices), ] < threshold)
                                & (v_distance[nrow(mask_vertices), ] > 0),
                                arr.ind = TRUE)

        if (length(close_selected) > 0){

          if(close_selected[1] == 1){
            new_mask_v[c(1, nrow(new_mask_v)), ] <- mask[2, ]
            close_selected <- head(close_selected[-1], -1)
          }
          if(length(close_selected) > 0){
            new_mask_v[close_selected[1], ] <- mask[2, ]
            if(length(close_selected) > 1){
              new_mask_v <- new_mask_v[-tail(close_selected, -1), ]
            }
          }

          g <- g + geom_path(data = new_mask_v, aes(X, Y), color = color_v[edit + 1],
                             linewidth = 1)
          mask_org[[n_mask]] <-
            st_polygon(list(as.matrix(new_mask_v)))
        }
      }

      print(g)
    }

    # return pressure frames for selected area
    pressure_data[[5]] <- mask_org
    return(pressure_data)
  } else {
    message("You must edit at least one vertex, change the value for 'n_edit'")
  }
}


# =============================================================================

#' @title CPEI
#' @description Determine Center of Pressure Excursion Index (CPEI) for
#' footprint pressure data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement. Not currently available for pedar.
#' @param foot_side String. "right" or "left". Required for automatic detection
#' of points
#' @param plot_result Logical. Plots pressure image with COP and CPEI overlaid
#' @return Numeric. CPEI value
#' @examplesIf interactive()
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' cpei(pressure_data, foot_side = "auto", plot_result = FALSE)
#' @importFrom sf st_convex_hull st_linestring st_distance st_coordinates
#' st_intersection
#' @importFrom dplyr pull summarise
#' @export

cpei <- function(pressure_data, foot_side, plot_result = TRUE) {
  # check set up
  if (interactive() == FALSE)
    stop("we recommend user reviews each measurement")

  # global variables
  x <- y <- me <- X <- Y <- NULL

  # side
  if (foot_side == "auto") {
    side <- auto_detect_side(pressure_data)
  } else {
    side <- foot_side
  }

  # footprint coordinates
  sc_df <- pressure_data[[7]]

  # determine medial line
  ## medial edge coords
  unq_y <- unique(sc_df$y)
  med_edge <- data.frame(x = rep(NA, length.out = length(unq_y)),
                         y = unq_y)
  for (i in 1:length(unq_y)) {
    if (side == "LEFT") {
      med_edge[i, 1] <- sc_df %>% filter(y == unq_y[i]) %>%
        summarise(me = max(x)) %>% pull(me)
    } else {
      med_edge[i, 1] <- sc_df %>% filter(y == unq_y[i]) %>%
        summarise(me = min(x)) %>% pull(me)
    }
  }

  ## convex hull of medial edge
  med_edge_sf <- med_edge %>% st_as_sf(coords = c("x", "y"))
  med_edge_chull <- st_convex_hull(st_combine(med_edge_sf))

  ## longest med edge line
  me_dis <- as.matrix(dist(st_coordinates(med_edge_chull)))
  me_dis <- me_dis[row(me_dis) == (col(me_dis) - 1)]
  me_max <- order(me_dis)[length(me_dis)]
  med_side <- st_coordinates(med_edge_chull)[c(me_max, me_max + 1), c(1, 2)]
  med_side_line <- st_linestring(med_side)
  med_side_line <- st_extend_line(med_side_line, 0.1)

  ## confirm no points to medial side of line


  # bounding box for whole foot
  ## Bounding box
  mbb <- getMinBBox(as.matrix(sc_df[, c(1, 2)]))
  side1 <- mbb[c(1:2), ]
  side2 <- mbb[c(3:4), ]

  ## Get distal trisection points and make line
  trisection_left <- side2[1, ] + ((side2[2, ] - side2[1, ]) * 0.667)
  trisection_right <- side1[1, ] + ((side1[2, ] - side1[1, ]) * 0.667)
  trisection_df <- data.frame(x = c(trisection_left[1], trisection_right[1]),
                              y = c(trisection_left[2], trisection_right[2]))
  tri_line <- st_linestring(as.matrix(trisection_df))
  tri_line <- st_extend_line(tri_line, 1)
  #tri_int <- st_intersection(st_cast(fp_chull, "MULTILINESTRING"), tri_line)

  # trisection of medial line
  med_edge_tri_pt <- st_intersection(med_side_line, tri_line)

  # perp med line
  med_slope <- (med_side_line[2, 2] - med_side_line[1, 2]) /
    (med_side_line[2, 1] - med_side_line[1, 1])
  perp_med_line_slope <- -1 / med_slope
  tri_med_int <- med_edge_tri_pt[2] - (perp_med_line_slope * med_edge_tri_pt[1])
  perp_med_line <- rbind(med_edge_tri_pt, c(0, tri_med_int))
  perp_med_line <- st_extend_line(st_linestring(perp_med_line), 0.2)

  # ff width
  ## convex hull
  df.sf <- sc_df %>%
    st_as_sf(coords = c( "x", "y" ))
  fp_chull <- st_convex_hull(st_combine(df.sf))
  ff_lat_int <- st_intersection(fp_chull, perp_med_line)
  ff_lat_int <- st_coordinates(ff_lat_int)[, c(1, 2)]
  if (side == "RIGHT") {
    ff_lat_int <- ff_lat_int[which.max(ff_lat_int[, 1]), ]
  } else {
    ff_lat_int <- ff_lat_int[which.min(ff_lat_int[, 1]), ]
  }
  ff_width_pts <- rbind(ff_lat_int, st_coordinates(med_edge_tri_pt))
  ff_width <- as.numeric(unname(dist(ff_width_pts)))

  # cop_line
  cop_df <- cop(pressure_data)
  cop_sf <- cop_df %>%
    st_as_sf(coords = c("x_coord", "y_coord"))
  cop_chull <- st_convex_hull(st_combine(cop_sf))

  # cop_straight
  cop_dis <- as.matrix(dist(st_coordinates(cop_chull)))
  cop_dis <- cop_dis[row(cop_dis) == (col(cop_dis) - 1)]
  cop_max <- order(cop_dis)[length(cop_dis)]
  cop_side <- st_coordinates(cop_chull)[c(cop_max, cop_max + 1), c(1, 2)]
  cop_side_line <- st_linestring(cop_side)
  cop_side_line <- st_extend_line(cop_side_line, 0.1)

  # cross point cop and per_med_line
  cop1_per_med_pt <- st_intersection(st_linestring(as.matrix(cop_df)),
                                     perp_med_line)

  # cross point cop straight and per_med_line
  cop2_per_med_pt <- st_intersection(cop_side_line, perp_med_line)

  # CPE distance
  CPE <- as.numeric(unname(st_distance(cop1_per_med_pt, cop2_per_med_pt)))

  # calculate CPEI
  CPEI <- (CPE / ff_width) * 100

  # plot
  if (plot_result == TRUE) {
    # make CPEI dfs for plotting
    med_side_df <- as.data.frame(med_side)
    cop_side_df <- as.data.frame(cop_side)
    cpe_df <- as.data.frame(rbind(cop1_per_med_pt, cop2_per_med_pt))
    colnames(cpe_df) <- c("X", "Y")

    # plot
    g <- plot_pressure(pressure_data, plot_colors = "custom",
                       break_values = c(100, 750),
                       break_colors = c("lightgrey", "grey", "darkgrey"),
                       plot_COP = TRUE)
    #g <- g + geom_line(data = med_side_df, aes(X, Y), linetype = "dashed",
    #                   color = "black", size = 1.5)
    g <- g + geom_line(data = cpe_df, aes(X, Y), color = "black",
                       size = 2)
    g <- g + geom_line(data = cop_side_df, aes(X, Y), colour = "blue",
                       alpha = 0.8)
    print(g)
  }

  # manual select function
  #manually_select <- function(n_points, mess) {
  #  message(mess)
  #  x <- gglocator(n_points)
  #  colnames(x) <- c("x_coord", "y_coord")
  #  return(x)
  #}

  ## check if automatic identification worked
  # plot max footprint with additional points
  #g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)
  #g <- g + geom_point(data = heel_coord, aes (x = x_coord, y = y_coord), shape = 6)
  #g <- g + geom_point(data = toe_coord, aes (x = x_coord, y = y_coord), shape = 6)
  #g <- g + geom_point(data = start_point, aes (x = x_coord, y = y_coord), shape = 2)
  #g <- g + geom_point(data = end_point, aes (x = x_coord, y = y_coord), shape = 2)
  #g <- g + geom_line(data = m_bor, aes(x = x_coord, y = y_coord), colour = "purple")
  #g <- g + geom_line(data = l_bor, aes(x = x_coord, y = y_coord), colour = "orange")
  #print(g)
  #auto_worked <- readline("Have points been correctly identified? c: manually select cop; a: manually select all")

  ## if automatic identification failed, redo manually
  #if (auto_worked == "a") {
  #  # plot footprint
  #  g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)

  # select points
  #  m_bor <- manually_select(2, "select medial border: proximal point first, then distal")
  #  l_bor <- manually_select(2, "select lateral border: proximal point first, then distal")
  #  heel_coord <- manually_select(1, "select most proximal point of heel")
  #  toe_coord <- manually_select(1, "select most distal point of toes")
  #  cop_start <- manually_select(1, "select the most medial point near the start of the COP")
  #  cop_end <- manually_select(1, "select the most medial point near the end of the COP")
  #}

  #if (auto_worked == "c") {
  #  # plot footprint
  #  g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)

  # select points
  #  start_point <- manually_select(1, "select the most medial point near the start of the COP")
  #  end_point <- manually_select(1, "select the most medial point near the end of the COP")
  #}

  # return CPEI
  return(CPEI)
}


# =============================================================================

#' Analyze masked regions of pressure data
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param partial_sensors Logical Defines how sensors that do not
#'   lie wholly within mask are dealt with. If FALSE, they will be excluded;
#'   if TRUE, for relevant variables their contribution will be weighted by the
#'   proportion of the sensor that falls within the mask border
#' @param variable String. Variable to be determined. "press_peak_sensor",
#' "press_peak_mask", "contact_area_peak", "pti_1", "pti_2",
#' "force_time_integral", "force_peak", "dpli", "press_peak_sensor_ts", "force_ts".
#' Variables ending in "_ts" produce time series data
#' @param pressure_units String. Default "kPa". Other options: "MPa", "Ncm2"
#' (Newtons per square centimeter)
#' @param area_units String. Default is "cm2" (square centimeters). Other
#' options "m2" (square meters); "mm2" (square millimeters)
#' @return Data frame. Contains values for each mask plus additional information
#' relevant to the data including cycle/step and foot side
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- create_mask_auto(pressure_data, "automask_simple", plot = FALSE)
#' mask_analysis(pressure_data, FALSE, variable = "press_peak_sensor")
#' @importFrom sf st_intersects st_geometry st_area st_intersection
#' @importFrom pracma trapz
#' @export

mask_analysis <- function(pressure_data, partial_sensors = FALSE,
                          variable = "press_peak_sensor",
                          pressure_units = "kPa", area_units = "cm2") {
  # set global variables
  pedar_insole_type <- sens_poly <- act_sens <- area <- max_df <-
    mask_sides <- overlap_list <- time <- side <- cycle <- NULL

  # masks
  masks <- pressure_data[[5]]

  # events
  events <- pressure_data[[6]]

  # set up mask/sensor areas
  ## Make active sensors into polygons
  sens_poly <- sens_df_2_polygon(pressure_data[[7]])

  ## For each region mask, find which polygons intersect
  sens_mask_df <- matrix(rep(0, length.out = (length(sens_poly) * length(masks))),
                         nrow = length(sens_poly), ncol = length(masks))
  sens_mask_df <- as.data.frame(sens_mask_df)
  colnames(sens_mask_df) <- names(masks)
  for (mask in 1:length(masks)){
    for (j in 1:length(sens_poly)) {
      x <- st_intersects(masks[[mask]], sens_poly[[j]])
      if (identical(x[[1]], integer(0)) == FALSE) {
        y <- st_intersection(masks[[mask]], sens_poly[[j]])
        if (pressure_data[[2]] == "pedar") {
          sens_mask_df[j, mask] <- st_area(y) / st_area(sens_poly[[j]])
        } else {
          sens_mask_df[j, mask] <- st_area(y) / pressure_data[[3]][j]
        }
      }
    }
  }

  if (pressure_data[[2]] == "pedar") {
    # get mask sides
    mask_sides <- rep(NA, length.out = length(masks))
    for (mask in seq_along(masks)) {
      mask_coords <- st_coordinates(masks[[mask]])
      if (mean(mask_coords[, 1]) < 0.5) {
        mask_sides[mask] <- "LEFT"
      } else {
        mask_sides[mask] <- "RIGHT"
      }
    }
  }

  # include partial sensors?
  if (partial_sensors == FALSE) {
    sens_mask_df[sens_mask_df < 0.999 & sens_mask_df > 0] <- 0
  }

  # create blank output dataframes
  output_df <- data.frame(time = double(), cycle = integer(), side = factor(),
                          pressure_variable = factor(), mask_name = factor(),
                          value = double())

  # analysis
  ## output_names
  col_names <- colnames(output_df)

  # side
  time <- side <- NA

  # process
  if (length(pressure_data[[6]]) == 0) {
    n_cycle = 1
    has_events <- FALSE
  } else {
    n_cycle = nrow(pressure_data[[6]])
    has_events <- TRUE
  }
  for (cycle in 1:n_cycle) {
    # get step data
    pressure_data_ <- pressure_data[[1]]
    if (n_cycle >= 1 & has_events == TRUE) {
      cyc_str <- unname(unlist(pressure_data[[6]][cycle, 2]))
      cyc_end <- unname(unlist(pressure_data[[6]][cycle, 3]))
      pressure_data_ <- pressure_data_[c(cyc_str:cyc_end), ]
    }

    # analysis
    for (mask in seq_along(pressure_data[[5]])) {
      # mask name
      mn <- names(pressure_data[[5]])[mask]

      # Analyse regions for maximum value of any sensor within region during trial
      if (variable == "press_peak_sensor") {
        pv <- "peak_press_sensor"
        val <- pressure_peak(pressure_data_, sens_mask_df, mask)
      }

      if (variable == "press_peak_sensor_ts") {
        pv <- "peak_press_sensor_ts"
        val <- pressure_peak_ts(pressure_data_, sens_mask_df, mask)
      }

      # Analyse regions for maximum force during the trial
      if (variable == "force_peak") {
        pv <- "peak_force"
        val <- force_peak(pressure_data_, pressure_data, sens_mask_df, mask)
      }

      if (variable == "force_ts") {
        pv <- "force_ts"
        val <- force_ts(pressure_data_, pressure_data, sens_mask_df, mask)
      }

      # Analyse regions for peak regional pressure (defined as the pressure
      # averaged over the mask)
      if (variable == "press_peak_mask") {
        pv <- "press_peak_mask"
        val <- pressure_mean(pressure_data_, pressure_data, sens_mask_df, mask)
      }

      # Analyse regions for maximum contact area of region during trial
      if (variable == "contact_area_peak") {
        pv <- "contact_area"
        val <- contact_area(pressure_data_, pressure_data, sens_mask_df, mask)
      }

      # Analyse regions for pressure time integral (Novel definition)
      if (variable == "pti_1") {
        pv <- "pti_novel"
        val <- pti_1(pressure_data_, pressure_data, sens_mask_df, mask)
      }

      # Analyse regions for pressure time integral (Melai definition)
      if (variable == "pti_2") {
        pv <- "pti_melai"
        val <- pti_2(pressure_data_, pressure_data, sens_mask_df, mask)
      }

      # Analyse regions for force time integral
      if (variable == "fti") {
        pv <- "fti"
        val <- fti(pressure_data_, pressure_data, sens_mask_df, mask)
      }

      if (variable == "dpli") {
        pv <- "dpli"
        val <- dpli(pressure_data_, pressure_data, sens_mask_df, mask, 30)
      }

      # pedar adjustments to output df
      if (pressure_data[[2]] == "pedar") {
        if (mask_sides[mask] == events[cycle, 1]) {
          side <- mask_sides[mask]} else {side <- "None"}
      }

      # add to output df
      if (!str_ends(variable, "_ts")) {
        output_df <- rbind(output_df, c(time, cycle, side, pv, mn, val))
      } else {
        ns <- length(val)
        df <- data.frame(time = seq(0, by = pressure_data[[4]], length.out = ns),
                         cycle = rep(cycle, ns), side = rep(side, ns),
                         pressure_variable = rep(pv, ns),
                         mask_name = rep(mn, ns), value = val)
        output_df <- rbind(output_df, df)
      }
    }
  }

  # fix col names
  colnames(output_df) <- col_names
  if (pressure_data[[2]] == "pedar") {
    output_df <- output_df %>% filter(side != "None")
  }

  # Analyse regions for peak mask pressure (defined as the maximum mean pressure
  # of active sensors in region during the trial). Outputs 101 point vector (kPa)
  if (variable == "press_peak_mask_ts") {
    for (mask in seq_along(masks)) {
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens] * sensor_area * 1000
        CA <- rep(sensor_area, length.out = length(P))
        force <- sum(P * sens_mask_df[, mask])
        contact_area <- sum(CA * sens_mask_df[, mask])
        output_mat[i, mask] <- force / contact_area / 1000
      }
    }
    if (pressure_units == "MPa") {output_mat <- output_mat * 0.001}
    if (pressure_units == "Ncm2") {output_mat <- output_mat * 0.1}
  }

  # Analyse regions for contact area throughout the trial (outputs vector)
  if (variable == "contact_area_ts") {
    for (mask in seq_along(masks)) {
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens]
        P[P > 0] <- sensor_area
        cArea <- sum(P * sens_mask_df[, mask])
        if (area_units == "cm2") {cArea <- cArea * 10000}
        if (area_units == "mm2") {cArea <- cArea * 1000000}
        output_mat[i, mask] <- cArea
      }
    }
  }

  # adjust units
  #if (pressure_units == "MPa") {output_df$value <- output_df$value * 0.001}
  #if (pressure_units == "Ncm2") {output_df$value <- output_df$value * 0.1}

  # adjust df
  output_df$time <- as.numeric(output_df$time)
  output_df$cycle <- as.integer(output_df$cycle)
  output_df$side <- as.factor(output_df$side)
  output_df$pressure_variable <- as.factor(output_df$pressure_variable)
  output_df$mask_name <- as.factor(output_df$mask_name)
  output_df$value <- as.numeric(output_df$value)
  if (!str_ends(variable, "_ts")) {output_df <- subset(output_df, select = -c(time))}
  if (pressure_data[[2]] != "pedar") {output_df <- subset(output_df, select = -c(side))}
  if (length(events) == 0) {output_df <- subset(output_df, select = -c(cycle))}

  # return
  return(output_df)
}


# =============================================================================

#' Calculate Arch Index.
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param plot Logical. Not implemented yet
#' @return Numeric. Arch index value
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' arch_index(pressure_data)
#' @importFrom sf st_bbox st_intersection
#' @export

arch_index <- function(pressure_data, plot = TRUE) {
  # set global variables
  id <- NULL

  # check emed
  if (!(pressure_data[[2]] == "emed" | pressure_data[[2]] == "footscan" |
        pressure_data[[2]] == "pliance"))
  stop("This function currently only works for pressure plate data")

  # side
  side <- auto_detect_side(pressure_data)

  # cop_line
  cop_df <- cop(pressure_data)
  cop_sf <- cop_df %>%
    st_as_sf(coords = c("x_coord", "y_coord"))
  cop_chull <- st_convex_hull(st_combine(cop_sf))

  # cop_straight
  cop_dis <- as.matrix(dist(st_coordinates(cop_chull)))
  cop_dis <- cop_dis[row(cop_dis) == (col(cop_dis) - 1)]
  cop_max <- order(cop_dis)[length(cop_dis)]
  cop_side <- st_coordinates(cop_chull)[c(cop_max, cop_max + 1), c(1, 2)]

  # remove toes
  sens_poly <- sens_df_2_polygon(pressure_data[[7]])
  mask_fp <- create_mask_auto(pressure_data, "automask_novel",
                              res_value = 100000, plot = FALSE)
  mask_fp <- mask_fp[[5]]
  not_in_mask <- c()
  for (j in 1:length(sens_poly)) {
    x <- st_intersects(mask_fp[[4]], sens_poly[[j]])
    if (identical(x[[1]], integer(0)) == FALSE) {not_in_mask <- c(not_in_mask, j)}
  }
  for (j in 1:length(sens_poly)) {
    x <- st_intersects(mask_fp[[5]], sens_poly[[j]])
    if (identical(x[[1]], integer(0)) == FALSE) {not_in_mask <- c(not_in_mask, j)}
  }
  not_in_mask <- unique(not_in_mask)

  # rotate sensor coords
  ## angle
  ang <-  atan((cop_side[1, 1] - cop_side[2, 1]) /
                 (cop_side[1, 2] - cop_side[2, 2])) * (180/pi)

  ## rotate points
  sens_df <- pressure_data[[7]]
  sens_df <- sens_df %>% filter(!id %in% not_in_mask)
  pts <- sens_df[, c(1, 2)]
  rotated_pts <- rot_pts(as.matrix(pts), ang)
  sens_df[, c(1, 2)] <- rotated_pts

  # bounding box of footprint (no toes)
  bb_poly <- as.data.frame(rotated_pts) %>% st_as_sf(coords = c("x", "y")) %>%
    st_bbox() %>% st_as_sfc()

  # split into thirds
  coords <- st_coordinates(bb_poly)
  bb_third <- (max(coords[, 2]) - min(coords[, 2])) / 3
  len_33 <- min(coords[, 2]) + bb_third
  len_66 <- min(coords[, 2]) + bb_third * 2
  x_min <- min(coords[, 1])
  x_max <- max(coords[, 1])
  line_33 <- st_linestring(matrix(c(x_min, x_max, len_33, len_33), ncol = 2))
  line_33 <- st_extend_line(line_33, 1)
  poly_33_up <- st_line2polygon(line_33, 1, "+Y")
  poly_33_down <- st_line2polygon(line_33, 1, "-Y")
  line_66 <- st_linestring(matrix(c(x_min, x_max, len_66, len_66), ncol = 2))
  line_66 <- st_extend_line(line_66, 1)
  poly_66_up <- st_line2polygon(line_66, 1, "+Y")
  poly_66_down <- st_line2polygon(line_66, 1, "-Y")

  ## ff
  ff_poly <- st_difference(bb_poly, poly_66_down)

  ## mf
  mf_poly <- st_difference(bb_poly, poly_66_up)
  mf_poly <- st_difference(mf_poly, poly_33_down)

  ## hf
  hf_poly <- st_difference(bb_poly, poly_33_up)

  # calculate areas
  sens_poly <- sens_df_2_polygon(sens_df)
  sens_area <- pressure_data[[3]][1]
  ff_area <- c()
  for (j in 1:length(sens_poly)) {
    x <- st_intersects(ff_poly, sens_poly[[j]])
      if (identical(x[[1]], integer(0)) == FALSE) {
        y <- st_intersection(ff_poly, sens_poly[[j]])
        ff_area <- c(ff_area, (st_area(y) / sens_area) * sens_area)
    }
  }
  ff_area <- sum(ff_area)
  mf_area <- c()
  for (j in 1:length(sens_poly)) {
    x <- st_intersects(mf_poly, sens_poly[[j]])
    if (identical(x[[1]], integer(0)) == FALSE) {
      y <- st_intersection(mf_poly, sens_poly[[j]])
      mf_area <- c(mf_area, (st_area(y) / sens_area) * sens_area)
    }
  }
  mf_area <- sum(mf_area)
  hf_area <- c()
  for (j in 1:length(sens_poly)) {
    x <- st_intersects(hf_poly, sens_poly[[j]])
    if (identical(x[[1]], integer(0)) == FALSE) {
      y <- st_intersection(hf_poly, sens_poly[[j]])
      hf_area <- c(hf_area, (st_area(y) / sens_area) * sens_area)
    }
  }
  hf_area <- sum(hf_area)

  # calculate arch index
  a_ind <- mf_area / (ff_area + mf_area + hf_area)

  # return arch index
  return(a_ind)
}


# =============================================================================
# =============================================================================

# helper functions

#' @title Get outline of pressure region
#' Determine outline (convex hull) of pressure measurement
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param pressure_image String. "max" = footprint of maximum sensors. "frame"
#' = an individual frame
#' @param frame Integer. Frame number to use
#' @return Polygon. sfg object representing convex hull outline
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_outline(pressure_data)
#' pressure_outline(pressure_data, frame = 50)
#' @importFrom magrittr "%>%"
#' @importFrom sf st_as_sf st_convex_hull st_combine st_buffer
#' @noRd

pressure_outline <- function(pressure_data, pressure_image = "max", frame) {
  # check inputs
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure data must contain array")

  # check that this is not pedar data
  if (pressure_data[[2]] == "pedar")
    stop("data cannot be from pedar")

  # determine active sensor coordinates
  if (pressure_image == "max") {
    coords <- sensor_coords(pressure_data, pressure_image = "max")
  } else {
    coords <- sensor_coords(pressure_data, "frame", frame)
  }

  # calculate convex hull
  sens_coords_df <- coords %>%
    st_as_sf(coords = c("x_coord", "y_coord"))
  fp_chull <- st_convex_hull(st_combine(sens_coords_df))

  # add buffer here
  buffer_size <- pressure_data[[3]][1] / 2
  fp_chull_buffer <- st_buffer(fp_chull, buffer_size)

  # return convex hull coordinates
  return(fp_chull_buffer)
}


#' interpolation function
#' @importFrom stats approx
#' @noRd

approxP <- function(x, interp_to) {
  y <- approx(x, n = interp_to)
  y$x <- NULL
  y <- unname(unlist(y))
  return(y)
}


#' @title Get minimum bounding box (angle)
#' @description Create a bounding box around sensor coordinate array. Adapted
#' from function in shotGroups package
#' @param xy Matrix. n row by 2 column with x-y coords of points to be
#' bounded
#' @return List
#' @importFrom grDevices chull
#' @noRd

getMinBBox <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)

  ## rotating calipers algorithm using the convex hull
  H    <- chull(xy)           # hull indices, vertices ordered clockwise
  n    <- length(H)           # number of hull vertices
  hull <- xy[H, ]             # hull vertices

  ## unit basis vectors for all subspaces spanned by the hull edges
  hDir  <- diff(rbind(hull, hull[1,])) # account for circular hull vertices
  hLens <- sqrt(rowSums(hDir^2))       # length of basis vectors
  huDir <- diag(1/hLens) %*% hDir      # scaled to unit length

  ## unit basis vectors for the orthogonal subspaces
  ## rotation by 90 deg -> y' = x, x' = -y
  ouDir <- cbind(-huDir[ , 2], huDir[ , 1])

  ## project hull vertices on the subspaces spanned by the hull edges, and on
  ## the subspaces spanned by their orthogonal complements - in subspace coords
  projMat <- rbind(huDir, ouDir) %*% t(hull)

  ## range of projections and corresponding width/height of bounding rectangle
  rangeH  <- matrix(numeric(n * 2), ncol = 2)   # hull edge
  rangeO  <- matrix(numeric(n * 2), ncol = 2)   # orth subspace
  widths  <- numeric(n)
  heights <- numeric(n)
  for(i in seq(along = H)) {
    rangeH[i, ] <- range(projMat[  i, ])
    rangeO[i, ] <- range(projMat[n + i, ])  # orth subspace is in 2nd half
    widths[i]   <- abs(diff(rangeH[i, ]))
    heights[i]  <- abs(diff(rangeO[i, ]))
  }

  ## extreme projections for min-area rect in subspace coordinates
  eMin  <- which.min(widths * heights)   # hull edge leading to minimum-area
  hProj <- rbind(   rangeH[eMin, ], 0)
  oProj <- rbind(0, rangeO[eMin, ])

  ## move projections to rectangle corners
  hPts <- sweep(hProj, 1, oProj[ , 1], "+")
  oPts <- sweep(hProj, 1, oProj[ , 2], "+")

  ## corners in standard coordinates, rows = x,y, columns = corners
  ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
  basis <- cbind(huDir[eMin, ], ouDir[eMin, ])  # basis formed by hull edge and orth
  hCorn <- basis %*% hPts
  oCorn <- basis %*% oPts
  pts   <- t(cbind(hCorn, oCorn[ , c(2, 1)]))

  # Which sides of BBox are longest?
  if (dist(pts[c(1, 2), ]) > dist(pts[c(2, 3), ])) {
    side1 <- pts[c(1, 2), ]
    side2 <- pts[c(3, 4), ]
  } else {
    side1 <- pts[c(2, 3), ]
    side2 <- pts[c(4, 1), ]
  }

  # order so lowest first
  side1 <- side1[order(side1[, 2]), ]
  side2 <- side2[order(side2[, 2]), ]

  # which side is L/R
  if (side1[1, 1] > side2[1, 1]) {
    side1_ <- side1 # side1_ is right
    side2_ <- side2
  } else {
    side1_ <- side2
    side2_ <- side1
  }

  mat <- rbind(side1_, side2_)

  # return
  return(mat)
}


#' @title masks (sf format) to dataframe
#' @description converts mask coordinates into a data frame
#' @param masks List
#' @returns Data frame
#' @importFrom dplyr bind_rows
#' @noRd

masks_2_df <- function(masks) {
  # no. of masks
  n_masks <- length(masks)

  # empty df
  df <- data.frame(mask = factor(),
                   x = double(),
                   y = double())

  # get coordinates of each mask and store in df
  for (mask in seq_along(masks)) {
    coords <- st_coordinates(masks[[mask]])[, c(1:2)]
    mask_name <- rep(names(masks)[mask], length.out = nrow(coords))
    df_ <- data.frame(mask = mask_name, x = coords[, 1], y = coords[, 2])
    df <- bind_rows(df, df_)
  }

  # return
  return(df)
}


#' @title Get coordinates of active sensors
#' @description Produces a data frame with coordinates of sensors
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param pressure_image. String. Which pressure image to use. Options are
#' "all_active", "all", or "frame".
#' @param frame Numeric. If pressure image is frame, the numeric value should be
#' provided here
#' @return Data frame. x and y coordinates of sensors
#' @examples
#' coords <- sensor_coords(pressure_data)
#' @noRd

sensor_coords <- function(pressure_data, pressure_image = "all_active", frame) {
  # pressure image
  if (pressure_image == "all" | pressure_image == "all_active" |
      pressure_image == "max" | pressure_image == "mean") {
    sens <- footprint(pressure_data, variable = "max")
  }
  if (pressure_image == "frame") {
    sens <- footprint(pressure_data, variable = "frame", frame = frame)
  }

  # dimensions
  sens_x <- abs(unique(pressure_data$sens_polygons$x)[2] - unique(pressure_data$sens_polygons$x)[1])
  sens_y <- abs(unique(pressure_data$sens_polygons$y)[2] - unique(pressure_data$sens_polygons$y)[1])

  # data frame with active sensors as coordinates
  x_cor <- seq(from = sens_x / 2, by = sens_x, length.out = ncol(sens))
  x_cor <- rep(x_cor, each = nrow(sens))
  y_cor <- seq(from = (sens_y / 2) + ((nrow(sens) - 1) * sens_y),
               by = (-1 * sens_y), length.out = nrow(sens))
  y_cor <- rep(y_cor, times = ncol(sens))
  coords <- data.frame(x_coord = x_cor, y_coord = y_cor)

  # remove inactive sensors
  P <- as.vector(sens)
  if (pressure_image != "all") {
    coords <- coords[which(P > 0), ]
  }

  # return sensor coordinates
  return(coords)
}


#' @title sensor array to polygon
#' @description converts sensor positions (coordinates) to polygons
#' @param pressure_data List
#' @param pressure_image String
#' @param frame Numeric
#' @param output String. "df" dataframe for plotting or "sf" shape poly for analysis
#' @importFrom sf st_polygon st_coordinates
#' @noRd

sensor_2_polygon <- function(pressure_data, pressure_image = "all_active",
                             frame = NA, output = "sf") {
  # set global variables
  pedar_insole_grid <- NULL

  # sensor coordinates
  sens_polygons <- list()
  if (pressure_data[[2]] != "pedar") {
    sens_coords <- sensor_coords(pressure_data, pressure_image, frame)

    # sensor dimensions
    width <- pressure_data[[3]][1] / 2
    height <- pressure_data[[3]][2] / 2

    # get corner points and make into polygon
    for (sens in 1:nrow(sens_coords)) {
      x1 <- sens_coords[sens, 1] - width
      y1 <- sens_coords[sens, 2] + height
      x2 <- sens_coords[sens, 1] + width
      y2 <- sens_coords[sens, 2] + height
      x3 <- sens_coords[sens, 1] + width
      y3 <- sens_coords[sens, 2] - height
      x4 <- sens_coords[sens, 1] - width
      y4 <- sens_coords[sens, 2] - height
      if (y1 < 0) {y1 <- 0}
      if (y2 < 0) {y2 <- 0}
      if (y3 < 0) {y3 <- 0}
      if (y4 < 0) {y4 <- 0}
      sens_polygons[[sens]] <- st_polygon(list(matrix(c(x1, x2, x3, x4, x1,
                                                        y1, y2, y3, y4, y1),
                                                      5, 2)))
    }
  } else {
    pedar_insole_grid <- pedar_insole_grids()
    rs <- c(1:99, 101:199)
    for (sens in 1:length(rs)) {
      x1 <- pedar_insole_grid[rs[sens], 1]
      y1 <- pedar_insole_grid[rs[sens], 2]
      x2 <- pedar_insole_grid[rs[sens], 3]
      y2 <- pedar_insole_grid[rs[sens], 4]
      x3 <- pedar_insole_grid[rs[sens], 5]
      y3 <- pedar_insole_grid[rs[sens], 6]
      x4 <- pedar_insole_grid[rs[sens], 7]
      y4 <- pedar_insole_grid[rs[sens], 8]
      sens_polygons[[sens]] <- st_polygon(list(matrix(c(x1, x2, x3, x4, x1,
                                                        y1, y2, y3, y4, y1),
                                                      5, 2)))
    }
  }

  # make df if required
  if (output == "df") {
    id_df <- data.frame(x = double(),
                        y = double(),
                        z = integer())
    # make into identity df
    for (i in 1:length(sens_polygons)) {
      mat <- st_coordinates(sens_polygons[[i]])[, c(1, 2)]
      mat_df <- data.frame(mat)
      id <- rep(i, length.out = nrow(mat_df))
      mat_df <- cbind(mat_df, id)
      colnames(mat_df) <- c("x", "y", "id")
      id_df <- rbind(id_df, mat_df)
    }
  }

  # check for minus values

  # return
  if (output == "sf") {return(sens_polygons)}
  if (output == "df") {return(id_df)}
}

sensor_2_polygon2 <- function(sensor_array) {
  # empty df
  id_df <- data.frame(x = double(),
                      y = double(),
                      z = integer())

  for (i in 1:nrow(sensor_array)) {
    mat <- matrix(c(sensor_array[i, ], sensor_array[i, c(1, 2)]), 5, 2, byrow = TRUE)
    mat_df <- data.frame(mat)
    id <- rep(i, length.out = nrow(mat_df))
    mat_df <- cbind(mat_df, id)
    colnames(mat_df) <- c("x", "y", "id")
    id_df <- rbind(id_df, mat_df)
  }

  # return
  return(id_df)
}

sensor_2_polygon3 <- function(max_mat, width, height, output = "df") {
  # empty list to store polygons
  sens_polygons <- list()

  # sensor coordinates
  ## centroids for all sensors
  x_cor <- seq(from = width / 2, by = width, length.out = ncol(max_mat))
  x_cor <- rep(x_cor, each = nrow(max_mat))
  y_cor <- seq(from = (height / 2) + ((nrow(max_mat) - 1) * height),
               by = (-1 * height), length.out = nrow(max_mat))
  y_cor <- rep(y_cor, times = ncol(max_mat))
  sens_coords <- data.frame(x_coord = x_cor, y_coord = y_cor)

  ## remove inactive
  P <- as.vector(max_mat)
  sens_coords <- sens_coords[which(P > 0), ]

  # make polygons
  wid <- width / 2
  ht <- height / 2
  # get corner points and make into polygon
  for (sens in 1:nrow(sens_coords)) {
    x1 <- sens_coords[sens, 1] - wid
    y1 <- sens_coords[sens, 2] + ht
    x2 <- sens_coords[sens, 1] + wid
    y2 <- sens_coords[sens, 2] + ht
    x3 <- sens_coords[sens, 1] + wid
    y3 <- sens_coords[sens, 2] - ht
    x4 <- sens_coords[sens, 1] - wid
    y4 <- sens_coords[sens, 2] - ht
    if (y1 < 0) {y1 <- 0}
    if (y2 < 0) {y2 <- 0}
    if (y3 < 0) {y3 <- 0}
    if (y4 < 0) {y4 <- 0}
    sens_polygons[[sens]] <- st_polygon(list(matrix(c(x1, x2, x3, x4, x1,
                                                      y1, y2, y3, y4, y1),
                                                    5, 2)))
  }

  # make into identity df
  ## empty df
  id_df <- data.frame(x = double(),
                      y = double(),
                      z = integer())

  for (i in 1:length(sens_polygons)) {
    mat <- st_coordinates(sens_polygons[[i]])[, c(1, 2)]
    mat_df <- data.frame(mat)
    id <- rep(i, length.out = nrow(mat_df))
    mat_df <- cbind(mat_df, id)
    colnames(mat_df) <- c("x", "y", "id")
    id_df <- rbind(id_df, mat_df)
  }

  # if sf required
  if (output == "sf") {id_df <- sens_polygons}

  # return sensor coordinates
  return(id_df)
}

sensor_2_polygon4 <- function() {
  sens_polygons <- list()
  pedar_insole_grid <- pedar_insole_grids()
  rs <- c(1:198)
  for (sens in 1:length(rs)) {
    x1 <- pedar_insole_grid[rs[sens], 1]
    y1 <- pedar_insole_grid[rs[sens], 2]
    x2 <- pedar_insole_grid[rs[sens], 3]
    y2 <- pedar_insole_grid[rs[sens], 4]
    x3 <- pedar_insole_grid[rs[sens], 5]
    y3 <- pedar_insole_grid[rs[sens], 6]
    x4 <- pedar_insole_grid[rs[sens], 7]
    y4 <- pedar_insole_grid[rs[sens], 8]
    sens_polygons[[sens]] <- st_polygon(list(matrix(c(x1, x2, x3, x4, x1,
                                                      y1, y2, y3, y4, y1),
                                                    5, 2)))
  }

  # reorder
  sens_polygons <- sens_polygons[c(100:198, 1:99)]

  # make into identity df
  ## empty df
  id_df <- data.frame(x = double(),
                      y = double(),
                      z = integer())

  for (i in 1:length(sens_polygons)) {
    mat <- st_coordinates(sens_polygons[[i]])[, c(1, 2)]
    mat_df <- data.frame(mat)
    id <- rep(i, length.out = nrow(mat_df))
    mat_df <- cbind(mat_df, id)
    colnames(mat_df) <- c("x", "y", "id")
    id_df <- rbind(id_df, mat_df)
  }

  # return sensor coordinates
  return(id_df)
}

sens_df_2_polygon <- function(sens_polygons) {
  id <- NULL
  sensor_polys <- list()
  sens_n <- unique(sens_polygons[, 3])
  for (sens in seq_along(sens_n)) {
    s_df <- sens_polygons %>% filter(id == sens_n[sens])
    sensor_polys[[sens]] <- st_polygon(list(matrix(c(s_df[1, 1], s_df[2, 1],
                                                     s_df[3, 1], s_df[4, 1],
                                                     s_df[1, 1], s_df[1, 2],
                                                     s_df[2, 2], s_df[3, 2],
                                                     s_df[4, 2], s_df[1, 2]),
                                                   5, 2)))
  }
  return(sensor_polys)
}

#' @title Sensor area
#' @description Calculates area of sensors from corner vertices
#' @param sensor_df Data frame. Contains sensor corner vertices
#' @return Vector
#' @noRd

sensor_area <- function(sensor_df) {
  id <- NULL

  # number of sensors
  n_sensors <- unique(sensor_df[, 3])

  # get areas
  areas <- rep(NA, length(n_sensors))
  for (sens in seq_along(n_sensors)) {
    sens_mat <- as.matrix(sensor_df %>% filter(id == n_sensors[sens]))[, c(1, 2)]
    areas[sens] <- st_area(st_polygon(list(sens_mat)))
  }

  # return
  return(areas)
}


#' @title pedar force
#' @description generates force curve from pedar data
#' @param pressure_data List.
#' @param variable "force", "area"
#' @param min_press Numeric. Threshold value for sensor to be considered active.
#' Currently only applies to insole data
#' @param threshold Numeric. For area calculations, minimum threshold for sensor to
#' be active
#' @return Vector
#' @noRd

force_pedar <- function(pressure_data, variable = "force", min_press = 10) {
  # check this is pedar data
  if (pressure_data[[2]] != "pedar")
    stop("must be pedar data")

  # split data
  pressure_data_R <- pressure_data[[1]][, 100:198]
  pressure_data_L <- pressure_data[[1]][, 1:99]

  # calculate force
  if (variable == "force") {
    force_right <- rowSums(pressure_data_R %*% diag(pressure_data[[3]][100:198])) * 1000
    force_left <- rowSums(pressure_data_L %*% diag(pressure_data[[3]][1:99])) * 1000
    output_df <- cbind(force_right, force_left)
  }

  # calculate area
  if (variable == "area") {
    area_right <- pressure_data_R > min_press
    area_right <- rowSums(area_right * pressure_data[[3]][100:198])
    area_left <- pressure_data_L > min_press
    area_left <- rowSums(area_left * pressure_data[[3]][1:99])
    output_df <- cbind(area_right, area_left)
  }

  # return
  return(output_df)
}


#' @title plot pedar
#' @description produces dataframe that plot_pressure uses to plot pedar data
#' @param pressure_data List.
#' @param pressure_image String. "max" = maximum values for full trial;
#' "mean" = mean value across trial; "step_max"; "step_mean"; "meanmax"
#' @param foot_side String. "both", "left", or "right"
#' @param step_n Numeric. If pressure_image is "step_max", the step number to
#' analyze
#' @return df
#' @noRd

plot_pedar <- function(pressure_data, pressure_image = "max",
                       foot_side = "left", step_n) {
  # set global variables
  pedar_insole_grid <- x <- y <- id <- NULL

  # check this is pedar (or other suitable) data
  if (!(pressure_data[[2]] == "pedar"))
    stop("data should be from pedar")

  # get L+ R data frames
  pressure_R_mat <- pressure_data[[1]][, 100:198]
  pressure_L_mat <- pressure_data[[1]][, 1:99]

  # pressure image
  if (pressure_image == "max") {
    R_data <- apply(pressure_R_mat, 2, max)
    L_data <- apply(pressure_L_mat, 2, max)
    comb_data <- c(L_data, R_data)
  }

  # separate into steps
  if (pressure_image == "step_max") {
    events <- pressure_data[[6]]
    pressure_R_mat <- pressure_R_mat[c(events[step_n, 2]:events[step_n, 3]), ]
    pressure_L_mat <- pressure_L_mat[c(events[step_n, 2]:events[step_n, 3]), ]
    R_data <- apply(pressure_R_mat, 2, max)
    L_data <- apply(pressure_L_mat, 2, max)
    comb_data <- c(L_data, R_data)
  }

  # average across steps
  if (pressure_image == "meanmax") {
    comb_data <- footprint(pressure_data, "meanmax")
  }

  # make ids
  ids <- 1:198

  # make df
  df <- data.frame(id = ids, value = comb_data)

  # add coordinates to df
  position <- pressure_data[[7]]
  df <- merge(df, position, by = c("id"))

  # add color
  df <- generate_colors(df)

  # side
  if (foot_side == "left") {
    df <- df[397:792,]
  }
  if (foot_side == "right") {
    df <- df[1:396,]
  }

  # return
  return(df)
}


#' @title Generate colors
#' @description Let's the user prescribe the color scale for plots
#' @param df Dataframe.
#' @param col_type String. "default": novel color scheme; "custom": user
#' supplied
#' @param break_values Vector. Vector with break values to be used. Should be
#' one longer than break_colors
#' @param break_colors Vector. Vector with colors to be used. Should be one
#' shorter than break_values
#' @return Data frame.
#' @noRd

generate_colors <- function(df, col_type = "default", break_values,
                            break_colors) {
  if (col_type == "default") {
    break_values <- c(0, 40, 60, 100, 150, 220, 300, 1000000)
    break_colors <- c("grey","lightblue", "darkblue","green","yellow",
                      "red","pink")
  } else {
    # check break_values and break_colors
    if (break_values[1] <= 0)
      stop("first break value should be >0")
    if(length(break_values) != (length(break_colors) - 1))
      stop("break_values should be one shorter than break_colors")
    if ("white" %in% break_colors)
      stop("break_colors cannot contain white")

    # add high final value to break_values
    break_values <- c(-1, 0, break_values, 1000000)
    break_colors <- c("white", break_colors)
  }

  # add col column
  df$cols <- cut(df$value, breaks = break_values, labels = break_colors)

  # return
  return(df)
}

#' @title extend st line
#' @param line sf linestring
#' @param distance Numeric.
#' @param end String
#' @importFrom sf st_sfc st_crs st_coordinates
#' @noRd

st_extend_line <- function(line, distance, end = "BOTH") {
  if (!(end %in% c("BOTH", "HEAD", "TAIL")) | length(end) != 1)
    stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")

  M <- st_coordinates(line)[, 1:2]
  keep <- !(end == c("TAIL", "HEAD"))

  ends <- c(1, nrow(M))[keep]
  i <- c(2, nrow(M) - 1)
  j <- c(1, -1)
  headings <- mapply(i, j, FUN = function(i, j) {
    Ax <- M[i - j, 1]
    Ay <- M[i - j, 2]
    Bx <- M[i, 1]
    By <- M[i, 2]
    unname(atan2(Ay - By, Ax - Bx))
  })
  #headings <- st_ends_heading(line)[keep]
  distances <- if (length(distance) == 1) rep(distance, 2) else rev(distance[1:2])

  M[ends,] <- M[ends,] + distances[keep] * c(cos(headings), sin(headings))
  newline <- st_linestring(M)

  # If input is sfc_LINESTRING and not sfg_LINESTRING
  if (is.list(line)) newline <- st_sfc(newline, crs = st_crs(line))

  return(newline)
}


#' @title line to polygon
#' @description Extrude line to form polygon
#' @param mat Matrix. xy points describing line or polyline
#' @param distance Numeric. Distance to extrude line
#' @param direction String. "+X" "-X", "+y", or "-Y"
#' @return Polygon. sf polygon object
#' @importFrom sf st_polygon
#' @noRd

st_line2polygon <- function(mat, distance, direction) {
  # get ends
  if (mat[1, 1] >= mat[nrow(mat), 1]) {dir1 <- "RL"} else {dir1 <- "LR"}
  if (mat[1, 2] >= mat[nrow(mat), 2]) {dir2 <- "TB"} else {dir2 <- "BT"}

  # +Y
  if (direction == "+Y" ) {
    mat_pts <- rbind(mat, c(mat[nrow(mat), 1], mat[nrow(mat), 2] + distance),
                     c(mat[1, 1], mat[1, 2] + distance), mat[1, ])
    if (dir1 == "RL") {mat <- mat[nrow(mat):1, ]}
  }
  # -Y
  if (direction == "-Y") {
    mat_pts <- rbind(mat, c(mat[nrow(mat), 1], mat[nrow(mat), 2] - distance),
                     c(mat[1, 1], mat[1, 2] - distance), mat[1, ])
    if (dir1 == "LR") {mat <- mat[nrow(mat):1, ]}
  }
  # +X
  if (direction == "+X") {
    mat_pts <- rbind(mat, c(mat[nrow(mat), 1] + distance, mat[nrow(mat), 2]),
                     c(mat[1, 1] + distance, mat[1, 2]), mat[1, ])
    if (dir2 == "TB") {mat <- mat[nrow(mat):1, ]}
  }

  # -X
  if (direction == "-X") {
    mat_pts <- rbind(mat, c(mat[nrow(mat), 1] - distance, mat[nrow(mat), 2]),
                     c(mat[1, 1] - distance, mat[1, 2]), mat[1, ])
    if (dir2 == "BT") {mat <- mat[nrow(mat):1, ]}
  }

  # make polygon
  mat_poly <- st_polygon(list(mat_pts))
}


#' @title low cost path
#' @description Calculates low cost path
#' @param mat Matrix.
#' @param offset_distance Numeric. If matrix does not start at 0, 0
#' @param start Vector. x and y coords
#' @param end Vector. x and y coords
#' @param base_x Numeric
#' @param base_y Numeric
#' @param res_scale Numeric
#' @param extend_line Logical. extends line horizontally
#' @importFrom sf st_cast st_join st_sf
#' @importFrom raster raster adjacent extent ncell
#' @importFrom gdistance transition geoCorrection shortestPath
#' @importFrom dplyr distinct group_by mutate n summarize across everything
#' @importFrom stats var
#' @noRd

shortest_path <- function(mat, offset_distance, start, end, base_x, base_y,
                          res_scale, extend_line = TRUE) {
  ## scale
  mat_norm <- mat / res_scale
  r <- raster(mat_norm)
  raster::extent(r) <- c(base_x * -1, base_x * (ncol(mat_norm) + 1),
                         offset_distance,
                         offset_distance + base_y * nrow(mat_norm))

  # line
  heightDiff <- function(x){x[2] - x[1]}
  hd <- suppressWarnings(transition(r, heightDiff, 8, symm = FALSE))
  slope <- geoCorrection(hd, type = "c")
  adj <- adjacent(r, cells = 1:ncell(r), pairs = TRUE, directions = 8)
  speed <- slope
  speed[adj] <- 6 * exp(-3.5 * abs(slope[adj]))
  Conductance <- geoCorrection(speed)
  AtoB <- shortestPath(Conductance, start, end, output = "SpatialLines")
  line_mat <- st_coordinates(st_as_sf(AtoB))[, c(1, 2)]

  # extend line
  if (extend_line == TRUE) {
    if (line_mat[1, 1] > line_mat[nrow(line_mat), 1]) {
      line_mat <- rbind(c(line_mat[1, 1] + 1, line_mat[1, 2]), line_mat,
                        c(line_mat[nrow(line_mat), 1] - 1,
                          line_mat[nrow(line_mat), 2]))
    } else {
      line_mat <- rbind(c(line_mat[1, 1] - 1, line_mat[1, 2]), line_mat,
                        c(line_mat[nrow(line_mat), 1] + 1,
                          line_mat[nrow(line_mat), 2]))
    }
  }

  # remove colinear points
  is_colinear <- function(p1, p2, p3) {
    area <- abs((p1[1] * p2[2] + p2[1] * p3[2] + p3[1] * p1[2]) -
      (p1[2] * p2[1] + p2[2] * p3[1] + p3[2] * p1[1])) / 2
    return(area < 0.0000001)
  }
  line_df <- data.frame(x = line_mat[, 1], y = line_mat[, 2])
  colinear_inds <- sapply(2:(nrow(line_df) - 1),
                          function(i) is_colinear(line_df[i - 1, ], line_df[i, ],
                                                  line_df[i + 1,]))
  line_mat <- line_mat[!(c(FALSE, colinear_inds)), ]

  # return
  return(line_mat)
}


#' @title Get toe line
#' @description Calculates line proximal to toes
#' @param pressure_data List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time
#' @importFrom sf st_cast st_join st_sf
#' @importFrom raster raster adjacent extent ncell
#' @importFrom gdistance transition geoCorrection shortestPath
#' @importFrom dplyr distinct group_by mutate n summarize across everything
#' @importFrom stats var
#' @noRd

toe_line <- function(pressure_data, side, res_scale = 10000) {
  # global variables
  id <- NULL

  # get max footprint and take top half
  pf_max <- pressure_data[[8]]

  # sensor dimensions
  sens_1 <- pressure_data[[7]] %>%
    filter(id == unique(pressure_data[[7]]$id)[50])
  base_x <- max(sens_1$x) - min(sens_1$x)
  base_y <- max(sens_1$y) - min(sens_1$y)

  # take top half
  pf_max_top <- pf_max[1:(floor(nrow(pf_max) / 2)), ]
  offset_row <- nrow(pf_max) - nrow(pf_max_top)
  offset_distance <- offset_row * base_y

  # remove islands
  ## get vertex coords
  verts <- pressure_data[[7]]

  ## remove duplicates (id included)
  verts <- verts %>% distinct()

  # toe line
  ## start and end points
  active_cols <- which(apply(pf_max_top, 2, var) != 0)
  active_rows <- which(rowSums(pf_max_top) > 0)
  active_row_1 <- active_rows[1]
  row_25 <- round((nrow(pf_max_top) - active_row_1) * 0.25)
  row_70 <- round((nrow(pf_max_top) - active_row_1) * 0.6)

  ## start point
  start_y <- offset_distance + (nrow(pf_max_top) * base_y) -
    (row_25 * base_y) - (active_row_1 * base_y)
  end_y <- offset_distance + (nrow(pf_max_top) * base_y) -
    (row_70 * base_y) - (active_row_1 * base_y)
  if (side == "LEFT") {
    start_x <- (active_cols[length(active_cols)] * base_x) - base_x / 2
    end_x <- (active_cols[1] * base_x) - base_x / 2
  }
  if (side == "RIGHT") {
    start_x <- ((active_cols[1] * base_x) - base_x / 2) - base_x
    end_x <- (((active_cols[length(active_cols)] + 1) * base_x) - base_x / 2) +
      base_x
  }
  start <- c(start_x, start_y)
  end <- c(end_x, end_y)

  # make raster
  ## buffer columns
  buffer_col <- rep(0, times = nrow(pf_max_top))
  pf_max_top <- cbind(buffer_col, pf_max_top, buffer_col)

  # line
  toe_line_mat <- shortest_path(mat = pf_max_top, offset_distance, start, end,
                                base_x, base_y, res_scale)

  # return
  return(toe_line_mat)
}


#' @title edge lines
#' @description create edge liens for footprint data
#' @param mat Matrix
#' @param side String
#' @importFrom sf st_as_sfc st_combine st_coordinates st_convex_hull
#' @importFrom utils head tail
#' @noRd

edge_lines <- function(pressure_data, side) {
  # global variables
  x <- y <- x_coord <- y_coord <- me <- sc_df <- mbb <- NULL

  # max pressure image
  max_fp <- pressure_data[[8]]

  # coordinates
  sens_coords <- pressure_data[[7]][, c(1, 2)]
  sens_coords_m <- as.matrix(sens_coords)
  mbb <- getMinBBox(sens_coords_m)

  # find chull
  # Define convex hull, expanding to include all sensors
  df_sf <- sens_coords %>%
    st_as_sf(coords = c("x", "y"))
  fp_chull <- st_convex_hull(st_union(df_sf))
  fp_chull <- st_buffer(fp_chull, 0.0025)

  # find foot angle
  side1 <- mbb[c(1:2), ]
  foot_angle_ <- atan((side1[2, 1] - side1[1, 1]) / (side1[2, 2] - side1[1, 2]))
  foot_ang <- foot_angle_ * 180 / pi

  # polys
  poly_08_up <- cx_poly(sens_coords_m, 0.8, "+Y")
  poly_075_up <- cx_poly(sens_coords, 0.75, "+Y")
  poly_06_down <- cx_poly(sens_coords, 0.6, "-Y")
  poly_055_down <- cx_poly(sens_coords, 0.55, "-Y")
  poly_03_up <- cx_poly(sens_coords, 0.3, "+Y")
  poly_01_down <- cx_poly(sens_coords, 0.1, "-Y")

  # top part medial (0.6 - 0.8)
  top_med_ <- st_difference(fp_chull, poly_06_down)
  top_med <- st_coordinates(st_difference(top_med_, poly_08_up))[, c(1, 2)]

  # top part lateral (0.55 - 0.75)
  top_lat_ <- st_difference(fp_chull, poly_075_up)
  top_lat <- st_coordinates(st_difference(top_lat_, poly_055_down))[, c(1, 2)]

  # bottom part (0.1 - 0.3)
  bot_ <- st_difference(fp_chull, poly_03_up)
  bot <- st_coordinates(st_difference(bot_, poly_01_down))[, c(1, 2)]

  # rotate coords
  top_med_rot <- rot_pts(top_med, foot_ang)
  top_lat_rot <- rot_pts(top_lat, foot_ang)
  bot_rot <- rot_pts(bot, foot_ang)

  # find extremes
  if (side == "RIGHT") {
    ff_med <- rot_pts(top_med_rot[which.min(top_med_rot[, 1]), ], foot_ang * -1)
    ff_lat <- rot_pts(top_lat_rot[which.max(top_lat_rot[, 1]), ], foot_ang * -1)
    heel_med <- rot_pts(bot_rot[which.min(bot_rot[, 1]), ], foot_ang * -1)
    heel_lat <- rot_pts(bot_rot[which.max(bot_rot[, 1]), ], foot_ang * -1)
  }
  if (side == "LEFT") {
    ff_med <- rot_pts(top_med_rot[which.max(top_med_rot[, 1]), ], foot_ang * -1)
    ff_lat <- rot_pts(top_lat_rot[which.min(top_lat_rot[, 1]), ], foot_ang * -1)
    heel_med <- rot_pts(bot_rot[which.max(bot_rot[, 1]), ], foot_ang * -1)
    heel_lat <- rot_pts(bot_rot[which.min(bot_rot[, 1]), ], foot_ang * -1)
  }

  # lines
  med_side_line_ <- st_linestring(rbind(heel_med, ff_med))
  med_side_line <- st_extend_line(med_side_line_, 3)
  lat_side_line_ <- st_linestring(rbind(heel_lat, ff_lat))
  lat_side_line <- st_extend_line(lat_side_line_, 3)

  # return lines
  return(list(med_side_line, lat_side_line))
}


#' poly as a percentage of min bbox
#' @noRd

cx_poly <- function(sens_coords, perc, direction) {
  # bounding box
  sens_coords_m <- as.matrix(sens_coords)
  mbb <- getMinBBox(sens_coords_m)

  ## Bounding box sides
  side1 <- mbb[c(1:2), ]
  side2 <- mbb[c(3:4), ]

  ## Get crossing lines
  side1_ <- side1[1, ] + ((side1[2, ] - side1[1, ]) * perc)
  side2_ <- side2[1, ] + ((side2[2, ] - side2[1, ]) * perc)
  line_ <- st_linestring(as.matrix(rbind(side1_, side2_)))
  line_ <- st_coordinates(st_extend_line(line_, 1))[, c(1:2)]
  poly <- st_line2polygon(line_[c(1, nrow(line_)), ], 1, direction)

  # return
  return(poly)
}


#' Color binning function
#' taken from ggplot2 codebase
#' @noRd

binned_pal <- function(palette) {
  function(x) {
    palette(length(x))
  }
}


#' Rotate line
#' @noRd
#'
rot_line <- function(line, ang, cnt) {
  ang_ <- ang * pi /180
  new_line <- (line - cnt) * matrix(c(cos(ang_), sin(ang_),
                                      -sin(ang_), cos(ang_)), 2, 2) + cnt
  return(new_line)
}

#' rotate points
#' @noRd
#'
rot_pts <- function (pts, ang) {
  radian <- ang * ( pi / 180)
  x <- c(cos(radian), -sin(radian))
  y <- c(sin(radian), cos(radian))
  new_pts = pts %*% cbind(x, y)
  return(new_pts)
}


#' @title icp mask
#' @description match two landmark configurations using iteratively closest
#' point search
#' @param x moving landmarks
#' @param y target landmarks
#' @param iterations integer: number of iterations
#' @param mindist restrict valid points to be within this distance
#' @param subsample use a subsample determined by kmean clusters to speed up
#' computation
#' @param type character: select the transform to be applied, can be "rigid",
#' "similarity" or "affine"
#' @param weights vector of length \code{nrow(x)} containing weights for each
#' row in \code{x}
#' @param centerweight logical: if weights are defined and centerweigths=TRUE,
#' the matrix will be centered according to these weights instead of the
#' barycenter.
#' @param threads integer: number of threads to use.
#' @return returns the rotated landmarks
#' @importFrom Rvcg vcgKDtree vcgSearchKDtree vcgCreateKDtree
#' @importFrom Morpho computeTransform applyTransform
#' @noRd
icp_mask <- function(x, y, mask_list, iterations = 100, mindist = 1e15,
                type = "similarity", threads = 1, centerweight = FALSE) {
  # add column if 2D
  m <- ncol(x)
  if (m == 2) {
    x <- cbind(x,0)
    y <- cbind(y,0)
  }

  # type
  type <- match.arg(type, c("rigid", "similarity", "affine"))

  # temporary moving coords
  xtmp <- x

  # KD tree
  yKD <- vcgCreateKDtree(y)

  # iterate to find best fit
  for (i in 1:iterations) {
    clost <- vcgSearchKDtree(yKD, xtmp, 1, threads = threads)
    good <- which(clost$distance < mindist)
    trafo <- computeTransform(y[clost$index[good],], xtmp[good,], type = type,
                              weights = NULL, centerweight = centerweight)
    xtmp <- applyTransform(xtmp[,], trafo)
  }

  # final transform
  fintrafo <- computeTransform(xtmp, x, type = type)

  # transform mask coordinates
  mask_t <- mask_list
  for (i in 1:length(mask_list)) {
    # mask
    mask_crds <- st_coordinates(mask_list[[i]])[, c(1, 2)]
    mask_crds <- cbind(mask_crds, 0)
    mask_crds_t <- applyTransform(mask_crds, inverse = TRUE, fintrafo)
    mask_crds_t <- rbind(mask_crds_t, mask_crds_t[nrow(mask_crds_t), ])
    mask_t[[i]] <- st_polygon(list(mask_crds_t))
  }

  # return
  return(mask_t)
}


#' @title Align mask
#' @description Align mask, usually to trials from the same person
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param mask List. List with masks to be transformed
#' @return List.
#' @importFrom sf st_coordinates st_polygon st_convex_hull
#' @noRd
align_mask <- function(pressure_data, masks) {
  # get outline of pressure
  outline_coords <- pressure_data$sens_polygons[, c(1, 2)] %>%
    st_as_sf(coords = c("x", "y"))
  outline_chull <- st_convex_hull(st_combine(outline_coords))
  outline_coords_mat <- st_coordinates(outline_chull)[, c(1, 2)]

  # align mask
  mask_coords <- data.frame(x = double(), y = double())
  for (i in 1:length(masks)) {
    crds <- st_coordinates(masks[[i]])[, c(1, 2)]
    mask_coords <- rbind(mask_coords, crds)
  }
  coords_df <- mask_coords %>%
    st_as_sf(coords = c("X", "Y"))
  fp_chull <- st_convex_hull(st_combine(coords_df))
  mask_coords_mat <- st_coordinates(fp_chull)[, c(1, 2)]

  # transform masks
  masks_trans <- icp_mask(outline_coords_mat, mask_coords_mat, masks)

  # return aligned mask
  pressure_data[[5]] <- masks_trans
  return(pressure_data)
}


#' @title Visualize masks
#' @description Visualize the existing masks
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param visual_list List. Mask numbers that want to be viewed. (Default is
#' all exisiting masks)
#' @param image String."max" = footprint of maximum sensors. "mean"
#'   average value of sensors over time
#' @return ggplot plot object
#' @examples
#' \dontrun{
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- automask(pressure_data, foot_side = "auto", plot = TRUE)
#' plot_masks(pressure_data, image = "max")
#' }
#' @noRd

plot_masks <- function(pressure_data,
                       visual_list = seq(1, length(pressure_data[[5]])),
                       image = "max"){
  # global variables
  X <- Y <- NULL

  # initialize plot
  #grDevices::x11()

  if (pressure_data[[2]] != "pedar"){
    x_max <- max(pressure_data[[7]]$x) + 0.01
    y_max <- max(pressure_data[[7]]$y) + 0.01
    g <- plot_pressure(pressure_data, image, plot = FALSE)
    g <- g + scale_x_continuous(expand = c(0, 0), limits = c(-0.01, x_max))
    g <- g + scale_y_continuous(expand = c(0, 0), limits = c(-0.01, y_max))
  } else {
    if (image == "max"){
      g <- plot_pressure(pressure_data, variable = "max", plot = FALSE)
    } else if (is.numeric(image)){
      g <- plot_pressure(pressure_data, frame = image, plot = FALSE)
    } else {
      stop("The variable image is invalid")
    }
  }
  #print(g)

  # plot original mask data
  for (n_mask in visual_list) {
    mask_plt_df <- data.frame(st_coordinates(pressure_data[[5]][[n_mask]])[,1:2])
    g <- g + geom_path(data = mask_plt_df, aes(X, Y), color = "red", linewidth = 1)
    g <- g + geom_point(data = mask_plt_df, aes(X, Y), color = "red", size = 2)
  }
  print(g)

  return(g)
}


#' @title Peak pressure
#' @description get highest value of sensor in mask
#' @noRd

pressure_peak <- function(pressure_data, sens_mask_df, mask) {
  # analysis
  val <- max(apply(pressure_data, 2, max)[which(sens_mask_df[, mask] > 0)])

  # return
  return(val)
}


#' @title Peak pressure time series
#' @description get highest value of sensor in mask per frame
#' @noRd

pressure_peak_ts <- function(pressure_data, sens_mask_df, mask) {
  val_ts <- rep(NA, nrow(pressure_data))
  # analysis
  for (i in 1:nrow(pressure_data)) {
    val_ts[i] <- max(pressure_data[i, which(sens_mask_df[, mask] > 0)])
  }

  # return
  return(val_ts)
}


#' @title Mean regional pressure
#' @description get average pressure across mask. Defined as highest force in
#' mask over maximum contact area
#' @noRd

pressure_mean <- function(pressure_data_, pressure_data, sens_mask_df, mask) {
  # process
  mean_pressure <- rep(NA, nrow(pressure_data_))
  active_sensors <- rep(0, ncol(pressure_data_))
  active_sensors[apply(pressure_data_, 2, max) > 0] <- 1
  total_ca_mask <- sum(active_sensors * sens_mask_df[, mask] * pressure_data[[3]])
  for (i in 1:nrow(pressure_data_)) {
    force <- sum(pressure_data_[i, ] * sens_mask_df[, mask] * pressure_data[[3]])
    mean_pressure[i] <- force / total_ca_mask
  }
  val <- max(mean_pressure)

  # return
  return(val)
}


#' @title Contact area by mask
#' @description get contact areas for masks
#' @noRd

contact_area <- function(pressure_data_, pressure_data, sens_mask_df, mask) {
  # process
  active_area <- rep(NA, length.out = nrow(pressure_data_))
  for (i in 1:nrow(pressure_data_)) {
    active_sensors <- rep(0, length(pressure_data_[i, ]))
    active_sensors[pressure_data_[i, ] > 0] <- 1
    active_area[i] <- sum(active_sensors * sens_mask_df[, mask] * pressure_data[[3]])
  }
  val <- max(active_area)

  # return
  return(val)
}


#' @title Pressure time integral (Novel)
#' @description calculate PTI (as defined by novel) for pressure data. Def:
#' Highest reading from any sensor in mask at each time point, then find area
#' under curve
#' @noRd

pti_1 <- function(pressure_data_, pressure_data, sens_mask_df, mask) {
  # process
  mask_ppt <- rep(NA, length.out = nrow(pressure_data_))
  for (i in 1:(nrow(pressure_data_))) {
    mask_ppt[i] <- max(pressure_data_[i, ][which(sens_mask_df[, mask] > 0)]) *
      pressure_data[[4]]
    }
  val <- sum(mask_ppt)

  # return
  return(val)
}


#' @title Pressure time integral (Melai)
#' @description calculate PTI (as defined by Melai) for pressure data
#' @noRd

pti_2 <- function(pressure_data_, pressure_data, sens_mask_df, mask) {
  # analysis
  force <- rep(NA, length.out = nrow(pressure_data_))
  active_sensors <- rep(0, ncol(pressure_data_))
  active_sensors[apply(pressure_data_, 2, max) > 0] <- 1
  total_ca_mask <- sum(active_sensors * sens_mask_df[, mask] * pressure_data[[3]])
  for (i in 1:nrow(pressure_data_)) {
    P <- pressure_data_[i, ] * pressure_data[[3]] * 1000
    force[i] <- sum(P * sens_mask_df[, mask])
  }
  val <- pracma::trapz(seq(0, by = pressure_data[[4]],
                           length.out = nrow(pressure_data_)),
                       force) / total_ca_mask / 1000

  # return
  return(val)
}


#' @title Force time integral
#' @description Calculate FTI
#' @noRd

fti <- function(pressure_data_, pressure_data, sens_mask_df, mask) {
  # analysis
  force <- rep(NA, length.out = nrow(pressure_data_))
  for (i in 1:(nrow(pressure_data_))) {
    P <- pressure_data_[i, ] * pressure_data[[3]] * 1000
    force[i] <- sum(P * sens_mask_df[, mask])
  }
  val <- pracma::trapz(seq(0, by = pressure_data[[4]],
                           length.out = nrow(pressure_data_)), force)

  # return
  return(val)
}


#' @title Force peak
#' @description Calculate force peak
#' @noRd

force_peak <- function(pressure_data_, pressure_data, sens_mask_df, mask) {
  # calculate peak force
  force <- rep(NA, length.out = nrow(pressure_data_))
  for (i in 1:nrow(pressure_data_)) {
    P <- pressure_data_[i, ] * pressure_data[[3]] * 1000
    force[i] <- sum(P * sens_mask_df[, mask])
  }
  val <- max(force)

  # return
  return(val)
}

#' @title Force time series
#' @description Calculate force per frame
#' @noRd

force_ts <- function(pressure_data_, pressure_data, sens_mask_df, mask) {
  force <- rep(NA, nrow(pressure_data_))
  for (i in 1:nrow(pressure_data_)) {
    P <- pressure_data_[i, ] * pressure_data[[3]] * 1000
    force[i] <- sum(P * sens_mask_df[, mask])
  }

  # return
  return(force)
}


#' @title threshold event
#' @description defines threshold based on force value
#' @noRd
threshold_event <- function(pressure_data, threshold = "auto", min_frames, side) {
  # make force vector
  if (pressure_data[[2]] != "pedar") {
    force <- whole_pressure_curve(pressure_data, "force")
  } else {
    if (side == "RIGHT") {force <- force_pedar(pressure_data)[, 1]}
    if (side == "LEFT") {force <- force_pedar(pressure_data)[, 2]}
  }

  # Adjust thresholds to avoid errors
  min_f <- min(force)
  max_f <- max(force)
  if (threshold == "auto") {
    thresh <- min_f + (0.1 * (max_f - min_f))
  } else {
    if (is.numeric(threshold) == TRUE) {
      thresh <- threshold + 0.01
    } else {
        stop("threshold can be auto or a numeric value")
    }
  }

  # throw error if threshold is less than min of trial
  if (min_f > thresh) {
    stop("threshold is less than minimum force recorded in trial")
  }

  # Get events
  FS_events <- which(force[-length(force)] < thresh &
                       force[-1] > thresh) + 1
  FO_events <- which(force[-length(force)] > thresh &
                       force[-1] < thresh) + 1

  # Create steps
  ## Adjust events to ensure stance phase is used
  if (FO_events[1] < FS_events[1]) {FO_events <- FO_events[-1]}
  if ((length(FS_events)) > (length(FO_events))) {
    FS_events <- FS_events[-length(FS_events)]
  }

  # remove short steps
  rm_stp <- c()
  for (i in 1:length(FS_events)) {
    if (FO_events[i] - FS_events[i] < min_frames) {rm_stp <- c(rm_stp, i)}
  }
  if (length(rm_stp > 0)) {
    FS_events <- FS_events[-rm_stp]
    FO_events <- FO_events[-rm_stp]
  }

  # Make steps into data frame
  df <- data.frame(step = integer(), frame = integer(), force = double())
  for (i in 1:length(FS_events)) {
    force_step <- force[FS_events[i]:FO_events[i]]
    step <- data.frame(step = rep(i, length_out = length(force_step)),
                       frame = c(1:length(force_step)),
                       force = force_step)
    df <- rbind(df, step)
  }

  # return
  return(list(df, FS_events, FO_events))
}


#' @title Sensor centroid
#' @description Calculate centroid cooridinates for sensor
#' @noRd
sensor_centroid <- function(pressure_data) {
  id <- NULL
  sensors <- unique(pressure_data[[7]][, 3])
  centroids <- data.frame(x = rep(NA, length(sensors)),
                          y = rep(NA, length(sensors)))
  for (sens in seq_along(sensors)) {
    mat <- pressure_data[[7]] %>% filter(id == sensors[sens])
    centroids[sens, ] <- colMeans(mat[1:(nrow(mat) - 1), c(1, 2)])
  }
  return(centroids)
}


#' @title automask footprint
#' @description Automatically creates mask of footprint data
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement. z dimension represents time
#' @param mask_scheme String.
#' @param res_scale Numeric vector. Defines low cost functions for lines
#' @param foot_side String. "RIGHT", "LEFT", or "auto". Auto uses
#' auto_detect_side function
#' @param plot Logical. Whether to play the animation
#' @return List. Masks are added to pressure data variable
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#' @importFrom zoo rollapply
#' @importFrom sf st_union st_difference st_intersection
#' @noRd

automask <- function(pressure_data, mask_scheme, res_scale, foot_side = "auto",
                     plot = TRUE) {
  # check data isn't from pedar
  if (pressure_data[[2]] == "pedar")
    stop("automask doesn't work with pedar data")

  # global variables
  x <- y <- mask <- x_coord <- y_coord <- me <- heel_cut_dist <-
    mfoot_cut_prox <- ffoot_cut_prox <- NULL

  # Find footprint (max)
  max_df <- pressure_data[[8]]

  # coordinates
  sens_coords <- pressure_data[[7]][, c(1, 2)]

  # side
  if (foot_side == "auto") {
    side <- auto_detect_side(pressure_data)
  } else {
    side <- foot_side
  }

  # Define minimum bounding box
  sens_coords_m <- as.matrix(sens_coords)
  mbb <- getMinBBox(sens_coords_m)

  # Define convex hull, expanding to include all sensors
  df_sf <- sens_coords %>%
    st_as_sf(coords = c("x", "y"))
  fp_chull <- st_convex_hull(st_union(df_sf))
  fp_chull <- st_buffer(fp_chull, 0.0025)

  # Simple 3 mask (forefoot, midfoot, and hindfoot)
  ## Bounding box sides
  side1 <- mbb[c(1:2), ]
  side2 <- mbb[c(3:4), ]

  ## Get distal 27% and 55% lines that divide into masks
  ### heel
  side1_27 <- side1[1, ] + ((side1[2, ] - side1[1, ]) * 0.27)
  side2_27 <- side2[1, ] + ((side2[2, ] - side2[1, ]) * 0.27)
  line_27 <- st_linestring(as.matrix(rbind(side1_27, side2_27)))
  line_27 <- st_extend_line(line_27, 1)
  line_27_ints <- st_coordinates(st_intersection(fp_chull, line_27))
  heel_line <- shortest_path(max_df, 0, line_27_ints[1, c(1:2)],
                             line_27_ints[2, c(1:2)], 0.005, 0.005,
                             res_scale[1])
  heel_line_dist_poly <- st_line2polygon(st_linestring(heel_line), 1, "+Y")
  heel_line_prox_poly <- st_line2polygon(st_linestring(heel_line), 1, "-Y")

  ### midfoot/forefoot
  side1_55 <- side1[1, ] + ((side1[2, ] - side1[1, ]) * 0.55)
  side2_55 <- side2[1, ] + ((side2[2, ] - side2[1, ]) * 0.55)
  line_55 <- st_linestring(as.matrix(rbind(side1_55, side2_55)))
  line_55 <- st_extend_line(line_55, 1)
  line_55_ints <- st_coordinates(st_intersection(fp_chull, line_55))
  midff_line <- shortest_path(max_df, 0, line_55_ints[1, c(1:2)],
                              line_55_ints[2, c(1:2)], 0.005, 0.005,
                              res_scale[2])
  midff_line_dist_poly <- st_line2polygon(st_linestring(midff_line), 1, "+Y")
  midff_line_prox_poly <- st_line2polygon(st_linestring(midff_line), 1, "-Y")

  ## forefoot, midfoot, and hindfoot masks
  hf_mask <- st_difference(fp_chull, heel_line_dist_poly)
  mf_mask <- st_difference(fp_chull, midff_line_dist_poly)
  mf_mask <- st_difference(mf_mask, heel_line_prox_poly)
  ff_mask <- st_difference(fp_chull, midff_line_prox_poly)

  # if more complex automask required
  if (mask_scheme == "automask_simple") {
    # Make mask list
    mask_list <- list(heel_mask = hf_mask, midfoot_mask = mf_mask,
                      forefoot_mask = ff_mask)
  } else if (mask_scheme == "automask_novel") {
    # Define angles for dividing lines between metatarsals
    ## get edge lines
    edges <- edge_lines(pressure_data, side)

    ## angle between lines
    med_pts <- st_coordinates(edges[[1]])
    lat_pts <- st_coordinates(edges[[2]])
    med_line_angle <- (atan((med_pts[2, 1] - med_pts[1, 1]) /
                              (med_pts[2, 2] - med_pts[1, 2]))) * 180 / pi
    lat_line_angle <- (atan((lat_pts[2, 1] - lat_pts[1, 1]) /
                              (lat_pts[2, 2] - lat_pts[1, 2]))) * 180 / pi
    alpha <- lat_line_angle - med_line_angle

    ## intersection point
    med_lat_int <- st_intersection(edges[[1]], edges[[2]])

    ## rotation angles for 2-4 lines
    MT_hx_alpha <- (0.33 * alpha)
    MT_12_alpha <- (0.3  * alpha)
    MT_23_alpha <- (0.47 * alpha)
    MT_34_alpha <- (0.64 * alpha)
    MT_45_alpha <- (0.81 * alpha)

    ## hallux and toe angle lines
    MT_hx_line <- rot_line(edges[[1]], MT_hx_alpha, med_lat_int)
    MT_12_line <- rot_line(edges[[1]], MT_12_alpha, med_lat_int)
    MT_23_line <- rot_line(edges[[1]], MT_23_alpha, med_lat_int)
    MT_34_line <- rot_line(edges[[1]], MT_34_alpha, med_lat_int)
    MT_45_line <- rot_line(edges[[1]], MT_45_alpha, med_lat_int)

    ## simplify toe line
    ### toe line
    toe_line_mat <- toe_line(pressure_data, side, res_scale[3])
    HX_toe_pt <- st_intersection(st_linestring(toe_line_mat), MT_hx_line)
    MT12_toe_pt <- st_intersection(st_linestring(toe_line_mat), MT_12_line)
    MT23_toe_pt <- st_intersection(st_linestring(toe_line_mat), MT_23_line)
    MT34_toe_pt <- st_intersection(st_linestring(toe_line_mat), MT_34_line)
    MT45_toe_pt <- st_intersection(st_linestring(toe_line_mat), MT_45_line)
    MT5_toe_pt <- st_intersection(st_linestring(toe_line_mat), edges[[2]])
    toe_line_mat_simple <- rbind(MT5_toe_pt, MT45_toe_pt, MT34_toe_pt,
                                 MT23_toe_pt, HX_toe_pt, MT12_toe_pt)
    if (side == "LEFT") {
      toe_line_mat_simple <- rbind(toe_line_mat[nrow(toe_line_mat), ],
                                   toe_line_mat_simple, toe_line_mat[1, ])
    }
    if (side == "RIGHT") {
      toe_line_mat_simple <- rbind(toe_line_mat[nrow(toe_line_mat), ],
                                   toe_line_mat_simple, toe_line_mat[1, ])
    }

    ## hallux line
    hx_int_top <- st_coordinates(st_intersection(fp_chull, MT_hx_line))
    hx_int_top_ <- unname(hx_int_top[which.max(hx_int_top[, 2]), c(1, 2)])
    hx_int_top_[2] <- hx_int_top_[2] - 0.003
    HX_toe_pt_ <- c(st_coordinates(HX_toe_pt))
    hx_l <- shortest_path(max_df, 0, HX_toe_pt_,
                          hx_int_top_, 0.005, 0.005, res_scale[4],
                          extend_line = FALSE)
    hx_line <- rbind(c(hx_l[1, 1], hx_l[1, 2] - 1), hx_l,
                     c(hx_l[nrow(hx_l), 1], hx_l[nrow(hx_l), 2] + 1))

    ## polys for met and hal cuts
    if (side == "RIGHT") {lat <- "+X"; med <- "-X"}
    if (side == "LEFT") {lat <- "-X"; med <- "+X"}
    MT_hx_poly_lat <- st_line2polygon(hx_line, 1, lat)
    MT_hx_poly_med <- st_line2polygon(hx_line, 1, med)
    MT_12_poly_lat <- st_line2polygon(st_coordinates(MT_12_line)[, 1:2], 1, lat)
    MT_12_poly_med <- st_line2polygon(st_coordinates(MT_12_line)[, 1:2], 1, med)
    MT_23_poly_lat <- st_line2polygon(st_coordinates(MT_23_line)[, 1:2], 1, lat)
    MT_23_poly_med <- st_line2polygon(st_coordinates(MT_23_line)[, 1:2], 1, med)
    MT_34_poly_lat <- st_line2polygon(st_coordinates(MT_34_line)[, 1:2], 1, lat)
    MT_34_poly_med <- st_line2polygon(st_coordinates(MT_34_line)[, 1:2], 1, med)
    MT_45_poly_lat <- st_line2polygon(st_coordinates(MT_45_line)[, 1:2], 1, lat)
    MT_45_poly_med <- st_line2polygon(st_coordinates(MT_45_line)[, 1:2], 1, med)

    ## toe cut polygons
    toe_poly_prox <- st_line2polygon(as.matrix(toe_line_mat_simple), 1, "-Y")
    toe_poly_dist <- st_line2polygon(as.matrix(toe_line_mat_simple), 1, "+Y")

    ## make masks
    toe_mask <- st_difference(fp_chull, toe_poly_prox)
    ff_mask <- st_difference(ff_mask, toe_poly_dist)

    ## make met masks
    hal_mask <- st_difference(toe_mask, MT_hx_poly_lat)
    l_toe_mask <- st_difference(toe_mask, MT_hx_poly_med)
    MT1_mask <- st_difference(ff_mask, MT_12_poly_lat)
    MT2_mask <- st_difference(ff_mask, MT_23_poly_lat)
    MT2_mask <- st_difference(MT2_mask, MT_12_poly_med)
    MT3_mask <- st_difference(ff_mask, MT_34_poly_lat)
    MT3_mask <- st_difference(MT3_mask, MT_23_poly_med)
    MT4_mask <- st_difference(ff_mask, MT_45_poly_lat)
    MT4_mask <- st_difference(MT4_mask, MT_34_poly_med)
    MT5_mask <- st_difference(ff_mask, MT_45_poly_med)

    # Make mask list
    mask_list <- list(heel_mask = hf_mask,
                      midfoot_mask = mf_mask,
                      forefoot_mask = ff_mask,
                      hal_mask = hal_mask,
                      l_toe_mask = l_toe_mask,
                      MT1_mask = MT1_mask,
                      MT2_mask = MT2_mask,
                      MT3_mask = MT3_mask,
                      MT4_mask = MT4_mask,
                      MT5_mask = MT5_mask)
  }

  # plot footprint and masks if required
  if (plot == TRUE) {
    # make masks into df format
    mask_df <- masks_2_df(mask_list)

    # Plot footprint and masks
    x_max <- max(pressure_data[[7]]$x) + 0.01
    y_max <- max(pressure_data[[7]]$y) + 0.01
    g <- plot_pressure(pressure_data, "max", plot = FALSE)
    g <- g + scale_x_continuous(expand = c(0, 0), limits = c(-0.01, x_max))
    g <- g + scale_y_continuous(expand = c(0, 0), limits = c(-0.01, y_max))
    g <- g + geom_path(data = mask_df, aes(x = x, y = y, group = mask),
                       color = "red", linewidth = 1)
    suppressMessages(print(g))
  }

  # Return masks for analysis
  pressure_data[[5]] <- mask_list
  return(pressure_data)
}


#' @title approve step
#' @description asks user to approve step
#' @noRd
approve_step <- function(df, on_v, off_v, n_steps, side = "none") {
  # global variables
  step <- frame <- NULL

  # Approve or discard steps
  include_stps <- c()
  for (stp in 1:length(on_v)) {
    # highlighted step df
    df1 <- df %>% filter(step %in% stp)

    # plot
    g <- ggplot()
    g <- g + geom_line(data = df, aes(x = frame, y = force, group = step),
                       linewidth = 1.5, color = "grey")
    g <- g + geom_line(data = df1, aes(x = frame, y = force),
                       linewidth = 1.5, color = "red")
    g <- g + xlab("Frame no") + ylab("Force (N)")
    g <- g + ggtitle(paste0(side, " ", stp))
    print(g)
    Sys.sleep(1.0)

    # get user to approve or reject step
    resp <- menu(c("Y", "N", "EXIT"),
                 title = "Do you want to keep (Y) or discard (N) this step?")
    if (resp == 3) {
      break
    } else {
      include_stps <- c(include_stps, resp)
    }

    # check if number of steps has been reached
    if (sum(include_stps == 1) >= n_steps) {break}
  }

  # make events df
  len_ <- length(which(include_stps == 1))
  event_df <- data.frame(side = rep(side, length.out = len_),
                         FON = on_v[which(include_stps == 1)],
                         FOFF = off_v[which(include_stps == 1)])

  # return
  return(event_df)
}


#' @title Pedar mask regions
#' @description Creates a mask region from pedar sensel coordinates
#' @param pressure_data List.
#' @param position Dataframe. A n x 3 dataframe of sensel coordinates
#' [sensel id, x, y]
#' @param sensel_list List. List of sensels to include
#' @return polygon. A mask polygon
#' @noRd

pedar_polygon <- function(pressure_data, sensel_list) {

  #polygon_list <- sensor_2_polygon(pressure_data, pressure_image = "all_active",
  #                                 frame = NA, output = "sf")
  sensor_polygons <- sens_df_2_polygon(pressure_data[[7]])
  # Left foot sensels are stored as 1:99, right foot senses are 101:99
  #if (foot_side == "RIGHT"){
  #  sensel_list <- sensel_list + 99
  #}
  sensel_polygon <- st_union(st_sfc(sensor_polygons[sensel_list]))
  if (length(st_geometry(sensel_polygon)[[1]]) > 1)
    stop("sensels must share a corner or an edge")
  #sensel_polygon <- polygon_list[[sensel_list[1]]]

  #for (sensel_idx in sensel_list[-1]) {
  #  if (class(st_intersection(sensel_polygon, polygon_list[[sensel_idx]]))[[2]]
  #      %in% c("LINESTRING", "POLYGON", "MULTILINESTRING")) {
  #    sensel_polygon <- st_union(sensel_polygon, polygon_list[[sensel_idx]])
  #  } else {
  #    stop("Sensels used to define a mask must share an edge")
  #  }
  #}

  return(sensel_polygon)
}


#' @title pedar_mask1
#' @description applies pedar mask 1
#' @noRd
pedar_mask1 <- function(pressure_data) {
  # define masks by sensel number
  hindfoot_L <- pedar_polygon(pressure_data, 1:26)
  midfoot_L <-pedar_polygon(pressure_data, 27:54)
  forefoot_L <- pedar_polygon(pressure_data, 55:82)
  toes_L <- pedar_polygon(pressure_data, 83:99)

  hindfoot_R <- pedar_polygon(pressure_data, 100:125)
  midfoot_R <-pedar_polygon(pressure_data, 126:153)
  forefoot_R <- pedar_polygon(pressure_data, 154:181)
  toes_R <- pedar_polygon(pressure_data, 182:198)

  mask_list <- list(L_hindfoot_mask = hindfoot_L,
                    L_midfoot_mask = midfoot_L,
                    L_forefoot_mask = forefoot_L,
                    L_toes_mask = toes_L,
                    R_hindfoot_mask = hindfoot_R,
                    R_midfoot_mask = midfoot_R,
                    R_forefoot_mask = forefoot_R,
                    R_toes_mask = toes_R)

  # return
  return(mask_list)
}


#' @title pedar_mask2
#' @description applies pedar mask 2
#' @noRd
pedar_mask2 <- function(pressure_data) {
  # foot outline and bounding box
  bbox_L <- sf::st_bbox(pedar_polygon(pressure_data, 1:99))
  outline_mask <- pedar_polygon(pressure_data, 1:99)

  # define percent lines
  hindfoot_line_Ly <- (bbox_L$ymax - bbox_L$ymin) / 2 + bbox_L$ymin
  forefoot_line_Ly <- (bbox_L$ymax - bbox_L$ymin) * 0.85 + bbox_L$ymin
  hallux_line_Lx <- bbox_L$xmax - (bbox_L$xmax - bbox_L$xmin) * 0.35

  # define masks
  hindfoot_line_L <- data.frame(x_coord = c(bbox_L$xmin, bbox_L$xmax),
                                y_coord =  c(hindfoot_line_Ly, hindfoot_line_Ly))
  hindfoot_line_L <- st_extend_line(st_linestring(as.matrix(hindfoot_line_L)), 1)
  hindfoot_line_L_dist <- st_line2polygon(hindfoot_line_L, 1, "+Y")
  hindfoot_L <- st_difference(outline_mask, hindfoot_line_L_dist)

  forefoot_line_L <- data.frame(x_coord = c(bbox_L$xmin,bbox_L$xmax),
                                y_coord =  c(forefoot_line_Ly,forefoot_line_Ly))
  forefoot_line_L <- st_extend_line(st_linestring(as.matrix(forefoot_line_L)),1)
  forefoot_line_L_dist <- st_line2polygon(forefoot_line_L, 1, "+Y")
  forefoot_line_L_prox <- st_line2polygon(forefoot_line_L, 1, "-Y")
  forefoot_L <- st_difference(outline_mask, forefoot_line_L_dist)

  toes_L <- st_difference(outline_mask, forefoot_line_L_prox)
  hallux_line_L <- data.frame(x_coord = c(hallux_line_Lx,hallux_line_Lx),
                              y_coord =  c(bbox_L$ymin,bbox_L$ymax))
  hallux_line_L <- st_extend_line(st_linestring(as.matrix(hallux_line_L)),1)
  hallux_line_L_lat <- st_line2polygon(hallux_line_L, 1, "-X")
  hallux_line_L_med <- st_line2polygon(hallux_line_L, 1, "+X")
  hallux_L <- st_difference(toes_L, hallux_line_L_lat)
  lesser_toes_L <- st_difference(toes_L, hallux_line_L_med)

  bbox_R <- sf::st_bbox(pedar_polygon(pressure_data, 100:198))
  outline_mask <- pedar_polygon(pressure_data, 100:198)

  hindfoot_line_Ry <- (bbox_R$ymax - bbox_R$ymin) / 2 + bbox_R$ymin
  forefoot_line_Ry <- (bbox_R$ymax - bbox_R$ymin) * 0.85 + bbox_R$ymin
  hallux_line_Rx <- bbox_R$xmin + (bbox_R$xmax - bbox_R$xmin) * 0.35

  hindfoot_line_R <- data.frame(x_coord = c(bbox_R$xmin, bbox_R$xmax),
                                y_coord =  c(hindfoot_line_Ry, hindfoot_line_Ry))
  hindfoot_line_R <- st_extend_line(st_linestring(as.matrix(hindfoot_line_R)), 1)
  hindfoot_line_R_dist <- st_line2polygon(hindfoot_line_R, 1, "+Y")
  hindfoot_R <- st_difference(outline_mask, hindfoot_line_R_dist)

  forefoot_line_R <- data.frame(x_coord = c(bbox_R$xmin, bbox_R$xmax),
                                y_coord =  c(forefoot_line_Ry, forefoot_line_Ry))
  forefoot_line_R <- st_extend_line(st_linestring(as.matrix(forefoot_line_R)), 1)
  forefoot_line_R_dist <- st_line2polygon(forefoot_line_R, 1, "+Y")
  forefoot_line_R_prox <- st_line2polygon(forefoot_line_R, 1, "-Y")
  forefoot_R <- st_difference(outline_mask, forefoot_line_R_dist)

  toes_R <- st_difference(outline_mask, forefoot_line_R_prox)
  hallux_line_R <- data.frame(x_coord = c(hallux_line_Rx, hallux_line_Rx),
                              y_coord =  c(bbox_R$ymin, bbox_R$ymax))
  hallux_line_R <- st_extend_line(st_linestring(as.matrix(hallux_line_R)), 1)
  hallux_line_R_lat <- st_line2polygon(hallux_line_R, 1, "+X")
  hallux_line_R_med <- st_line2polygon(hallux_line_R, 1, "-X")
  hallux_R <- st_difference(toes_R, hallux_line_R_lat)
  lesser_toes_R <- st_difference(toes_R, hallux_line_R_med)

  mask_list <- list(L_hindfoot_mask = hindfoot_L,
                    L_forefoot_mask = forefoot_L,
                    L_hallux_mask = hallux_L,
                    L_lesser_toes_mask = lesser_toes_L,
                    R_hindfoot_mask = hindfoot_R,
                    R_forefoot_mask = forefoot_R,
                    R_hallux_mask = hallux_R,
                    R_lesser_toes_mask = lesser_toes_R)
}


#' @title pedar_mask3
#' @description applies pedar mask 3
#' @noRd
pedar_mask3 <- function(pressure_data) {
  # define masks by sensel numbers
  med_rf_L <- pedar_polygon(pressure_data, c(1:2, 6:8, 13:15, 20:22))
  lat_rf_L <-pedar_polygon(pressure_data, c(9:12, 16:19, 23:26))
  med_mf_L <- pedar_polygon(pressure_data, c(27:29, 34:36, 41:43, 48:50, 55:57))
  lat_mf_L <- pedar_polygon(pressure_data, c(30:33, 37:40, 44:47, 51:54, 58:59))
  MTPJ1_L <- pedar_polygon(pressure_data, c(62:63, 69:70, 76:77))
  MTPJ23_L <- pedar_polygon(pressure_data, c(64:66, 71:73, 78:80))
  MTPJ45_L <- pedar_polygon(pressure_data, c(60:61, 67:68, 74:75))
  hallux_L <- pedar_polygon(pressure_data, c(83:84, 90:91, 96))
  lesser_toes_L <- pedar_polygon(pressure_data, c(81:82, 85:89, 92:95, 97:99))

  med_rf_R <- pedar_polygon(pressure_data, c(100:101, 105:107, 112:114,
                                             119:121))
  lat_rf_R <-pedar_polygon(pressure_data, c(108:111, 115:118, 122:125))
  med_mf_R <- pedar_polygon(pressure_data, c(126:128, 133:135, 140:142, 147:149,
                                             154:156))
  lat_mf_R <- pedar_polygon(pressure_data, c(129:132, 136:139, 143:146, 150:153,
                                             157:158))
  MTPJ1_R <- pedar_polygon(pressure_data, c(161, 162, 168, 169, 175, 176))
  MTPJ23_R <- pedar_polygon(pressure_data, c(163:165, 170:172, 177:179))
  MTPJ45_R <- pedar_polygon(pressure_data, c(159, 160, 166, 167, 173, 174))
  hallux_R <- pedar_polygon(pressure_data, c(182, 183, 189, 190, 195))
  lesser_toes_R <- pedar_polygon(pressure_data, c(180:181, 184:188, 191:194,
                                                  196:198))

  mask_list <- list(L_medial_hindfoot_mask = med_rf_L,
                    L_lateral_hindfoot_mask = lat_rf_L,
                    L_medial_midfoot_mask = med_mf_L,
                    L_lateral_midfoot_mask = lat_mf_L,
                    L_MTPJ1_mask = MTPJ1_L,
                    L_MTPJ23_mask = MTPJ23_L,
                    L_MTPJ45_mask = MTPJ45_L,
                    L_hallux_mask = hallux_L,
                    L_lesser_toes_mask = lesser_toes_L,
                    R_medial_hindfoot_mask = med_rf_R,
                    R_lateral_hindfoot_mask = lat_rf_R,
                    R_medial_midfoot_mask = med_mf_R,
                    R_lateral_midfoot_mask = lat_mf_R,
                    R_MTPJ1_mask = MTPJ1_R,
                    R_MTPJ23_mask = MTPJ23_R,
                    R_MTPJ45_mask = MTPJ45_R,
                    R_hallux_mask = hallux_R,
                    R_lesser_toes_mask = lesser_toes_R)
}



#' @title Dynamic Plantar Loading Index
#' @description Determine Dynamic plantar loading index (DPLI) for
#' footprint pressure data
#' @param pressure_data List. First item is a 2D array covering each timepoint
#' of the measurement.
#' @param n_bins Numeric. Number of bins used to calculate DPLI
#' @returns Data frame. Contains DPLI for each step
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- create_mask_auto(pressure_data, "automask_novel")
#' pressure_data_ <- pressure_data[[1]]
#' dpli(pressure_data_, pressure_data, sens_mask_df, 1)
#' @importFrom graphics hist
#' @importFrom stats dnorm lm sd
#' @noRd

dpli <- function(pressure_data_, pressure_data, sens_mask_df, mask,
                 n_bins = 30) {
  # row_max
  press_mat <- pressure_data_[, which(sens_mask_df[mask] > 0)]
  press_vec <- apply(press_mat, 1, max)

  # remove empty values from begining and end
  press_vec <- press_vec[min(which(press_vec != 0)):max(which(press_vec != 0))]

  # interpolate to 101
  press_vec <- approxP(press_vec, 101)

  # bin data
  binned <- hist(press_vec, breaks = seq(from = min(press_vec),
                                         to = max(press_vec),
                                         length.out = n_bins + 1),
                 plot = FALSE)

  # get distribution
  sd_press <- sd(press_vec)
  mean_press <- mean(press_vec)
  norm_dist <- dnorm(binned$mids, mean_press, sd_press)
  norm_dist <- norm_dist * (max(binned$counts) / max(norm_dist))

  # get r^2 (dpli)
  val <- summary(lm(binned$counts ~ norm_dist))$adj.r.squared

  # return
  return(val)
}


# =============================================================================

# data
pedar_insole_area <- function() {
  pedar_insole_areas <- data.frame(
    a = c(239, 246, 246, 247, 245, 238, 241, 241, 238, 241, 239, 239, 245, 246,
          246, 243, 246, 246, 245, 240, 239, 240, 238, 240, 239, 240, 244, 243,
          244, 243, 244, 234, 244, 248, 249, 249, 248, 249, 248, 249, 239, 240,
          239, 239, 240, 239, 240, 241, 241, 240, 242, 241, 241, 241, 245, 244,
          245, 245, 245, 244, 245, 242, 241, 242, 241, 244, 241, 242, 246, 247,
          246, 247, 247, 246, 247, 243, 245, 244, 244, 244, 244, 245, 239, 239,
          239, 239, 239, 239, 239, 234, 233, 232, 232, 232, 234, 240, 241, 243,
          242),
    u = c(126, 129, 129, 129, 129, 126, 126, 125, 126, 126, 126, 126, 129, 129,
          129, 129, 129, 129, 129, 126, 126, 127, 126, 126, 126, 127, 127, 127,
          128, 127, 127, 127, 127, 131, 130, 131, 130, 131, 131, 130, 126, 126,
          125, 126, 126, 126, 126, 126, 126, 125, 126, 126, 126, 126, 129, 129,
          129, 129, 129, 129, 128, 127, 127, 128, 127, 127, 127, 127, 129, 129,
          130, 129, 129, 129, 129, 128, 127, 128, 127, 128, 127, 128, 126, 126,
          126, 126, 126, 125, 126, 123, 122, 123, 123, 123, 123, 126, 126, 127,
          127),
    v = c(139, 143, 143, 142, 143, 139, 139, 139, 139, 139, 139, 139, 143, 143,
          143, 143, 143, 143, 143, 139, 139, 139, 139, 139, 139, 139, 141, 141,
          141, 141, 141, 141, 141, 144, 144, 144, 144, 144, 144, 144, 139, 139,
          139, 139, 139, 139, 139, 140, 140, 140, 140, 140, 140, 140, 143, 143,
          143, 143, 143, 143, 143, 140, 140, 140, 140, 140, 140, 140, 144, 144,
          144, 144, 144, 144, 144, 141, 141, 141, 141, 141, 141, 141, 139, 139,
          139, 139, 139, 139, 139, 136, 136, 136, 136, 136, 136, 140, 140, 141,
          140),
    w = c(155, 160, 159, 159, 159, 155, 155, 155, 155, 156, 155, 155, 158, 160,
          158, 158, 160, 158, 160, 155, 156, 155, 155, 156, 154, 156, 158, 159,
          158, 159, 158, 158, 158, 160, 161, 160, 161, 160, 161, 160, 155, 155,
          155, 155, 155, 155, 155, 157, 157, 156, 157, 157, 157, 157, 158, 159,
          158, 159, 159, 159, 159, 156, 157, 158, 157, 157, 157, 157, 160, 159,
          160, 160, 159, 160, 159, 159, 158, 158, 158, 158, 158, 158, 155, 155,
          155, 155, 155, 155, 155, 151, 152, 151, 151, 151, 151, 156, 156, 157,
          156),
    x = c(172, 176, 176, 176, 176, 172, 172, 172, 172, 172, 172, 172, 176, 176,
          176, 176, 176, 176, 176, 172, 172, 172, 172, 172, 172, 172, 175, 175,
          175, 175, 175, 175, 175, 178, 178, 178, 178, 178, 178, 178, 172, 172,
          172, 172, 172, 172, 172, 173, 173, 173, 173, 173, 173, 173, 176, 176,
          176, 176, 176, 176, 176, 173, 173, 173, 173, 173, 173, 173, 178, 178,
          178, 178, 178, 178, 178, 174, 174, 174, 174, 174, 174, 174, 171, 171,
          171, 171, 171, 171, 171, 168, 168, 168, 168, 168, 168, 173, 173, 174,
          174),
    y = c(192, 195, 197, 197, 196, 191, 191, 191, 191, 191, 192, 191, 196, 196,
          196, 196, 196, 197, 195, 192, 190, 189, 192, 191, 191, 189, 195, 193,
          194, 195, 194, 195, 193, 197, 197, 197, 197, 197, 199, 196, 191, 191,
          190, 191, 191, 191, 190, 192, 192, 192, 192, 191, 192, 192, 196, 196,
          196, 196, 196, 196, 196, 191, 193, 192, 191, 193, 192, 192, 197, 198,
          197, 197, 197, 197, 198, 195, 194, 194, 194, 193, 194, 194, 191, 190,
          191, 190, 191, 190, 190, 187, 186, 187, 186, 187, 186, 192, 191, 193,
          192),
    z = c(214, 221, 220, 221, 219, 213, 216, 216, 213, 216, 214, 214, 219, 220,
          220, 218, 220, 220, 219, 215, 214, 215, 213, 215, 214, 215, 219, 217,
          219, 217, 218, 217, 219, 222, 223, 223, 222, 223, 222, 223, 214, 215,
          214, 214, 215, 214, 215, 216, 216, 215, 217, 216, 216, 216, 219, 218,
          219, 219, 219, 218, 219, 217, 216, 217, 216, 218, 216, 217, 221, 222,
          221, 222, 222, 221, 222, 218, 219, 218, 219, 218, 218, 219, 214, 214,
          214, 214, 214, 214, 214, 209, 208, 208, 208, 208, 209, 215, 216, 217,
          217),
    uw = c(114, 101, 107, 101, 114, 103, 115, 114, 115, 114, 113, 105, 121, 123,
           122, 122, 122, 121, 123, 127, 129, 128, 127, 128, 128, 128, 135, 136,
           136, 135, 136, 137, 135, 141, 142, 141, 142, 141, 143, 141, 149, 150,
           149, 150, 149, 150, 148, 155, 155, 156, 155, 155, 156, 154, 156, 159,
           162, 162, 164, 167, 167, 160, 164, 165, 168, 171, 173, 174, 163, 166,
           167, 170, 170, 173, 176, 158, 160, 162, 164, 165, 170, 170, 137, 147,
           148, 150, 151, 154, 155, 210, 132, 133, 135, 248, 67, 98, 110, 108,
           104),
    xw = c(158, 140, 149, 140, 158, 145, 158, 158, 159, 158, 158, 145, 170, 168,
           169, 169, 169, 168, 170, 179, 177, 179, 178, 179, 177, 179, 188, 187,
           188, 187, 188, 187, 188, 196, 196, 196, 195, 196, 196, 196, 207, 206,
           206, 207, 206, 206, 207, 217, 215, 216, 216, 216, 215, 217, 233, 229,
           229, 225, 222, 218, 218, 242, 238, 235, 233, 230, 226, 224, 242, 240,
           238, 235, 232, 229, 226, 235, 234, 230, 227, 224, 222, 218, 216, 213,
           210, 208, 206, 205, 189, 343, 186, 185, 182, 290, 144, 149, 154, 134,
           93),
    vw = c(128, 113, 120, 113, 128, 117, 129, 127, 128, 129, 127, 117, 136, 137,
           136, 136, 137, 136, 136, 144, 144, 143, 144, 144, 143, 144, 152, 153,
           151, 152, 153, 151, 152, 158, 159, 158, 158, 159, 158, 158, 168, 168,
           167, 168, 168, 167, 168, 175, 176, 174, 175, 176, 174, 175, 188, 186,
           183, 181, 180, 177, 175, 196, 193, 191, 189, 187, 184, 181, 197, 195,
           193, 190, 187, 185, 183, 190, 189, 186, 183, 182, 181, 177, 175, 172,
           170, 169, 167, 165, 153, 276, 152, 149, 147, 235, 116, 122, 124, 109,
           76))
  return(pedar_insole_areas)
}

pedar_insole_grids <- function() {
  pedar_insole_grid <- data.frame(
    V1 = c(0.643, 0.671, 0.696, 0.726, 0.778, 0.608, 0.636, 0.664, 0.696, 0.726,
           0.758, 0.786, 0.604, 0.632, 0.664, 0.696, 0.728, 0.758, 0.788, 0.6,
           0.632, 0.664, 0.692, 0.724, 0.758, 0.79, 0.596, 0.628, 0.66, 0.692,
           0.726, 0.758, 0.79, 0.592, 0.624, 0.656, 0.692, 0.726, 0.758, 0.79,
           0.588, 0.62, 0.656, 0.688, 0.722, 0.758, 0.79, 0.568, 0.604, 0.64,
           0.676, 0.714, 0.754, 0.79, 0.552, 0.592, 0.632, 0.672, 0.71, 0.75,
           0.79, 0.546, 0.584, 0.624, 0.664, 0.71, 0.75, 0.79, 0.536, 0.576,
           0.62, 0.66, 0.702, 0.742, 0.782, 0.536, 0.572, 0.612, 0.652, 0.692,
           0.728, 0.77, 0.536, 0.572, 0.608, 0.64, 0.676, 0.714, 0.744, 0.554,
           0.584, 0.616, 0.652, 0.684, 0.708, 0.572, 0.604, 0.636, 0.68, 0.157,
           0.129, 0.104, 0.074, 0.022, 0.192, 0.164, 0.136, 0.104, 0.074, 0.042,
           0.014, 0.196, 0.168, 0.136, 0.104, 0.072, 0.042, 0.012, 0.2, 0.168,
           0.136, 0.108, 0.076, 0.042, 0.01, 0.204, 0.172, 0.14, 0.108, 0.074,
           0.042, 0.01, 0.208, 0.176, 0.144, 0.108, 0.074, 0.042, 0.01, 0.212,
           0.18, 0.144, 0.112, 0.078, 0.042, 0.01, 0.232, 0.196, 0.16, 0.124,
           0.086, 0.046, 0.01, 0.248, 0.208, 0.168, 0.128, 0.09, 0.05, 0.01,
           0.254, 0.216, 0.176, 0.136, 0.09, 0.05, 0.01, 0.264, 0.224, 0.18,
           0.14, 0.098, 0.058, 0.018, 0.264, 0.228, 0.188, 0.148, 0.108, 0.072,
           0.03, 0.264, 0.228, 0.192, 0.16, 0.124, 0.086, 0.056, 0.246, 0.216,
           0.184, 0.148, 0.116, 0.092, 0.228, 0.196, 0.164, 0.12),
    V2 = c(0.154, 0.154, 0.154, 0.154, 0.154, 0.227, 0.227, 0.227, 0.227, 0.227,
           0.227, 0.227, 0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.36,
           0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.424, 0.424, 0.424, 0.424,
           0.424, 0.424, 0.424, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492,
           0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.612, 0.612, 0.612,
           0.612, 0.612, 0.612, 0.612, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664,
           0.664, 0.716, 0.716, 0.716, 0.716, 0.716, 0.716, 0.716, 0.77, 0.77,
           0.77, 0.77, 0.77, 0.77, 0.77, 0.822, 0.822, 0.822, 0.822, 0.822,
           0.822, 0.822, 0.878, 0.878, 0.878, 0.878, 0.878, 0.878, 0.878, 0.931,
           0.931, 0.931, 0.931, 0.931, 0.931, 0.99, 0.995, 0.99, 0.967, 0.154,
           0.154, 0.154, 0.154, 0.154, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227,
           0.227, 0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.36, 0.36,
           0.36, 0.36, 0.36, 0.36, 0.36, 0.424, 0.424, 0.424, 0.424, 0.424,
           0.424, 0.424, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.552,
           0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.612, 0.612, 0.612, 0.612,
           0.612, 0.612, 0.612, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664,
           0.716, 0.716, 0.716, 0.716, 0.716, 0.716, 0.716, 0.77, 0.77, 0.77,
           0.77, 0.77, 0.77, 0.77, 0.822, 0.822, 0.822, 0.822, 0.822, 0.822,
           0.822, 0.878, 0.878, 0.878, 0.878, 0.878, 0.878, 0.878, 0.931, 0.931,
           0.931, 0.931, 0.931, 0.931, 0.99, 0.995, 0.99, 0.967),
    V3 = c(0.592, 0.643, 0.671, 0.696, 0.726, 0.575, 0.608, 0.636, 0.664, 0.696,
           0.726, 0.758, 0.572, 0.604, 0.632, 0.664, 0.696, 0.728, 0.758, 0.568,
           0.6, 0.632, 0.664, 0.692, 0.724, 0.758, 0.564, 0.596, 0.628, 0.66,
           0.692, 0.726, 0.758, 0.56, 0.592, 0.624, 0.656, 0.692, 0.726, 0.758,
           0.556, 0.588, 0.62, 0.656, 0.688, 0.722, 0.758, 0.532, 0.568, 0.604,
           0.64, 0.676, 0.714, 0.754, 0.518, 0.552, 0.592, 0.632, 0.672, 0.71,
           0.75, 0.504, 0.546, 0.584, 0.624, 0.664, 0.71, 0.75, 0.495, 0.536,
           0.576, 0.62, 0.66, 0.702, 0.742, 0.495, 0.536, 0.572, 0.612, 0.652,
           0.692, 0.728, 0.503, 0.536, 0.572, 0.608, 0.64, 0.676, 0.714, 0.52,
           0.554, 0.584, 0.616, 0.652, 0.684, 0.544, 0.572, 0.604, 0.636, 0.208,
           0.157, 0.129, 0.104, 0.074, 0.225, 0.192, 0.164, 0.136, 0.104, 0.074,
           0.042, 0.228, 0.196, 0.168, 0.136, 0.104, 0.072, 0.042, 0.232, 0.2,
           0.168, 0.136, 0.108, 0.076, 0.042, 0.236, 0.204, 0.172, 0.14, 0.108,
           0.074, 0.042, 0.24, 0.208, 0.176, 0.144, 0.108, 0.074, 0.042, 0.244,
           0.212, 0.18, 0.144, 0.112, 0.078, 0.042, 0.268, 0.232, 0.196, 0.16,
           0.124, 0.086, 0.046, 0.282, 0.248, 0.208, 0.168, 0.128, 0.09, 0.05,
           0.296, 0.254, 0.216, 0.176, 0.136, 0.09, 0.05, 0.305, 0.264, 0.224,
           0.18, 0.14, 0.098, 0.058, 0.305, 0.264, 0.228, 0.188, 0.148, 0.108,
           0.072, 0.297, 0.264, 0.228, 0.192, 0.16, 0.124, 0.086, 0.28, 0.246,
           0.216, 0.184, 0.148, 0.116, 0.256, 0.228, 0.196, 0.164),
    V4 = c(0.154, 0.154, 0.154, 0.154, 0.154, 0.227, 0.227, 0.227, 0.227, 0.227,
           0.227, 0.227, 0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.36,
           0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.424, 0.424, 0.424, 0.424,
           0.424, 0.424, 0.424, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492,
           0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.612, 0.612, 0.612,
           0.612, 0.612, 0.612, 0.612, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664,
           0.664, 0.716, 0.716, 0.716, 0.716, 0.716, 0.716, 0.716, 0.77, 0.77,
           0.77, 0.77, 0.77, 0.77, 0.77, 0.822, 0.822, 0.822, 0.822, 0.822,
           0.822, 0.822, 0.878, 0.878, 0.878, 0.878, 0.878, 0.878, 0.878, 0.931,
           0.931, 0.931, 0.931, 0.931, 0.931, 0.967, 0.99, 0.995, 0.99, 0.154,
           0.154, 0.154, 0.154, 0.154, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227,
           0.227, 0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.36, 0.36,
           0.36, 0.36, 0.36, 0.36, 0.36, 0.424, 0.424, 0.424, 0.424, 0.424,
           0.424, 0.424, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.552,
           0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.612, 0.612, 0.612, 0.612,
           0.612, 0.612, 0.612, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664,
           0.716, 0.716, 0.716, 0.716, 0.716, 0.716, 0.716, 0.77, 0.77, 0.77,
           0.77, 0.77, 0.77, 0.77, 0.822, 0.822, 0.822, 0.822, 0.822, 0.822,
           0.822, 0.878, 0.878, 0.878, 0.878, 0.878, 0.878, 0.878, 0.931, 0.931,
           0.931, 0.931, 0.931, 0.931, 0.967, 0.99, 0.995, 0.99),
    V5 = c(0.614, 0.631, 0.664, 0.695, 0.73, 0.592, 0.616, 0.643, 0.671, 0.696,
           0.726, 0.75, 0.575, 0.608, 0.636, 0.664, 0.696, 0.726, 0.758, 0.572,
           0.604, 0.632, 0.664, 0.696, 0.728, 0.758, 0.568, 0.6, 0.632, 0.664,
           0.692, 0.724, 0.758, 0.564, 0.596, 0.628, 0.66, 0.692, 0.726, 0.758,
           0.56, 0.592, 0.624, 0.656, 0.692, 0.726, 0.758, 0.556, 0.588, 0.62,
           0.656, 0.688, 0.722, 0.758, 0.532, 0.568, 0.604, 0.64, 0.676, 0.714,
           0.754, 0.518, 0.552, 0.592, 0.632, 0.672, 0.71, 0.75, 0.504, 0.546,
           0.584, 0.624, 0.664, 0.71, 0.75, 0.495, 0.536, 0.576, 0.62, 0.66,
           0.702, 0.742, 0.495, 0.536, 0.572, 0.612, 0.652, 0.692, 0.728, 0.503,
           0.544, 0.584, 0.62, 0.66, 0.702, 0.52, 0.584, 0.616, 0.652, 0.186,
           0.169, 0.136, 0.105, 0.07, 0.208, 0.184, 0.157, 0.129, 0.104, 0.074,
           0.05, 0.225, 0.192, 0.164, 0.136, 0.104, 0.074, 0.042, 0.228, 0.196,
           0.168, 0.136, 0.104, 0.072, 0.042, 0.232, 0.2, 0.168, 0.136, 0.108,
           0.076, 0.042, 0.236, 0.204, 0.172, 0.14, 0.108, 0.074, 0.042, 0.24,
           0.208, 0.176, 0.144, 0.108, 0.074, 0.042, 0.244, 0.212, 0.18, 0.144,
           0.112, 0.078, 0.042, 0.268, 0.232, 0.196, 0.16, 0.124, 0.086, 0.046,
           0.282, 0.248, 0.208, 0.168, 0.128, 0.09, 0.05, 0.296, 0.254, 0.216,
           0.176, 0.136, 0.09, 0.05, 0.305, 0.264, 0.224, 0.18, 0.14, 0.098,
           0.058, 0.305, 0.264, 0.228, 0.188, 0.148, 0.108, 0.072, 0.297, 0.256,
           0.216, 0.18, 0.14, 0.098, 0.28, 0.216, 0.184, 0.148),
    V6 = c(0.11, 0.091, 0.079, 0.078, 0.09, 0.154, 0.154, 0.154, 0.154, 0.154,
           0.154, 0.154, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227, 0.295,
           0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.36, 0.36, 0.36, 0.36,
           0.36, 0.36, 0.36, 0.424, 0.424, 0.424, 0.424, 0.424, 0.424, 0.424,
           0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.552, 0.552, 0.552,
           0.552, 0.552, 0.552, 0.552, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612,
           0.612, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.716, 0.716,
           0.716, 0.716, 0.716, 0.716, 0.716, 0.77, 0.77, 0.77, 0.77, 0.77,
           0.77, 0.77, 0.822, 0.822, 0.822, 0.822, 0.822, 0.822, 0.822, 0.878,
           0.878, 0.878, 0.878, 0.878, 0.878, 0.931, 0.931, 0.931, 0.931, 0.11,
           0.091, 0.079, 0.078, 0.09, 0.154, 0.154, 0.154, 0.154, 0.154, 0.154,
           0.154, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227, 0.295, 0.295,
           0.295, 0.295, 0.295, 0.295, 0.295, 0.36, 0.36, 0.36, 0.36, 0.36,
           0.36, 0.36, 0.424, 0.424, 0.424, 0.424, 0.424, 0.424, 0.424, 0.492,
           0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.552, 0.552, 0.552, 0.552,
           0.552, 0.552, 0.552, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612,
           0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.716, 0.716, 0.716,
           0.716, 0.716, 0.716, 0.716, 0.77, 0.77, 0.77, 0.77, 0.77, 0.77, 0.77,
           0.822, 0.822, 0.822, 0.822, 0.822, 0.822, 0.822, 0.878, 0.878, 0.878,
           0.878, 0.878, 0.878, 0.931, 0.931, 0.931, 0.931),
    V7 = c(0.631, 0.664, 0.696, 0.73, 0.76, 0.616, 0.643, 0.671, 0.696, 0.726,
           0.75, 0.778, 0.608, 0.636, 0.664, 0.696, 0.726, 0.758, 0.786, 0.604,
           0.632, 0.664, 0.696, 0.728, 0.758, 0.788, 0.6, 0.632, 0.664, 0.692,
           0.724, 0.758, 0.79, 0.596, 0.628, 0.66, 0.692, 0.726, 0.758, 0.79,
           0.592, 0.624, 0.656, 0.692, 0.726, 0.758, 0.79, 0.588, 0.62, 0.656,
           0.688, 0.722, 0.758, 0.79, 0.568, 0.604, 0.64, 0.676, 0.714, 0.754,
           0.79, 0.552, 0.592, 0.632, 0.672, 0.71, 0.75, 0.79, 0.546, 0.584,
           0.624, 0.664, 0.71, 0.75, 0.79, 0.536, 0.576, 0.62, 0.66, 0.702,
           0.742, 0.782, 0.536, 0.572, 0.612, 0.652, 0.692, 0.728, 0.77, 0.544,
           0.584, 0.62, 0.66, 0.702, 0.744, 0.584, 0.616, 0.652, 0.708, 0.169,
           0.136, 0.104, 0.07, 0.04, 0.184, 0.157, 0.129, 0.104, 0.074, 0.05,
           0.022, 0.192, 0.164, 0.136, 0.104, 0.074, 0.042, 0.014, 0.196, 0.168,
           0.136, 0.104, 0.072, 0.042, 0.012, 0.2, 0.168, 0.136, 0.108, 0.076,
           0.042, 0.01, 0.204, 0.172, 0.14, 0.108, 0.074, 0.042, 0.01, 0.208,
           0.176, 0.144, 0.108, 0.074, 0.042, 0.01, 0.212, 0.18, 0.144, 0.112,
           0.078, 0.042, 0.01, 0.232, 0.196, 0.16, 0.124, 0.086, 0.046, 0.01,
           0.248, 0.208, 0.168, 0.128, 0.09, 0.05, 0.01, 0.254, 0.216, 0.176,
           0.136, 0.09, 0.05, 0.01, 0.264, 0.224, 0.18, 0.14, 0.098, 0.058,
           0.018, 0.264, 0.228, 0.188, 0.148, 0.108, 0.072, 0.03, 0.256, 0.216,
           0.18, 0.14, 0.098, 0.056, 0.216, 0.184, 0.148, 0.092),
    V8 = c(0.091, 0.079, 0.078, 0.09, 0.116, 0.154, 0.154, 0.154, 0.154, 0.154,
           0.154, 0.154, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227, 0.295,
           0.295, 0.295, 0.295, 0.295, 0.295, 0.295, 0.36, 0.36, 0.36, 0.36,
           0.36, 0.36, 0.36, 0.424, 0.424, 0.424, 0.424, 0.424, 0.424, 0.424,
           0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.552, 0.552, 0.552,
           0.552, 0.552, 0.552, 0.552, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612,
           0.612, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.716, 0.716,
           0.716, 0.716, 0.716, 0.716, 0.716, 0.77, 0.77, 0.77, 0.77, 0.77,
           0.77, 0.77, 0.822, 0.822, 0.822, 0.822, 0.822, 0.822, 0.822, 0.878,
           0.878, 0.878, 0.878, 0.878, 0.878, 0.931, 0.931, 0.931, 0.931, 0.091,
           0.079, 0.078, 0.09, 0.116, 0.154, 0.154, 0.154, 0.154, 0.154, 0.154,
           0.154, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227, 0.227, 0.295, 0.295,
           0.295, 0.295, 0.295, 0.295, 0.295, 0.36, 0.36, 0.36, 0.36, 0.36,
           0.36, 0.36, 0.424, 0.424, 0.424, 0.424, 0.424, 0.424, 0.424, 0.492,
           0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.552, 0.552, 0.552, 0.552,
           0.552, 0.552, 0.552, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612,
           0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.664, 0.716, 0.716, 0.716,
           0.716, 0.716, 0.716, 0.716, 0.77, 0.77, 0.77, 0.77, 0.77, 0.77, 0.77,
           0.822, 0.822, 0.822, 0.822, 0.822, 0.822, 0.822, 0.878, 0.878, 0.878,
           0.878, 0.878, 0.878, 0.931, 0.931, 0.931, 0.931))
  return(pedar_insole_grid)
}
