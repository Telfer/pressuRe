# Tests for the pressuRe package
# Install:  install.packages("pressuRe")
#           remotes::install_github("Telfer/pressuRe")
# Run with: testthat::test_file("test-pressuRe.R")
#
# Current pressuRe data structure (list with named fields):
#   $pressure_array   – 2D numeric matrix, rows = timepoints, cols = sensors
#   $pressure_system  – character string identifying the hardware system
#   $sens_size        – character string identifying the sensor type/size
#   $time             – numeric, time interval between measurements (seconds)
#   $masks            – list of masks
#   $events           – list of event markers
#   $sensor_polygons  – data frame with sensor corner coordinates
#   $max_matrix       – matrix of per-sensor maximum values
#
# NOTE: if any function name or argument has drifted from what is tested here,
# check ls("package:pressuRe") and update accordingly.

library(testthat)
library(pressuRe)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Number of fake sensors and timepoints used across tests
N_SENSORS   <- 48L   # columns of pressure_array
N_TIMEPOINTS <- 100L  # rows of pressure_array

# Build a minimal valid pressuRe list with uniform pressure
make_pressure_obj <- function(n_timepoints = N_TIMEPOINTS,
                              n_sensors    = N_SENSORS,
                              value        = 200,
                              time         = 0.01,   # 100 Hz
                              pressure_system = "emed",
                              sens_size    = 0.005 * 0.005) {

  pressure_array <- matrix(value,
                           nrow = n_timepoints,
                           ncol = n_sensors)

  # Minimal sensor_polygons: one row per sensor with an x/y centroid
  sensor_polygons <- data.frame(
    sensor_id = seq_len(n_sensors),
    x         = rep(seq_len(n_sensors / 8), each = 8),
    y         = rep(seq_len(8), times = n_sensors / 8)
  )

  # max_matrix: per-sensor max (for uniform data, equals value)
  max_matrix <- matrix(value, nrow = 1L, ncol = n_sensors)

  list(
    pressure_array  = pressure_array,
    pressure_system = pressure_system,
    sens_size       = sens_size,
    time            = time,
    masks           = list(),
    events          = list(),
    sensor_polygons = sensor_polygons,
    max_matrix      = max_matrix
  )
}

# Ramp: pressure rises from 0 to peak and back to 0 (sinusoidal loading curve)
make_ramp_obj <- function(n_timepoints = N_TIMEPOINTS,
                          n_sensors    = N_SENSORS,
                          peak         = 300,
                          time         = 0.01) {
  ramp <- outer(sin(pi * seq_len(n_timepoints) / n_timepoints),
                rep(peak, n_sensors))
  obj <- make_pressure_obj(n_timepoints, n_sensors, value = 0, time = time)
  obj$pressure_array <- ramp
  obj$max_matrix     <- matrix(apply(ramp, 2, max), nrow = 1L)
  obj
}

make_zero_obj <- function(...) make_pressure_obj(..., value = 0)

# ---------------------------------------------------------------------------
# 1.  Package loads and key symbols are exported
# ---------------------------------------------------------------------------

test_that("pressuRe package loads without errors", {
  expect_true("package:pressuRe" %in% search())
})

test_that("core functions are exported", {
  exported <- ls("package:pressuRe")
  expected <- c(
    "load_emed", "load_pedar",
    "pressure_interp", "auto_detect_side",
    "cop", "cpei",
    "whole_pressure_curve", "mask_analysis",
    "create_mask_auto", "create_mask_manual"
  )
  for (fn in expected) {
    expect_true(fn %in% exported,
                info = paste("Expected export missing:", fn))
  }
})

# ---------------------------------------------------------------------------
# 2.  Data structure integrity helpers
# ---------------------------------------------------------------------------

test_that("make_pressure_obj produces the expected list structure", {
  obj <- make_pressure_obj()
  expect_true(is.list(obj))
  for (field in c("pressure_array", "pressure_system", "sens_size",
                  "time", "masks", "events", "sensor_polygons", "max_matrix")) {
    expect_true(field %in% names(obj), info = paste("Field missing:", field))
  }
})

test_that("pressure_array is a 2D matrix with correct dimensions", {
  obj <- make_pressure_obj(n_timepoints = 80, n_sensors = 32)
  expect_true(is.matrix(obj$pressure_array))
  expect_equal(nrow(obj$pressure_array), 80)
  expect_equal(ncol(obj$pressure_array), 32)
})

test_that("max_matrix columns match column-wise maxima of pressure_array", {
  obj <- make_ramp_obj()
  expected_max <- apply(obj$pressure_array, 2, max)
  expect_equal(as.numeric(obj$max_matrix), expected_max, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 3.  pressure_interp – temporal resampling
# ---------------------------------------------------------------------------

test_that("pressure_interp returns object with resampled pressure_array row count", {
  obj <- make_pressure_obj(n_timepoints = 60)
  out <- pressure_interp(pressure_data = obj, interp_to = 101)
  expect_equal(nrow(out$pressure_array), 101)
})

test_that("pressure_interp preserves the number of sensors (columns)", {
  obj <- make_pressure_obj(n_timepoints = 60, n_sensors = N_SENSORS)
  out <- pressure_interp(pressure_data = obj, interp_to = 101)
  expect_equal(ncol(out$pressure_array), N_SENSORS)
})

test_that("pressure_interp preserves a constant signal exactly", {
  obj <- make_pressure_obj(value = 150, n_timepoints = 60)
  out <- pressure_interp(pressure_data = obj, interp_to = 101)
  expect_true(all(abs(out$pressure_array - 150) < 1e-6))
})

test_that("pressure_interp preserves non-array metadata fields", {
  obj <- make_pressure_obj(n_timepoints = 60)
  out <- pressure_interp(pressure_data = obj, interp_to = 101)
  expect_equal(out$pressure_system, obj$pressure_system)
  expect_equal(out$sens_size,       obj$sens_size)
})

test_that("pressure_interp errors on non-positive n_frames", {
  obj <- make_pressure_obj()
  expect_error(pressure_interp(pressure_data = obj, interp_to = 0))
  expect_error(pressure_interp(pressure_data = obj, interp_to = -5))
})

test_that("pressure_interp errors when pressure_array field is missing", {
  bad <- make_pressure_obj()
  bad$pressure_array <- NULL
  expect_error(pressure_interp(pressure_data = bad, interp_to = 101))
})

# ---------------------------------------------------------------------------
# 4.  whole_pressure_curve – whole-foot summary curves
# ---------------------------------------------------------------------------

test_that("whole_pressure_curve (peak_pressure) returns a numeric vector", {
  obj <- make_pressure_obj()
  out <- whole_pressure_curve(pressure_data = obj, variable = "peak_pressure",
                              plot = FALSE)
  expect_true(is.numeric(out))
})

test_that("whole_pressure_curve vector length equals number of timepoints", {
  obj <- make_pressure_obj(n_timepoints = 80)
  out <- whole_pressure_curve(pressure_data = obj, variable = "peak_pressure",
                              plot = FALSE)
  expect_equal(length(out), 80)
})

test_that("whole_pressure_curve (contact_area) returns a numeric vector", {
  obj <- make_pressure_obj()
  out <- whole_pressure_curve(pressure_data = obj, variable = "area",
                              plot = FALSE)
  expect_true(is.numeric(out))
})

test_that("whole_pressure_curve (force) returns a numeric vector", {
  obj <- make_pressure_obj()
  out <- whole_pressure_curve(pressure_data = obj, variable = "force",
                              plot = FALSE)
  expect_true(is.numeric(out))
})

test_that("whole_pressure_curve peak_pressure is 0 everywhere for zero data", {
  obj <- make_zero_obj()
  out <- whole_pressure_curve(pressure_data = obj, variable = "peak_pressure",
                              plot = FALSE)
  expect_true(all(out == 0))
})

test_that("whole_pressure_curve contact_area is 0 everywhere for zero data", {
  obj <- make_zero_obj()
  out <- whole_pressure_curve(pressure_data = obj, variable = "area",
                              plot = FALSE)
  expect_true(all(out == 0))
})

test_that("whole_pressure_curve peak_pressure matches known constant value", {
  obj <- make_pressure_obj(value = 350)
  out <- whole_pressure_curve(pressure_data = obj, variable = "peak_pressure",
                              plot = FALSE)
  expect_true(all(abs(out - 350) < 1e-6))
})

test_that("whole_pressure_curve peak_pressure is higher for higher pressure", {
  low  <- make_pressure_obj(value = 100)
  high <- make_pressure_obj(value = 400)
  out_low  <- whole_pressure_curve(low,  variable = "peak_pressure", plot = FALSE)
  out_high <- whole_pressure_curve(high, variable = "peak_pressure", plot = FALSE)
  expect_true(all(out_high > out_low))
})

test_that("whole_pressure_curve errors on unknown metric", {
  obj <- make_pressure_obj()
  expect_error(
    whole_pressure_curve(pressure_data = obj, variable = "not_a_metric",
                         plot = FALSE)
  )
})

# ---------------------------------------------------------------------------
# 5.  cop – center of pressure
# ---------------------------------------------------------------------------

test_that("cop returns a two-column data frame or matrix", {
  emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
  obj <- load_emed(emed_data)
  result <- cop(pressure_data = obj)
  expect_true(is.data.frame(result) || is.matrix(result))
  expect_equal(ncol(result), 2)
})

test_that("cop has one row per timepoint", {
  emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
  obj <- load_emed(emed_data)
  result <- cop(pressure_data = obj)
  expect_equal(nrow(result), 106)
})

test_that("cop is stable (same result on repeated calls)", {
  emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
  obj <- load_emed(emed_data)
  expect_equal(cop(obj), cop(obj))
})

test_that("cop errors on non-list input", {
  expect_error(cop(pressure_data = matrix(1, 10, 10)))
})

# ---------------------------------------------------------------------------
# 6.  cpei – center of pressure excursion index
# ---------------------------------------------------------------------------

# test_that("cpei returns a single non-negative numeric scalar", {
#   emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#   obj <- load_emed(emed_data)
#   result <- cpei(pressure_data = obj)
#   expect_length(result, 1)
#   expect_true(is.numeric(result))
#   expect_gte(result, 0)
# })
#
# test_that("cpei is reproducible", {
#   emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#   obj <- load_emed(emed_data)
#   expect_equal(cpei(obj), cpei(obj))
# })

# ---------------------------------------------------------------------------
# 7.  mask_analysis – regional pressure metrics
# ---------------------------------------------------------------------------

# Build a simple logical mask: TRUE for first half of sensors
# make_sensor_mask <- function(n_sensors = N_SENSORS) {
#   mask <- rep(FALSE, n_sensors)
#   mask[seq_len(n_sensors %/% 2)] <- TRUE
#   mask
# }

# test_that("mask_analysis returns a data frame or named list", {
#   obj  <- make_pressure_obj()
#   mask <- make_sensor_mask()
#   out  <- mask_analysis(pressure_data = obj, mask = mask)
#   expect_true(is.data.frame(out) || is.list(out))
# })

# test_that("mask_analysis output contains expected metric names", {
#   obj  <- make_pressure_obj()
#   mask <- make_sensor_mask()
#   out  <- mask_analysis(pressure_data = obj, mask = mask)
#   for (metric in c("peak_pressure", "contact_area", "pti")) {
#     expect_true(metric %in% names(out),
#                 info = paste("Expected metric missing:", metric))
#   }
# })

# test_that("mask_analysis peak_pressure is 0 for all-zero data", {
#   obj  <- make_zero_obj()
#   mask <- make_sensor_mask()
#   pp   <- mask_analysis(pressure_data = obj, mask = mask)$peak_pressure
#   expect_equal(pp, 0)
# })

# test_that("mask_analysis contact_area is 0 for all-zero data", {
#   obj  <- make_zero_obj()
#   mask <- make_sensor_mask()
#   ca   <- mask_analysis(pressure_data = obj, mask = mask)$contact_area
#   expect_equal(ca, 0)
# })

# test_that("mask_analysis peak_pressure matches known uniform value", {
#   val  <- 400
#   obj  <- make_pressure_obj(value = val)
#   mask <- make_sensor_mask()
#   pp   <- mask_analysis(pressure_data = obj, mask = mask)$peak_pressure
#   expect_equal(pp, val)
# })

# test_that("mask_analysis pti increases with higher uniform pressure", {
#   low  <- make_pressure_obj(value = 100)
#   high <- make_pressure_obj(value = 500)
#   mask <- make_sensor_mask()
#   expect_lt(mask_analysis(low,  mask)$pti,
#             mask_analysis(high, mask)$pti)
# })

# test_that("mask_analysis sub-region contact_area <= whole-foot contact_area", {
#   obj       <- make_pressure_obj()
#   full_mask <- rep(TRUE, N_SENSORS)
#   sub_mask  <- make_sensor_mask()
#   ca_full   <- mask_analysis(obj, full_mask)$contact_area
#   ca_sub    <- mask_analysis(obj, sub_mask)$contact_area
#   expect_gte(ca_full, ca_sub)
# })


# ---------------------------------------------------------------------------
# 8.  auto_detect_side – left/right foot classification
# ---------------------------------------------------------------------------

test_that("auto_detect_side returns 'RIGHT'", {
  emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
  obj <- load_emed(emed_data)
  result <- auto_detect_side(pressure_data = obj)
  expect_equal(result, "RIGHT")
})

# ---------------------------------------------------------------------------
# 9.  Input validation (shared)
# ---------------------------------------------------------------------------

test_that("whole_pressure_curve errors when pressure_array is missing", {
  bad <- make_pressure_obj()
  bad$pressure_array <- NULL
  expect_error(
    whole_pressure_curve(pressure_data = bad, metric = "peak_pressure",
                         plot = FALSE)
  )
})

test_that("cop errors when pressure_array is missing", {
  bad <- make_pressure_obj()
  bad$pressure_array <- NULL
  expect_error(cop(pressure_data = bad))
})

test_that("pressure_interp errors on non-list input", {
  expect_error(
    pressure_interp(pressure_data = matrix(1, 10, 10), n_frames = 101)
  )
})

# ---------------------------------------------------------------------------
# 10.  Integration / pipeline tests
# ---------------------------------------------------------------------------

# test_that("interp → mask_analysis pipeline is deterministic", {
#   obj  <- make_ramp_obj(n_timepoints = 60)
#   out  <- pressure_interp(obj, interp_to = 101)
#   mask <- make_sensor_mask()
#   res1 <- mask_analysis(out, mask)
#   res2 <- mask_analysis(out, mask)
#   expect_equal(res1$peak_pressure, res2$peak_pressure)
#   expect_equal(res1$pti,           res2$pti)
# })

# test_that("max of whole_pressure_curve equals mask_analysis peak for full mask", {
#   obj       <- make_ramp_obj()
#   pp_curve  <- whole_pressure_curve(obj, metric = "peak_pressure", plot = FALSE)
#   full_mask <- rep(TRUE, N_SENSORS)
#   pp_mask   <- mask_analysis(obj, full_mask)$peak_pressure
#   expect_equal(max(pp_curve), pp_mask, tolerance = 1e-6)
# })

# test_that("pti is proportional to pressure magnitude for uniform data", {
#   dat_1x <- make_pressure_obj(value = 100)
#   dat_3x <- make_pressure_obj(value = 300)
#   mask   <- make_sensor_mask()
#   pti_1x <- mask_analysis(dat_1x, mask)$pti
#   pti_3x <- mask_analysis(dat_3x, mask)$pti
#   expect_equal(pti_3x / pti_1x, 3, tolerance = 0.01)
# })

# test_that("sub-region pti <= whole-foot pti for the same data", {
#   obj       <- make_ramp_obj()
#   full_mask <- rep(TRUE, N_SENSORS)
#   sub_mask  <- make_sensor_mask()
#   pti_full  <- mask_analysis(obj, full_mask)$pti
#   pti_sub   <- mask_analysis(obj, sub_mask)$pti
#   expect_gte(pti_full, pti_sub)
# })

test_that("pressure_interp followed by pressure_interp back is near-identical for smooth data", {
  obj    <- make_ramp_obj(n_timepoints = 101)
  down   <- pressure_interp(obj, interp_to = 51)
  up     <- pressure_interp(down, interp_to = 101)
  # Round-tripping a smooth signal should be close (not exact due to interpolation)
  expect_lt(max(abs(up$pressure_array - obj$pressure_array)), 15)
})
