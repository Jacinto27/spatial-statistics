gd_to_sf <- function(gd, value_name = "z", crs = NA) {
  stopifnot(inherits(gd, "geodata"))
  df <- data.frame(
    x = gd$coords[, 1],
    y = gd$coords[, 2]
  )
  df[[value_name]] <- gd$data
  st_as_sf(df, coords = c("x", "y"), crs = crs)
}
