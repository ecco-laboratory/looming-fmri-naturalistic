# version of write_csv that returns the file path as expected when called by targets
write_csv_target <- function (x, file, ...) {
  write_csv(x = x, file = file, ...)
  return (file)
}