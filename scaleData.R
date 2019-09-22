scaleData <- function(data.use, do.center = T) {
    data.use <- t(scale(t(data.use), center = do.center, scale = TRUE))
    return(data.use)
}
