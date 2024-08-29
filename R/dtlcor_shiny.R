
#' @export
#'
dtl_shiny <- function(appname = "shiny") {

    appDir <- system.file(appname, package = "dtlcor")
    if (appDir == "") {
        stop("Could not find Shiny directory. Please try re-installing 'dtlcor'.",
             call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal");
}
