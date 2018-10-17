# waveclock.auto.R
#
# T.S.Price 2008
#
# run waveclock in another instance of R
# this can be useful when running clockwave over many datasets in a loop
#
waveclock.auto <-
function( ... )
{
  args <- list( ... )
  results <- NULL
  dir0 <- getwd()
  dir1 <- gsub( "[\\]", "/", tempdir() )
  dir1 <- gsub( "((?:[^/]*/)*).*$", "\\1", dir1, perl = TRUE )
  Rfile <- "temp-Rcode.Rdata"
  outfile <- "temp-output.Rdata"
  code <- character( 5 )
  code[ 1 ] <- "require( Rwave )\n"
  code[ 2 ] <- "require( waveclock )\n"
  code[ 3 ] <- paste( "results <- do.call( \"waveclock\", ", paste( deparse( args ), collapse = "" ) ," )\n", sep = "" )
  code[ 4 ] <- paste( "setwd( \"", dir1, "\" )\n", sep = "" )
  code[ 5 ] <- paste( "save( results, file = \"", outfile, "\" )\n", sep = "" )
  setwd( dir1 )
  cat( code, file = Rfile )
  results <- system( paste( R.home(), "/bin/R --vanilla --file=", Rfile, sep = "" ), intern = TRUE )
  load( outfile )
  unlink( outfile )
  setwd( dir0 )
  return( results )
}
