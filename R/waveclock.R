# waveclock.R
#
# T.S.Price 2008
#
# continuous wavelet transform (cwt) using Morlet wavelet
# followed by ridge reconstruction using "crazy climbers" simulated annealing algorithm
# (Carmona, Hwang & Torresani: Practical Time-Frequency Series Analysis, Academic Press, 1998)
#
# function clockwave takes the following arguments:
#   x              time series data (class = "ts")
#   period         range defining lower and upper period limits
#   time.limits    time range for truncation of series
#   extend         ameliorates edge effects by reflecting data series at time limits (extend="reflect")
#                    or repeating time series (extend="repeat")
#   noctave        number of octaves ( default = floor( log2( length( x ) ) ) )
#   nvoice         number of voices per octave ( default = 96 )
#   mask.coi       mask "cone of influence" from output
#   crc.args       list of arguments for crazy climbers algorithm. Defaults:
#     seed = 0       seed for random number generator
#     nbclimb = 50   number of crazy climbers
#   cfamily.args   list of arguments for chaining algorithm. Defaults:
#     ptile = 0.005  relative threshold for the ridges
#     bstep = 5      maximal length for a gap in a ridge
#     nbchain = 200  maximal number of chains produced by the function
#   crcrec.args   list of arguments for chaining algorithm. Defaults:
#     compr = 3      compression rate for sampling the ridges (set to 1 for a slower, but more precise analysis)
#     epsilon = 0    constant in front of the smoothness term in penalty function
#     para = 3       scale parameter for extrapolating the ridges
#   xlab           x axis label for plot of cwt scalogram
#   ylab           y axis label for plot of cwt scalogram
#   png            name of png filename for plot output
#   color.palette  color palette used in scalogram plot
#   mode.col       color of line marking modal frequency
#   mode.lty       type of line marking modal frequency
#   mode.lwd       width of line marking modal frequency
#   ...            additional parameters passed to filled.contour plot function
#
# function clockwave gives the following output:
#   cwt          continuous wavelet transform
#   crc          crazy climber output
#   cfamily      output of procedure to find ridges (chains) in crazy climber output
#   modes        modes of signal corresponding to reconstructed ridges:
#                  index     which ridge, corresponding to row of cfamily$chain
#                  median    median of modal signal, excluding cone of influence
#                  mean      mean of modal signal, excluding cone of influence
#   rec          reconstruction of signal from ridges
#   amp          estimated instantaneous amplitudes of modes
#   per          estimated instantaneous periods of modes
#   phase        estimated instantaneous phase of modes
#   mask         mask indicating boundaries of "cone of influence" for modes
#
#
#
waveclock <-
function(
  x,
  period = c( 6, 48 ),
  time.limits = NULL,
  extend = "reflect",
  noctave = NULL,
  nvoice = 96,
  mask.coi = TRUE,
  crc.args = list(
    seed = 0,
    nbclimb = 50
  ),
  cfamily.args = list(
    ptile = 0.005,
    bstep = 5,
    nbchain = 400
  ),
  crcrec.args = list(
    compr = 3,
    epsilon = 0,
    para = 3,
    plot = FALSE
  ),
  xlab = 'Time (h)',
  ylab = 'Period (h)',
  png = NULL,
  color.palette = heat.colors,
  mode.col = "green",
  mode.lty = "solid",
  mode.lwd = 2,
  ...
)
{
  extend <- pmatch( extend, c( "reflect", "repeat" ) )
  do.reflect <-
  function( x )
  {
    dt <- deltat( x )
    rt <- range( time( x ) )
    x <- as.matrix( x )
    d <- dim( x )[ 1 ]
    ##z <- 2^ceiling( log2( d ) ) * 2
    ##zeros <- ts( matrix( 0, nr = z - d, nc = dim( x )[ 2 ] ), deltat = dt )
    x <- rbind( x[ d:2, , drop = FALSE ], x, x[ ( d - 1 ):1, , drop = FALSE ] )
    return( ts( x, deltat = dt, start = rt[ 1 ] - rt[ 2 ] ) )
  }
  do.repeat <-
  function( x )
  {
    dt <- deltat( x )
    rt <- range( time( x ) )
    x <- as.matrix( x )
    d <- dim( x )[ 1 ]
    x <- rbind( x, x, x )
    return( ts( x, deltat = dt, start = rt[ 1 ] - rt[ 2 ] - 1 ) )
  }
  orig <- x
  if ( class( x ) != "ts" )
  {
    x <- as.ts( x )
  }
  dt <- deltat( x )
  if ( is.null( time.limits ) )
  {
    time.limits <- range( time( x ) )
  }
  if ( length( time.limits ) == 1 )
  {
    time.limits <- c( time.limits, max( time( x ) ) )
  }
  if ( length( time.limits ) == 2 && is.na( time.limits[ 1 ] ) )
  {
    time.limits[ 1 ] <- min( time( x ) )
  }
  if ( length( time.limits ) == 2 && is.na( time.limits[ 2 ] ) )
  {
    time.limits[ 2 ] <- max( time( x ) )
  }

  # first time point = time.limits[ 1 ]
  x <- ts( as.matrix( x[ time( x ) >= time.limits[ 1 ] & time( x ) <= time.limits[ 2 ] ] ), deltat = dt, start = 0 )

  N    <- length( x )
  T    <- dt * N
  s0   <- dt * 2

  # do this before extending the series rather than after
  if ( is.null( noctave ) )
  {
    noctave <- floor( log2( length( x ) ) )
  }
  noct <- noctave

  # reflect x at ends of time series
  if( length( extend ) && extend == 1 )
  {
    x <- do.reflect( x )
  }
  # or repeat series
  else if( length( extend ) && extend == 2 )
  {
    x <- do.repeat( x )
  }

  # continuous wavelet transform
  t1 <- do.call( "cwt", list( input = x, noctave = noct, nvoice = nvoice, w0 = 6.203608, plot = FALSE ) )

  # ridge reconstruction
  if ( !is.null( crc.args$seed ) )
  {
    set.seed( crc.args$seed )
  }
  crc.args$tfrep <- t1
  c1 <- do.call( "my.crc", crc.args )
  cfamily.args$ccridge <- c1
  f1 <- do.call( "cfamily", cfamily.args )
  crcrec.args$siginput <- x
  crcrec.args$inputwt <- t1
  crcrec.args$beemap <- c1
  crcrec.args$cfam <- f1
  crcrec.args$noct <- noct
  crcrec.args$nvoice <- nvoice
  crcrec.args$w0 <- 6.203608
  crcrec.args$plot <- FALSE
  r1 <- do.call( "my.crcrec", crcrec.args )

  # plot transforms
  i    <- seq( 0, noct * nvoice ) / nvoice
  y    <- i * nvoice + 1
  p    <- s0 * 2^i
  q1   <- seq( 1, noct * nvoice + 1 ) * ( ( p > period[ 1 ] ) & ( p <= period[ 2 ] ) )
  s1   <- max( min( q1[ q1 > 0 ] ) - 1, 1 )
  s2   <- max( q1[ q1 > 0 ] )
  p1   <- p[ s1 ]
  p2   <- p[ s2 ]
  noct <- log2( p2 ) - log2( p1 )
  o    <- 0.5 / nvoice # offset on plot
  lt   <- unique( ( time( x )[ time( x ) >= 0 & time( x ) <= diff( time.limits ) ] + time.limits[ 1 ] ) %/% 12 ) * 12
  lp   <- c( 48, 24, 12, 6, 3, 1.5, 0.75 )
  pp   <- ( log2( lp ) - log2( p1 ) + 0 ) / ( noct + 1 / nvoice )
  bt   <- 1 / diff( time.limits )
  at   <- 0
  Keps <- 1000 * .Machine$double.eps

  # plot
  if ( !is.null( png ) )
  {
    png( paste( png, "png", sep = "." ) )
  }
  filled.contour(
    Mod( t1 )[ time( x ) >= 0 & time( x ) <= diff( time.limits ), y[ s1:s2 ] ],
    xlab = xlab,
    ylab = ylab,
    plot.axes =
    {
      axis( 1, at + bt * ( lt - time.limits[ 1 ] ), lab = lt )
      axis( 2, pp, lab = lp )

      # circadian period
      scale.24h <- ( log2( 24 ) - log2( p1 ) + 0 ) / ( noct + 1 / nvoice )
      abline( h = scale.24h, lty = 2 )

      # cone of influence
      s <- 2 ^ seq( 0, noct, 0.5 / nvoice ) * p1                            # scales ( = period, since fourier factor = 1 )
      if ( mask.coi )
      {
        cy1 <- ( log2( s ) - log2( p1 ) + o ) / ( noct + 1 / nvoice )
        cx1 <- sqrt( 2 ) * s / T
        sel1 <- ( cx1 < 1 / 2 )
        cy2 <- cy1
        cx2 <- 1 - cx1
        sel2 <- ( cx2 > 1 / 2 )
        ox <- sqrt( 2 ) * 2^( log2( p1 ) - o ) / T
        cx <- c( cx1[ sel1 ], 0.5, cx2[ sel2 ] )
        cy <- c( cy1[ sel1 ], ( log2( T ) - log2( p1 ) - 1.5 + o ) / ( noct + 1 / nvoice ), cy2[ sel2 ] )
        ord <- order( cx )
        cx <- cx[ ord ]
        cy <- cy[ ord ]
        lines( cx, cy )
        polygon( c( 0, 0, ox, cx, 1 - ox, 1, 1 ), c( 1, 0, 0, cy, 0, 0, 1 ), density = 5, angle = -45 )
        polygon( c( 0, 0, ox, cx, 1 - ox, 1, 1 ), c( 1, 0, 0, cy, 0, 0, 1 ), density = 5, angle =  45 )
      }

      # draw ridges
      chains <- r1$idx[ r1$idx > 0 ]
      rec <- matrix( NA, nc = max( c( 0, chains ) ), nr = N )
      amp <- matrix( NA, nc = max( c( 0, chains ) ), nr = N )
      per <- matrix( NA, nc = max( c( 0, chains ) ), nr = N )
      phase <- matrix( NA, nc = max( c( 0, chains ) ), nr = N )
      mask <- matrix( NA, nc = max( c( 0, chains ) ), nr = N )
      modes <- matrix( NA, nr = max( c( 0, chains + 1 ) ), nc = 8 )
      colnames( modes ) <- c(
        "index",
        "median (midpoint)",
        "median (lower limit)",
        "median (upper limit)",
        "mean (midpoint)",
        "mean (lower limit)",
        "mean (upper limit)",
        "var"
      )
      tsel <- time( x ) >= -Keps & time( x ) <= diff( time.limits ) + Keps
      d1 <- sum( tsel )
      d2 <- s2 - s1 + 1
      t0 <- which( abs( time( x ) ) < Keps )
      mm <- numeric( 0 )
      for( j in chains )
      {
        m  <- numeric( 0 )
        oy <- NA
        x1 <- seq( f1$chain[ j, 1 ], f1$chain[ j, 1 ] + f1$chain[ j, 2 ] - 1 )
        y1 <- f1$chain[ j, 3:( f1$chain[ j, 2 ] + 2 ) ]
        if ( min( y1 ) >= s1 && max( y1 ) <= s2 && max( x1 ) >= min( which( tsel ) ) && min( x1 ) <= max( which( tsel ) ) )
        {
          rec[ , j ] <- r1$comp[ j, tsel ]
          for( i in 1:length( x1 ) )
          {
            x2 <- ( x1[ i ] - t0 ) / d1
            y2 <- ( y1[ i ] - s1 ) / d2 + o / ( noct + 1 / nvoice )
            if ( is.na( oy ) )
            {
              oy <- y2
            }
            coi <- ( x2 >= sqrt( 2 ) * 2 ^ ( ( y1[ i ] - s1 ) / nvoice ) * p1 / T ) && ( ( 1 - x2 ) >= sqrt( 2 ) * 2 ^ ( ( y1[ i ] - s1 ) / nvoice ) * p1 / T )
            lines(
              x = c( x2, x2, x2 + 1 / d1 ),
              y = c( oy, y2, y2 ),
              col = mode.col,
              lty = mode.lty,
              lwd = mode.lwd
            )
            oy <- y2
            if ( ( !mask.coi && tsel[ x1[ i ] ] ) || ( !is.na( coi ) && coi ) )
            {
              m <- c( m, y1[ i ] )
              mask[ x1[ i ] - t0, j ] <- y1[ i ]
              per[ x1[ i ] - t0, j ] <- p[ y1[ i ] ]
              amp[ x1[ i ] - t0, j ] <- Mod( t1[ x1[ i ], y1[ i ] ] )
              phase[ x1[ i ] - t0, j ] <- Arg( t1[ x1[ i ], y1[ i ] ] )
            }
          }
        }
        print( j )
        if ( f1$chain[ j, 1 ] - t0 > 1 )
        {
          y1 <- f1$chain[ j, 3 ]
          for( x1 in 1:min( f1$chain[ j, 1 ] - t0, N ) )
          {
            x2 <- ( x1 - 1 ) / d1
            coi <- ( x2 >= sqrt( 2 ) * 2 ^ ( ( y1 - s1 ) / nvoice ) * p1 / T ) && ( ( 1 - x2 ) >= sqrt( 2 ) * 2 ^ ( ( y1 - s1 ) / nvoice ) * p1 / T )
            if( !mask.coi || ( is.na( mask[ x1, j ] ) && !is.na( coi ) && coi ) )
            {
              mask[ x1, j ] <- 0
              per[ x1, j ] <- 0
              amp[ x1, j ] <- 0
              phase[ x1, j ] <- 0
            }
          }
        }
        if ( ( f1$chain[ j, 1 ] + f1$chain[ j, 2 ] - t0 - 1 ) < N )
        {
          y1 <- f1$chain[ j, f1$chain[ j, 2 ] + 2 ]
          for( x1 in max( f1$chain[ j, 1 ] + f1$chain[ j, 2 ] - t0, 1 ):N )
          {
            x2 <- ( x1 - 1 ) / d1
            coi <- ( x2 >= sqrt( 2 ) * 2 ^ ( ( y1 - s1 ) / nvoice ) * p1 / T ) && ( ( 1 - x2 ) >= sqrt( 2 ) * 2 ^ ( ( y1 - s1 ) / nvoice ) * p1 / T )
            if( !mask.coi || ( is.na( mask[ x1, j ] ) && !is.na( coi ) && coi ) )
            {
              mask[ x1, j ] <- 0
              per[ x1, j ] <- 0
              amp[ x1, j ] <- 0
              phase[ x1, j ] <- 0
            }
          }
        }
        if ( length( m ) > 0 )
        {
          modes[ j, ] <- c(
            j,
            2 ^ ( ( median( m ) - 1 ) / nvoice ) * s0,
            2 ^ ( ( median( m ) - 1.5 ) / nvoice ) * s0,
            2 ^ ( ( median( m ) - 0.5 ) / nvoice ) * s0,
            2 ^ ( ( mean( m ) - 1 ) / nvoice ) * s0,
            2 ^ ( ( mean( m ) - 1.5 ) / nvoice ) * s0,
            2 ^ ( ( mean( m ) - 0.5 ) / nvoice ) * s0,
            var( rec[ , j ] * ( mask[ , j ] > 0 ), na.rm = TRUE )
          )
          mm <- c( mm, m )
        }
      }
      if ( length( mm ) > 0 )
      {
        v <- rowSums( rec * ( mask > 0 ), na.rm = TRUE )
        v[ is.na( rowMeans( rec * ( mask > 0 ), na.rm = TRUE ) ) ] <- NA
        modes[ j + 1, ] <- c(
          0,
          2 ^ ( ( median( mm ) - 1 ) / nvoice ) * s0,
          2 ^ ( ( median( mm ) - 1.5 ) / nvoice ) * s0,
          2 ^ ( ( median( mm ) - 0.5 ) / nvoice ) * s0,
          2 ^ ( ( mean( mm ) - 1 ) / nvoice ) * s0,
          2 ^ ( ( mean( mm ) - 1.5 ) / nvoice ) * s0,
          2 ^ ( ( mean( mm ) - 0.5 ) / nvoice ) * s0,
          var( v, na.rm = TRUE )
        )
      }
      if ( length( mm ) > 0 )
      {
        modes <- as.data.frame( modes )
        if ( sum( !is.na( modes[ , 1 ] ) ) > 0 )
        {
          modes <- modes[ !is.na( modes[ , 1 ] ), , drop = FALSE ]
          modes[ dim( modes )[ 1 ], 1 ] <- "Total"
        }
        else
        {
          modes <- NULL
        }
      }
      else
      {
        modes <- NULL
      }
      if ( !is.null( modes ) )
      {
        sel <- as.numeric( modes[ , 1 ][ seq( dim( modes )[ 1 ] - 1 ) ] )
        rec <- ts( rec[ , sel, drop = FALSE ], deltat = dt, start = time.limits[ 1 ] )
        per <- ts( per[ , sel, drop = FALSE ], deltat = dt, start = time.limits[ 1 ] )
        amp <- ts( amp[ , sel, drop = FALSE ], deltat = dt, start = time.limits[ 1 ] )
        phase <- ts( phase[ , sel, drop = FALSE ], deltat = dt, start = time.limits[ 1 ] )
        mask <- ts( mask[ , sel, drop = FALSE ], deltat = dt, start = time.limits[ 1 ] )
      }
      else
      {
        rec <- NULL
        per <- NULL
        amp <- NULL
        phase <- NULL
        mask <- NULL
      }
    },
    color.palette = color.palette,
    ...
  )
  if ( !is.null( png ) )
  {
    dev.off()
  }
  return(
    invisible(
      list(
        original.signal = orig,
        modified.signal = x,
        cwt = t1,
        crc = c1,
        cfamily = f1,
        crcrec = r1,
        modes = modes,
        rec = rec,
        per = per,
        amp = amp,
        phase = phase,
        mask = mask
      )
    )
  )
}
