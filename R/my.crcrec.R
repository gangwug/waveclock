# my.crcrec.R
#
# T.S.Price 2008
#
# my.crcrec is an internal function used by clockwave to reconstruct
# the original signal from 'ridges' (modal frequencies) in the continuous
# wavelet decomposition.
#
# it is a slight alteration of the function crcrec in package Rwave-1.24-2
#
my.crcrec <-
function (
  siginput,
  inputwt,
  beemap,
  cfam,
  noct,
  nvoice,
  compr,
  minnbnodes = 2,
  w0 = 6.203608,
  epsilon = 0,
  fast = FALSE,
  para = 3,
  real = FALSE,
  plot = FALSE
)
{
  chain <- cfam$chain
  nbchain <- cfam$nbchain
  ordered <- cfam$ordered
  sigsize <- length(siginput)
  rec <- numeric(sigsize)
  plnb <- 0
  if ( plot != FALSE )
  {
    par( mfrow = c( 2, 1 ) )
    plot.ts( siginput, main = "Original signal" )
    image( cfam$ordered, main = "Chained Ridges" )
  }
  tmp <- matrix( 0, nbchain, length( siginput ) )
  totnbnodes <- 0
  idx <- numeric( nbchain )
  p <- 0
  for ( j in 1:nbchain )
  {
    phi.x.min <- 2 * 2^( chain[ j, 3 ] / nvoice )
    if ( chain[ j, 2 ] > ( para * phi.x.min ) )
    {
      cat( "Chain number", j )
      phi.x.max <- 2 * 2^( chain[ j, ( 2 + chain[ j, 2 ] ) ] / nvoice )
      x.min <- chain[ j, 1 ]
      x.max <- chain[ j, 1 ] + chain[ j, 2 ] - 1
      x.min <- x.min - round( para * phi.x.min )
      x.max <- x.max + round( para * phi.x.max )
      tmp2 <- regrec(
        siginput[
          chain[ j, 1 ]:( chain[ j, 1 ] + chain[ j, 2 ] - 1 ) ],
          inputwt[ chain[ j, 1 ]:( chain[ j, 1 ] + chain[j, 2] - 1),
        ],
        chain[ j, 3:( chain[ j, 2 ] + 2 ) ],
        compr,
        noct,
        nvoice,
        epsilon,
        w0 = w0,
        fast = fast,
        para = para,
        minnbnodes = minnbnodes,
        real = real
      )
      if ( is.list( tmp2 ) == TRUE )
      {
        totnbnodes <- totnbnodes + tmp2$nbnodes
        np <- length( tmp2$sol )
        start <- max( 1, x.min )
        end <- min( sigsize, x.min + np - 1 )
        start1 <- max( 1, 2 - x.min )
        end1 <- min( np, sigsize + 1 - x.min )
        end <- end1 - start1 + start
        rec[ start:end ] <- rec[ start:end ] + tmp2$sol[ start1:end1 ]
        tmp[ j, start:end ] <- Re( tmp2$sol[ start1:end1 ] )
      }
      plnb <- plnb + 1
      p <- p + 1
      idx[ p ] <- j
    }
  }
  if ( plot == 1 )
  {
    par( mfrow = c( 2, 1 ) )
    par( cex = 1.1 )
    plot.ts( siginput, main = "Original signal" )
    plot.ts( Re( rec ), main = "Reconstructed signal" )
  }
  else {
    if ( plot == 2 )
    {
      par( mfrow = c( plnb + 2, 1 ) )
      par( mar = c( 2, 0, 0, 0 ) )
      par( cex = 1.1 )
      par( err = -1 )
      plot.ts( siginput, main = "Original signal" )
      for ( j in 1:p ) plot.ts( tmp[ idx[ j ], ] )
      plot.ts( Re( rec ), main = "Reconstructed signal" )
    }
  }
  cat( "Total number of ridge samples used: ", totnbnodes, "\n" )
  par( mfrow = c( 1, 1 ) )
  list(
    rec = rec,
    ordered = ordered,
    chain = chain,
    comp = tmp,
    idx = idx
  )
}
