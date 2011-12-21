zic <- function( formula, data, bbar, dbar, ebar, fbar, n.burnin, n.mcmc, n.thin )
{
  # data matrices and dimensions
  mdl <- model.frame( formula, data )
  y <- model.response( mdl )
  X <- model.matrix( formula, mdl )
  n <- dim(X)[1]
  k <- dim(X)[2]

  # initialize memories
  n.save = n.mcmc / n.thin
  betamem <- matrix( 0.0, n.save, k )
  deltamem <- matrix( 0.0, n.save, k )
  sigma2mem <- matrix( 0.0, n.save, 1 )
  gammamem <- matrix( 0.0, n.save, k ) 
  kappamem <- matrix( 0.0, n.save, k )
    
  # set up prior
  m_bbar <- matrix( 0.0, k, 2 )
  m_bbar[,2] <- bbar
  m_dbar <- matrix( 0.0, k, 2 )
  m_dbar[,2] <- dbar
  pbar <- vector( "numeric", k )
  qbar <- vector( "numeric", k )

  # call C++ code
  draws <- .C( "zic",
               as.integer( y ), as.double( X ),
               as.integer( n ), as.integer( k ),
               as.integer( 0 ),
               as.integer( n.burnin ), as.integer( n.mcmc ), as.integer( n.thin ),
               as.double( m_bbar ), as.double( m_dbar ), as.double( ebar ), as.double( fbar ), as.double( pbar ), as.double( qbar ),
               betamem = as.double( betamem ), deltamem = as.double( deltamem ), sigma2mem = as.double( sigma2mem ),
               gammamem = as.double( gammamem ), kappamem = as.double( kappamem ),
               PACKAGE="zic" )

  # format results
  res <- list( beta = matrix( draws$betamem, n.save, k ), delta = matrix( draws$deltamem, n.save, k ), sigma2 = draws$sigma2mem )
  colnames(res$beta) <- all.vars( formula )
  colnames(res$beta)[1] <- "const"
  colnames(res$delta) <- all.vars( formula )
  colnames(res$delta)[1] <- "const"
  
  # return results
  return( res )
}


zic.ssvs <- function( formula, data, tausq0bar, tausq1bar, omegasq0bar, omegasq1bar, ebar, fbar, pbar, qbar, n.burnin, n.mcmc, n.thin )
{
  # data matrices and dimensions
  mdl <- model.frame( formula, data )
  y <- model.response( mdl )
  X <- model.matrix( formula, mdl )
  n <- dim(X)[1]
  k <- dim(X)[2]

  # initialize memories
  n.save = n.mcmc / n.thin
  betamem <- matrix( 0.0, n.save, k )
  deltamem <- matrix( 0.0, n.save, k )
  sigma2mem <- matrix( 0.0, n.save, 1 )
  gammamem <- matrix( 0.0, n.save, k ) 
  kappamem <- matrix( 0.0, n.save, k )

  # set up prior
  bbarm <- matrix( 0.0, k, 2 )
  bbarm[,1] <- tausq0bar
  bbarm[,2] <- tausq1bar
  dbarm <- matrix( 0.0, k, 2 )
  dbarm[,1] <- omegasq0bar
  dbarm[,2] <- omegasq1bar

  # call C++ code
  draws <- .C( "zic",
               as.integer( y ), as.double( X ),
               as.integer( n ), as.integer( k ),
               as.integer( 1 ),
               as.integer( n.burnin ), as.integer( n.mcmc ), as.integer( n.thin ),
               as.double( bbarm ), as.double( dbarm ), as.double( ebar ), as.double( fbar ), as.double( pbar ), as.double( qbar ),
               betamem = as.double( betamem ), deltamem = as.double( deltamem ), sigma2mem = as.double( sigma2mem ),
               gammamem = as.double( gammamem ), kappamem = as.double( kappamem ),
               PACKAGE="zic" )

  # return results
  res <- list( beta = matrix( draws$betamem, n.save, k ), delta = matrix( draws$deltamem, n.save, k ), sigma2 = draws$sigma2mem, gamma = matrix( draws$gammamem, n.save, k ), kappa = matrix( draws$kappamem, n.save, k ) )

  colnames(res$beta) <- all.vars( formula )
  colnames(res$beta)[1] <- "const"
  colnames(res$delta) <- all.vars( formula )
  colnames(res$delta)[1] <- "const"
  colnames(res$gamma) <- all.vars( formula )
  colnames(res$gamma)[1] <- "const"
  colnames(res$kappa) <- all.vars( formula )
  colnames(res$kappa)[1] <- "const"

  #return results
  return( res )
}

