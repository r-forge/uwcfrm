#' Estimate covariance matrix
#' 
#' @details
#' Add detail for each method
#' 
#' @param R xts or matrix of asset returns
#' @param method The method used to compute the covariance estimate
#' @param control named list of arguments
#' @author Ross Bennett
#' @export
covEstimate <- function(R, 
                        method=c("classic", "mcd", "m_estimate", "mm_estimate", 
                                 "mve", "ogk", "sde", "s_estimate", "lw_shrinkage", "js_shrinkage"),
                        control=list())
  {
  R <- as.matrix(R)
  # Match the method
  method <- match.arg(method)
  
  # Match the function call
  # call <- match.call()
  
  # Switch to select the method
  switch(method,
         # Classic
         classic = {
           tmp_out <- CovClassic(x=R)
           cov <- tmp_out@cov
           },
         # MCD
         mcd = {
           ctrl_fun <- match.fun("CovControlMcd")
           .formals <- formals(ctrl_fun)
           .formals <- modify.args(.formals, control)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_ctrl <- do.call(ctrl_fun, .formals)
           tmp_out <- CovMcd(x=R, control=tmp_ctrl)
           cov <- tmp_out@cov
           },
         # M Estimate
         m_estimate = {
           ctrl_fun <- match.fun("CovControlMest")
           .formals <- formals(ctrl_fun)
           .formals <- modify.args(.formals, control)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_ctrl <- do.call(ctrl_fun, .formals)
           tmp_out <- CovMest(x=R, control=tmp_ctrl)
           cov <- tmp_out@cov
           },
         # MM Estimate
         mm_estimate = {
           ctrl_fun <- match.fun("CovControlMMest")
           .formals <- formals(ctrl_fun)
           .formals <- modify.args(.formals, control)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_ctrl <- do.call(ctrl_fun, .formals)
           tmp_out <- CovMMest(x=R, control=tmp_ctrl)
           cov <- tmp_out@cov
           },
         # MVE
         mve = {
           ctrl_fun <- match.fun("CovControlMve")
           .formals <- formals(ctrl_fun)
           .formals <- modify.args(.formals, control)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_ctrl <- do.call(ctrl_fun, .formals)
           tmp_out <- CovMve(x=R, control=tmp_ctrl)
           cov <- tmp_out@cov
           },
         # OGK
         ogk = {
           if(is.null(control$vrob)) control$vrob <- rrcov:::.vrobGK
           ctrl_fun <- match.fun("CovControlOgk")
           .formals <- formals(ctrl_fun)
           .formals <- modify.args(.formals, control)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_ctrl <- do.call(ctrl_fun, .formals)
           tmp_out <- CovOgk(x=R, control=tmp_ctrl)
           cov <- tmp_out@cov
           },
         # SDE
         sde = {
           ctrl_fun <- match.fun("CovControlSde")
           .formals <- formals(ctrl_fun)
           .formals <- modify.args(.formals, control)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_ctrl <- do.call(ctrl_fun, .formals)
           tmp_out <- CovSde(x=R, control=tmp_ctrl)
           cov <- tmp_out@cov
           },
         # S Estimate
         s_estimate = {
           ctrl_fun <- match.fun("CovControlSest")
           .formals <- formals(ctrl_fun)
           .formals <- modify.args(.formals, control)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_ctrl <- do.call(ctrl_fun, .formals)
           tmp_out <- CovSest(x=R, control=tmp_ctrl)
           cov <- tmp_out@cov
           },
         # Ledoit Wolf Shrinkage
         lw_shrinkage = {
           fun <- match.fun("lwShrink")
           .formals <- formals(fun)
           .formals <- modify.args(.formals, control, x=R)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_out <- do.call(fun, .formals)
           cov <- tmp_out$cov
           },
         # James-Stein Shrinkage
         js_shrinkage = {
           fun <- match.fun("jsShrink")
           .formals <- formals(fun)
           .formals <- modify.args(.formals, control, x=R)
           if(is.pairlist(.formals)) .formals <- as.list(.formals)
           tmp_out <- do.call(fun, .formals)
           cov <- tmp_out$cov
         }
         ) # end switch
  return(structure(list(cov=cov,
                        object=tmp_out)))
}

