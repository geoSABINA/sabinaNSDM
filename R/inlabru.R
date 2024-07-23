#' @name NSDM.inlabru
#'
#' @title Nested species distribution modelling by means of hierarchical Bayesian models
#'
#' @description Estimation of presence probability by means of hierarchical 
#' Bayesian models.
#'
#' @param nsdm_obj An object of class \class{nsdm.vinput} with all the required data to fit the model.
#'
#' @param output A character. Either "intensity" or "probability".
#'
#' @return A list with three elements: \code{fit} (mode fit), \code{pred}
#' (prediction from the model) and \code{pred_sp} (prediction of the
#' spatial effect.
#'
#' @details TBC.
#'
#' @seealso \code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}
#'
#' @examples
#'
#' # FIXME: Add these packages to DESCRIPTION?
#' library("glmnet")
#' library("stringi")
#'
#' # Example from GitHub
#' 
#' SpeciesName <- "Fagus.sylvativa"
#' 
#' # Species occurrences
#' data(Fagus.sylvatica.xy.global, package = "sabinaNSDM")
#' spp.data.global <- Fagus.sylvatica.xy.global
#' data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM")
#' spp.data.regional <- Fagus.sylvatica.xy.regional
#' 
#' data(expl.var.global, package = "sabinaNSDM")
#' data(expl.var.regional, package = "sabinaNSDM")
#' expl.var.global <- terra::unwrap(expl.var.global)
#' expl.var.regional <- terra::unwrap(expl.var.regional)
#'
#' # new escenarios
#' data(new.env, package = "sabinaNSDM")
#' new.env <- terra::unwrap(new.env)
#'
#' nsdm_input <- NSDM.InputData(SpeciesName = SpeciesName,
#'   spp.data.global = Fagus.sylvatica.xy.global,
#'   spp.data.regional = Fagus.sylvatica.xy.regional,
#'   expl.var.global = expl.var.global,
#'   expl.var.regional = expl.var.regional,
#'   new.env = new.env,
#'   new.env.names = "scenario1",
#'   Background.Global = NULL,
#'   Background.Regional = NULL,
#'   Absences.Global = NULL,
#'   Absences.Regional = NULL)
#' 
#' nsdm_finput <- NSDM.FormattingData(nsdm_input,
#'   nPoints = 100, # number of background points
#'   Min.Dist.Global = "resolution",
#'   Min.Dist.Regional = "resolution",
#'   Background.method = "random", # method “random" or "stratified” to generate background points 
#'   save.output = TRUE) #save outputs locally
#' 
#' nsdm_selvars <- NSDM.SelectCovariates(nsdm_finput,
#'   maxncov.Global = 2,   # Max number of covariates to be selected at the global scale
#'   maxncov.Regional = 2, # Max number of covariates to be selected at the regional scale
#'   corcut = 0.7, #  correlation threshold
#'   algorithms = c("glm","gam","rf"),
#'   ClimaticVariablesBands = NULL, # covariate bands to be excluded in the covariate selection at the regional scale
#'   save.output = TRUE)
#' 
#' # Probability of presence
#' nsdm_inlabru_prob <- NSDM.inlabru(nsdm_selvars)
#' # Intensity
#' nsdm_inlabru_int <- NSDM.inlabru(nsdm_selvars, "intensity")
#' 
#' @export

NSDM.inlabru <- function(nsdmobj, output = "probability") {

  # Checks
  if(!(output %in% c("probability", "intensity")))
    stop("Wrong value for 'output'.")

  # Regional point pattern
  pp_regional <- st_as_sf(nsdmobj$SpeciesData.XY.Regional,
    coords = c("x", "y"))

  # Global point pattern
  pp_global <- st_as_sf(nsdmobj$SpeciesData.XY.Global,
    coords = c("x", "y"))

  # Covariates
  sp_covglo <- unwrap(nsdmobj$IndVar.Global.Selected)
  sp_covreg <- unwrap(nsdmobj$IndVar.Regional.Selected)

  # SET CRS
  # FIXME: Set "correct" CRS
  st_crs(pp_global) <- st_crs(sp_covglo)
  st_crs(pp_regional) <- st_crs(sp_covglo)

  # -- Using the points plus the points in the regional raster (so that
  #      prediction can be done).
  aux <- st_as_sf(st_transform(st_as_sfc(st_as_sf(as.points(sp_covreg))),
    st_crs(pp_global)))
  st_geometry(aux) <- "geometry"
  aux <- rbind(pp_global, aux)
  bdy_global <- st_convex_hull(st_union(aux))


  bdy_regional <- st_as_sf(rasterToPolygons(raster(sp_covreg)))
  bdy_regional <- st_union(st_make_valid(bdy_regional))


  # Enlarge a bit
  # TODO: THis can be set to a fraction of the study region
  bdy_global <- st_buffer(bdy_global, 0.01)


  # SET CRS
  # FIXME: Set "correct" CRS
  st_crs(bdy_global) <- st_crs(sp_covglo)
  st_crs(bdy_regional) <- st_crs(sp_covreg)

  # Display covariates
  ggplot() + geom_spatraster(data = sp_covglo, aes(fill = bio4)) +
    geom_sf(data = bdy_regional) + geom_sf(data = bdy_global) #+ gg(mesh)

  # Display raster and covariates
  # Display covariates
  ggplot() + geom_spatraster(data = sp_covglo, aes(fill = bio4)) +
    xlim(-10, 5) + ylim(35, 45) +
    geom_sf(data = bdy_regional) + geom_sf(data = bdy_global) + 
    geom_sf(data = pp_global, color = "black") +
    geom_sf(data = pp_regional, color = "grey")

  # Define mesh
  # FIXME: Compute mesh and pass it ti this function? 
  mesh <- fm_mesh_2d(
    boundary = bdy_global,
    max.edge = 4 * c(0.5, 1),
    offset = c(0.25, 0.5)
  )

  # FIXME: Set "correct" CRS
  fm_crs(mesh) <- st_crs(sp_covglo)

  # Display mesh
  ggplot() + geom_spatraster(data = sp_covglo, aes(fill = bio4)) +
    #xlim(-10, 5) + ylim(35, 45) +
    geom_sf(data = bdy_regional) + 
    gg(mesh) 

  # Define spde
  matern <- inla.spde2.pcmatern(mesh,
    prior.range = c(5, 0.01),
    prior.sigma = c(1, 0.01)
  )


  # Create bit of covariates for model formula
  # obj: Object with required data
  # spobjglo: Name of spatial object with layers (global)
  # spobjreg: Name of spatial object with layers (regional)
  # FIXME: Check whether regional variables are included also at regional 
  #   level when these are not avaialble at the regional level.
  fcov <- function(obj, spobjglo, spobjreg) {

    # Components
    vars <- unique(c(obj$Selected.Variables.Global, obj$Selected.Variables.Regional))
    covars <- paste0("covars", vars)

    cmp1 <- paste(paste(vars, "(1)", sep = ""), collapse = " + ")

    cmpglobal <- sapply(obj$Selected.Variables.Global, function(X) {
      paste0(X, "GL", "(main = ", spobjglo, ", main_layer = \"", X,
        "\", model = \"const\")")
    })
    cmpglobal <- paste(cmpglobal, collapse = " + ")

    cmpregional <- sapply(obj$Selected.Variables.Regional, function(X) {
      paste0(X, "RE", "(main = ", spobjreg, ", main_layer = \"", X,
        "\", model = \"const\")")
    })
    cmpregional <- paste(cmpregional, collapse = " + ")

    cmp <- paste(c(cmp1, cmpglobal, cmpregional), collapse = " + ")

 
    fglobal <- sapply(obj$Selected.Variables.Global, function(X) {
      paste0(X, " * ", X, "GL")
    })
    fglobal <- paste(fglobal, collapse = " + ")

    fregional <- sapply(obj$Selected.Variables.Regional, function(X) {
      paste0(X, " * ", X, "RE")
    })
    fregional <- paste(fregional, collapse = " + ")

    return(list(cmp = cmp,
      like = list(fglobal = fglobal, fregional = fregional))
    )
  }

  cmp_cov <- fcov(nsdmobj, "sp_covglo", "sp_covreg")

  # Fit model
  cmp <- as.formula(
    paste0("~ IGlobal(1) + IRegional(1) + spatial(geometry, model = matern)",
      " + ",
      cmp_cov$cmp
    )
  )


  #Update 
  if( output == "probability") {
    # Add presence-absence data
    pp_regional$presence <- 1

    pseudo_regional <- st_as_sf(nsdm_global$Background.XY.Regional,
      coords = c("x", "y"))
    pseudo_regional$presence <- 0
    st_crs(pseudo_regional) <- st_crs(pp_regional)

    pp_regional <- rbind(pp_regional, pseudo_regional)

    # Global point pattern
    pp_global$presence <- 1

    pseudo_global <- st_as_sf(nsdm_global$Background.XY.Global,
      coords = c("x", "y"))
    pseudo_global$presence <- 0
    st_crs(pseudo_global) <- st_crs(pp_global)

    pp_global <- rbind(pp_global, pseudo_global)

    # Likelihoods

    lik_global <- inlabru::like(family = "binomial",
      formula = as.formula(paste0("presence ~ IGlobal + spatial + ",
        cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = st_as_sf(bdy_global),
      domain = list(geometry = mesh)
    )

    lik_regional <- inlabru::like(family = "binomial",
      formula = as.formula(paste0("presence ~ IRegional + spatial + ",
        cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = st_as_sf(bdy_regional),
      domain = list(geometry = mesh)
    )

    prerd_formula <- as.formula(
      paste0(" ~ 1/ (1 + exp(-(IRegional + spatial + ",
        cmp_cov$like$fregional, ")))"))

  } else {
    # FIXME: Check all the parts with output == "intensity"
    lik_global <- inlabru::like(family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal + spatial + ",
        cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = st_as_sf(bdy_global),
      domain = list(geometry = mesh)
    )

    lik_regional <- inlabru::like(family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional + spatial + ",
        cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = st_as_sf(bdy_regional),
      domain = list(geometry = mesh)
    )

    pred_formula <- as.formula(
        paste0(" ~ exp(IRegional + spatial + ", cmp_cov$like$fregional, ")")
      )
  }

  # Fit model
  fit1 <- bru(cmp,
    lik_global, lik_regional
  )

  # Predict results (?)
  #pred.df <- fm_pixels(mesh, mask = bdy_regional, dims = c(20, 30))
  pred.df <- st_as_sf(as.points(sp_covreg))
  pred1 <- predict(fit1, pred.df, pred_formula)
  

  # Display prediction
  library("ggplot2")
  ggplot() +
    gg(pred1, geom = "tile") + # gg() with sf points and geom = "tile" plots a raster
    gg(st_as_sf(bdy_regional), alpha = 0, lwd = 2) +
    gg(pp_regional, color = "DarkGreen")

  # Spatial effect
  pred1_sp <- predict(fit1, pred.df, ~ spatial)

  ggplot() +
    gg(pred1_sp, geom = "tile") + # gg() with sf points and geom = "tile" plots a raster
    gg(st_as_sf(bdy_regional), alpha = 0, lwd = 2) +
    gg(pp_regional, color = "DarkGreen")

  return(list(fit = fit1, pred = pred1, pred_sp = pred1_sp))
}
