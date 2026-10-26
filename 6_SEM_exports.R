#############################################################################################
# Exports from the SEM for plotting the path and for the QTL analysis                      #
#############################################################################################
  # Reads the results of 5_SEM_path_analysis.R (SEM_results.rds) and writes, for one model:
  #   1. The path: every variable (nodes) and every arrow (edges) with its direct effect, raw and standardized,
  #      the regression terms behind it, and a Graphviz .dot file of the path with grey levels by strength
  #   2. The QTL trait sheets: every variable in the model (including the ones only used to calculate others),
  #      one row per strike and one row per fish (fish means), plus a variable list
  #   3. Total effects: every variable's total effect on each performance outcome (all paths), standardized
  #      across strikes and across fish means, for combining with the QTL results
  # Each export has a Markdown file that explains it (SEM_export_path.md, SEM_export_QTL_traits.md, SEM_export_total_effects.md),
  # copied from docsFolder into the export folder, plus export_info.md with this model's details
  #
  # Nothing is refitted: everything comes from the saved results (run script 5 first)

rm(list=ls())


##############
#  Settings  #
##############

  pathData    <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/Data'   ### CHANGE
  resultsFile <- file.path(pathData, 'SEM_results', 'SEM_results.rds')   ### CHANGE  (SEM_results_byFish.rds for the per-fish run)
  outFolder   <- file.path(pathData, 'SEM_exports')                      ### CHANGE  a subfolder per model is made in here
  # Folder with the three Markdown explanations (the folder with the R scripts). NA = don't copy them
  docsFolder  <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/R'   ### CHANGE

  # The model to export (a label from model_comparison.csv). NULL = the bootstrapped model if there is exactly one,
  # otherwise the model with the highest cross-validated R^2 of pressure
  exportModel <- NULL   ### CHANGE  e.g. "tmaxVel_noWidth_staticK+mouthFrontCorrK_static_flowAdditive"

  # Grey levels in the .dot file: weight = |standardized direct effect| / its greyRefQuantile quantile over the fitted arrows
    # (capped at 1). Weight 1 = black, 0 = white. Identity arrows (calculations, no error) are drawn dashed
  greyRefQuantile <- 0.95   ### CHANGE

  # One value per fish for the QTL sheet (and the fish-level standardized total effects): each trait's maximum over the fish's strikes
    # ("max"; the largest value in the trait's usual direction: the most negative for traits that are mostly negative, such as
    # Max_sub_pressure) or the mean ("mean"). Both sheets are written; this one is QTL_traits_per_fish.csv
  perFishSummary <- "max"   ### CHANGE  "max" or "mean"

  # Groups of variables (anatomical structures and other parts of the strike) for the partial effects: each phenotype's effect on
    # each outcome through the paths that pass through any variable of the group. Name = group, value = a regular expression
    # matched against variable names. A phenotype that is itself in a group has all its effect counted as going through that group
    # A path can pass through several groups, so the groups overlap and their effects don't add up to the total
  partialGroups <- c(neurocranium = "Neurocranium|Epaxial|Supraoccipital",
                     hyoid        = "Hyoid|Sternohyoideus|Hypaxial|Urohyal|Ceratohyal",
                     premaxilla   = "Protrusion|premaxilla",
                     maxilla      = "Maxilla|^Oral_",               # maxilla and the oral four-bar (mandible -> maxilla -> premaxilla)
                     mandible     = "Mandible|Protractor_hyoideus",
                     opercular    = "Opercular|Operculum|IOP|RA_process",   # opercular four-bar (neurocranium -> mandible)
                     gape         = "^R_|^A_|^R$|^A$",            # gape radius and area: downstream of mandible and premaxilla (overlaps them)
                     ram          = "U_ram|Cov_R_Uram",
                     timing       = "^T$|^ttpg$|^Time_|^tmax_from|^Relative_|Kinetic_Synchronization",
                     effort       = "^Effort|^Behavior_fish$")   ### CHANGE


#####################
#  Read the results #
#####################

  results  <- readRDS(resultsFile)
  settings <- results$settings
  comparison <- results$comparison
  modelLabels <- setdiff(names(results), c("comparison", "comparisonEquations", "settings"))

  if(is.null(exportModel)) {
    bootLabels <- modelLabels[sapply(modelLabels, function(lab) {
      bt <- results[[lab]]$bootTotal
      !is.null(bt) && any(is.finite(bt))
    })]
    exportModel <- if(length(bootLabels) == 1) bootLabels else
      comparison$label[which.max(comparison$R2_cv_common_Max_sub_pressure)]
  }
  if(!(exportModel %in% modelLabels)) {
    stop("No model '", exportModel, "' in ", basename(resultsFile), ". Models: ", paste(modelLabels, collapse = ", "))
  }
  print(paste0("Exporting model: ", exportModel))

  res  <- results[[exportModel]]
  main <- res$main
  d    <- main$d
  outcomes <- if(!is.null(settings$outcomes)) settings$outcomes else c("Max_sub_pressure", "Total_ingested_volume")
  effortIndicators <- if(!is.null(settings$effortIndicators)) settings$effortIndicators else c("Neurocranium_Adj_vel_mean", "Hyoid_Adj_vel_mean", "T")
  ciLevel <- if(!is.null(settings$ciLevel)) settings$ciLevel else 0.95
  lowQ <- (1 - ciLevel) / 2
  upQ  <- 1 - lowQ

  outNow <- file.path(outFolder, gsub("[^A-Za-z0-9_.+-]", "_", exportModel))
  dir.create(outNow, recursive = TRUE, showWarnings = FALSE)


#############################################
#  Path structure (saved, or rebuilt)       #
#############################################

  # Variables (columns) in an equation's right side, as in script 5
  constantNames <- if(!is.null(settings$constantNames)) settings$constantNames else
    c("rho", "pi", "oralOpenSign", "opercularOpenSign", "maxillaSign", "oralStaticAngle", "opercularStaticAngle")
  varsOf <- function(rhs) {
    v <- setdiff(all.vars(parse(text = rhs)), constantNames)
    v[!sapply(v, exists, mode = "function")]
  }

  pathEqs <- res$pathEqs
  if(!is.null(res$structure)) {
    activeIdentities <- res$structure$activeIdentities
    nodeParents      <- res$structure$nodeParents
    nodeOrder        <- res$structure$nodeOrder
  } else {
    # Older results: rebuild from the equations, as script 5's buildStructure does
    nodeOrder <- rownames(main$B)
    activeIdentities <- settings$identities[intersect(names(settings$identities), setdiff(nodeOrder, names(pathEqs)))]
    nodeParents <- setNames(lapply(nodeOrder, function(v) {
      if(v %in% names(pathEqs)) varsOf(pathEqs[[v]]) else if(v %in% names(activeIdentities)) varsOf(activeIdentities[[v]]) else character(0)
    }), nodeOrder)
  }
  exogenous <- names(nodeParents)[sapply(nodeParents, length) == 0]

  # Rough group of each variable, for colouring the plot (by name; check and edit in path_nodes.csv)
  groupOf <- function(v) {
    if(v %in% outcomes) return("performance")
    if(v %in% c("Effort_score", "Effort_strike", "Behavior_fish")) return("effort")
    if(grepl("_node_zR", v)) return("measurement covariate")
    if(grepl("^U_ram", v)) return("ram")
    if(grepl("^U_flow|^Q_flow|^Cov_", v)) return("flow")
    if(grepl("V_buccal|^Q_buccal|U_flow_ff_predicted", v)) return("buccal volume")
    if(v %in% c("T", "ttpg", "Time_hyoid", "Time_cranial", "tmax_from_Tstart", "Relative_tmax", "Relative_ttpg", "Kinetic_Synchronization")) return("timing")
    if(grepl("_tmax$|_t_mean$|_Max$|_mean$|_vel|_Adj_", v) && !grepl("_rest$|4bar", v)) return("kinematics")
    "morphology"
  }
  roleOf <- function(v) {
    if(v %in% outcomes) return("outcome")
    if(v %in% names(pathEqs)) return("fitted")
    if(v %in% names(activeIdentities)) return("identity")
    if(v %in% exogenous) return("exogenous")
    "other"
  }

  sdOf <- function(x) if(is.numeric(x) && sum(is.finite(x)) > 1) sd(x, na.rm = TRUE) else NA_real_
  fishMeanOf <- function(x) if(is.numeric(x)) tapply(x, d$Fish_ID, function(v) if(all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)) else NULL


##########################
#  1. The path           #
##########################

  # Nodes
  eqSummary <- res$eqSummary
  nodes <- data.frame(variable = nodeOrder, role = sapply(nodeOrder, roleOf), group_guess = sapply(nodeOrder, groupOf),
                      stringsAsFactors = FALSE, row.names = NULL)
  nodes$equation <- sapply(nodeOrder, function(v) if(v %in% names(pathEqs)) pathEqs[[v]] else if(v %in% names(activeIdentities)) activeIdentities[[v]] else "")
  if(!is.null(eqSummary)) {
    nodes$n_rows        <- eqSummary$n[match(nodes$variable, eqSummary$equation)]
    nodes$R2_marginal   <- eqSummary$R2_marginal[match(nodes$variable, eqSummary$equation)]
    nodes$R2_conditional <- eqSummary$R2_conditional[match(nodes$variable, eqSummary$equation)]
  }
  nodes$mean_strikes <- sapply(nodeOrder, function(v) if(v %in% names(d) && is.numeric(d[[v]])) mean(d[[v]], na.rm = TRUE) else NA_real_)
  nodes$sd_strikes   <- sapply(nodeOrder, function(v) if(v %in% names(d)) sdOf(d[[v]]) else NA_real_)

  # Edges: every arrow parent -> child, direct effect at the reference point (raw units of child per unit of parent) and standardized
  edges <- do.call(rbind, lapply(names(nodeParents), function(v) {
    if(length(nodeParents[[v]]) == 0) return(NULL)
    data.frame(from = nodeParents[[v]], to = v, stringsAsFactors = FALSE)
  }))
  edges$edge_type  <- ifelse(edges$to %in% names(pathEqs), "fitted", "identity")
  edges$direct     <- main$B[cbind(edges$to, edges$from)]
  edges$direct_std <- main$directStd[cbind(edges$to, edges$from)]
  if(!is.null(res$bootDirectStd) && any(is.finite(res$bootDirectStd))) {
    keysNow <- paste(edges$from, edges$to, sep = " -> ")
    bootNow <- res$bootDirectStd[, match(keysNow, colnames(res$bootDirectStd)), drop = FALSE]
    edges$direct_std_low  <- apply(bootNow, 2, quantile, lowQ, na.rm = TRUE)
    edges$direct_std_high <- apply(bootNow, 2, quantile, upQ, na.rm = TRUE)
  } else {
    edges$direct_std_low <- edges$direct_std_high <- NA_real_
  }

  # The regression terms of the child's equation that contain the parent (estimate and p, from script 5's coefficient table)
  coefsNow <- res$coefs
  termsOf <- function(from, to) {
    rowsNow <- coefsNow[coefsNow$equation == to & grepl(paste0("(^|[^A-Za-z0-9_.])", from, "([^A-Za-z0-9_.]|$)"), coefsNow$term), , drop = FALSE]
    if(nrow(rowsNow) == 0) return(c(terms = "", min_p = NA_real_))
    c(terms = paste(paste0(rowsNow$term, " = ", signif(rowsNow$estimate, 3), " (p ", signif(rowsNow$p, 2), ")"), collapse = "; "),
      min_p = min(rowsNow$p, na.rm = TRUE))
  }
  termInfo <- t(mapply(termsOf, edges$from, edges$to))
  edges$regression_terms <- termInfo[, "terms"]
  edges$min_p <- suppressWarnings(as.numeric(termInfo[, "min_p"]))

  # Effort measurement model: the latent Effort score and its indicators (standardized loadings), and its fish / strike parts
  loadNow <- main$effort$loadings
  if(!is.null(loadNow) && nrow(loadNow) > 0) {
    effortEdges <- data.frame(from = "Effort_score", to = loadNow$rhs, edge_type = "measurement (loading)", direct = NA_real_,
                              direct_std = loadNow$est.std, direct_std_low = NA_real_, direct_std_high = NA_real_,
                              regression_terms = "", min_p = loadNow$pvalue, stringsAsFactors = FALSE)
    if(any(c("Behavior_fish", "Effort_strike") %in% nodeOrder)) {
      effortEdges <- rbind(effortEdges, data.frame(from = "Effort_score", to = intersect(c("Behavior_fish", "Effort_strike"), nodeOrder),
                                                   edge_type = "decomposition (fish mean / strike deviation)", direct = 1, direct_std = NA_real_,
                                                   direct_std_low = NA_real_, direct_std_high = NA_real_, regression_terms = "", min_p = NA_real_))
    }
    edges <- rbind(edges, effortEdges)
    nodes[nrow(nodes) + 1, ] <- NA
    nodes[nrow(nodes), c("variable", "role", "group_guess", "equation")] <-
      c("Effort_score", "latent", "effort", paste("measurement model:", paste(loadNow$rhs, collapse = " + ")))
  }

  # Strength for plotting: |standardized direct effect| relative to the greyRefQuantile quantile of the fitted arrows, capped at 1
  refStrength <- quantile(abs(edges$direct_std[edges$edge_type == "fitted"]), greyRefQuantile, na.rm = TRUE)
  edges$abs_std <- abs(edges$direct_std)
  edges$plot_weight <- pmin(edges$abs_std / refStrength, 1)
  edges$sign <- ifelse(is.na(edges$direct_std), NA, ifelse(edges$direct_std >= 0, "+", "-"))

  write.csv(nodes, file.path(outNow, "path_nodes.csv"), row.names = FALSE)
  write.csv(edges, file.path(outNow, "path_edges.csv"), row.names = FALSE)
  write.csv(coefsNow, file.path(outNow, "path_equation_terms.csv"), row.names = FALSE)

  # Graphviz .dot (render with dot -Tpdf path_model.dot -o path_model.pdf, or paste into an online Graphviz viewer)
  greyOf <- function(w) if(is.na(w)) "gray70" else paste0("gray", round(100 * (1 - w) * 0.85))   # gray0 = black; up to gray85 for the weakest
  dotLines <- c("digraph SEM {", "  rankdir = LR;", "  node [shape = box, fontsize = 10];",
                sprintf('  "%s" [style = "%s"];', nodes$variable,
                        ifelse(nodes$role == "outcome", "bold", ifelse(nodes$role == "latent", "rounded,dashed", ifelse(nodes$role == "identity", "dashed", "solid")))),
                sprintf('  "%s" -> "%s" [color = "%s", penwidth = %.2f, style = "%s"];', edges$from, edges$to, sapply(edges$plot_weight, greyOf),
                        0.5 + 3 * ifelse(is.na(edges$plot_weight), 0, edges$plot_weight),
                        ifelse(edges$edge_type == "fitted", "solid", "dashed")),
                "}")
  writeLines(dotLines, file.path(outNow, "path_model.dot"))
  print(paste0("Path: ", nrow(nodes), " variables, ", nrow(edges), " arrows (", sum(edges$edge_type == "fitted"), " fitted)"))


#################################
#  2. QTL trait sheets          #
#################################

  # Every variable in the model, including the inputs of identities and the Effort indicators and scores
  traitVars <- unique(c(nodeOrder, effortIndicators, intersect(c("Effort_score", "Effort_strike", "Behavior_fish"), names(d))))
  traitVars <- traitVars[traitVars %in% names(d)]
  # Standard length (script 4's SL, not in the SEM): next to the IDs, for the size corrections
  idCols <- intersect(c("vidName", "Fish_ID", "Family", "Animal.ID", "SL"), names(d))
  if(!("SL" %in% names(d)) || all(is.na(d$SL))) print("No standard length (SL) in the data: set slColumn in script 4 and rerun scripts 4 and 5")

  # Per strike
  perStrike <- d[, c(idCols, traitVars), drop = FALSE]
  write.csv(perStrike, file.path(outNow, "QTL_traits_per_strike.csv"), row.names = FALSE)

  # Per fish: each trait's maximum over the fish's strikes, or its mean (perFishSummary)
    # Maximum = the largest value in the trait's usual direction: the most negative for traits whose median is negative (e.g. pressure)
    # Fish-level variables (morphology) are the same in every strike, so both give the value. Each trait's maximum can come
    # from a different strike
  fishIDs <- sort(unique(d$Fish_ID))
  fishMaxOf <- function(x) {
    if(!is.numeric(x)) return(NULL)
    negNow <- isTRUE(median(x, na.rm = TRUE) < 0)
    tapply(x, d$Fish_ID, function(v) { v <- v[is.finite(v)]; if(length(v) == 0) NA_real_ else if(negNow) min(v) else max(v) })
  }
  perFishOf <- function(summaryNow) {
    out <- data.frame(Fish_ID = fishIDs, stringsAsFactors = FALSE)
    for(v in intersect(c("Family", "Animal.ID", "SL"), names(d))) out[[v]] <- d[[v]][match(fishIDs, d$Fish_ID)]
    out$n_strikes <- as.vector(table(d$Fish_ID)[fishIDs])
    for(v in traitVars) {
      if(is.numeric(d[[v]])) out[[v]] <- as.vector((if(summaryNow == "max") fishMaxOf(d[[v]]) else fishMeanOf(d[[v]]))[fishIDs])
    }
    out
  }
  if(!(perFishSummary %in% c("max", "mean"))) stop("perFishSummary must be 'max' or 'mean'")
  perFish      <- perFishOf(perFishSummary)
  perFishOther <- perFishOf(setdiff(c("max", "mean"), perFishSummary))
  write.csv(perFish, file.path(outNow, "QTL_traits_per_fish.csv"), row.names = FALSE)
  write.csv(perFishOther, file.path(outNow, paste0("QTL_traits_per_fish_", setdiff(c("max", "mean"), perFishSummary), ".csv")), row.names = FALSE)

  # Variable list: role, level (fish = one value per fish; strike = varies between strikes), coverage and spread
  varList <- data.frame(variable = traitVars, stringsAsFactors = FALSE)
  varList$role <- sapply(traitVars, function(v) if(v %in% c("Effort_score")) "latent score" else if(v %in% effortIndicators && !(v %in% nodeOrder)) "Effort indicator" else roleOf(v))
  varList$group_guess <- sapply(traitVars, groupOf)
  varList$level <- sapply(traitVars, function(v) {
    x <- d[[v]]
    if(!is.numeric(x)) return("strike")
    within <- tapply(x, d$Fish_ID, function(z) { z <- z[is.finite(z)]; length(unique(z)) <= 1 })
    if(all(within)) "fish" else "strike"
  })
  varList$n_strikes_with_value <- sapply(traitVars, function(v) sum(!is.na(d[[v]])))
  varList$n_fish_with_value    <- sapply(traitVars, function(v) sum(!is.na(perFish[[v]])))
  varList$sd_strikes    <- sapply(traitVars, function(v) sdOf(d[[v]]))
  varList$fish_value <- ifelse(varList$level == "fish", "value", ifelse(perFishSummary == "max", "maximum over strikes", "mean over strikes"))
  varList$max_direction <- sapply(traitVars, function(v) if(is.numeric(d[[v]]) && isTRUE(median(d[[v]], na.rm = TRUE) < 0)) "most negative" else "largest")
  varList$sd_fish <- sapply(traitVars, function(v) sdOf(perFish[[v]]))
  varList$equation <- sapply(traitVars, function(v) if(v %in% names(pathEqs)) pathEqs[[v]] else if(v %in% names(activeIdentities)) activeIdentities[[v]] else "")
  if("SL" %in% names(d)) {
    varList[nrow(varList) + 1, c("variable", "role", "group_guess", "level", "fish_value")] <- c("SL", "size (not in the SEM)", "size", "fish", "value")
    varList[nrow(varList), c("n_strikes_with_value", "n_fish_with_value", "sd_strikes", "sd_fish")] <-
      c(sum(!is.na(d$SL)), sum(!is.na(perFish$SL)), sdOf(d$SL), sdOf(perFish$SL))
  }
  write.csv(varList, file.path(outNow, "QTL_variable_list.csv"), row.names = FALSE)
  print(paste0("QTL sheets: ", length(traitVars), " variables, ", nrow(perStrike), " strikes, ", nrow(perFish), " fish (per-fish values: ", perFishSummary, ")",
               if("SL" %in% names(d)) paste0("; SL for ", sum(!is.na(perFish$SL)), " fish") else ""))


#################################
#  3. Total effects             #
#################################

  # Total effect of every variable on each outcome = sum over all paths of the products of the direct effects (linearized)
    # total: raw units (outcome units per unit of the variable)
    # total_std: x SD of the variable / SD of the outcome, both across strikes (as in script 5), with bootstrap interval
    # total_std_fish: the same with the SDs across fish of the per-fish values in QTL_traits_per_fish.csv (maxima or means,
      # perFishSummary): the scale of a QTL analysis on those values
    # var_share_multiplier(_fish) = total_std(_fish)^2: share of the outcome's variance per unit share of the variable's variance,
      # if the variable acts on the outcome only through the model's paths (see SEM_export_total_effects.md)
  sdFish <- sapply(nodeOrder, function(v) if(v %in% names(perFish)) sdOf(perFish[[v]]) else NA_real_)
  sdStrike <- sapply(nodeOrder, function(v) if(v %in% names(d)) sdOf(d[[v]]) else NA_real_)
  bootTS <- res$bootTotalStd
  hasBoot <- !is.null(bootTS) && any(is.finite(bootTS))
  totalTable <- do.call(rbind, lapply(outcomes, function(o) {
    causes <- setdiff(nodeOrder, o)
    totalNow <- main$totalE[o, causes]
    directNow <- main$B[o, causes]
    data.frame(outcome = o, variable = causes, role = sapply(causes, roleOf), group_guess = sapply(causes, groupOf),
               has_path = totalNow != 0,
               direct = directNow, indirect = totalNow - directNow, total = totalNow,
               total_std = main$totalStd[o, causes],
               total_std_low  = if(hasBoot) apply(bootTS[, o, causes, drop = FALSE], 3, quantile, lowQ, na.rm = TRUE) else NA_real_,
               total_std_high = if(hasBoot) apply(bootTS[, o, causes, drop = FALSE], 3, quantile, upQ, na.rm = TRUE) else NA_real_,
               total_std_fish = totalNow * sdFish[causes] / sdFish[[o]],
               sd_variable_strikes = sdStrike[causes], sd_variable_fish = sdFish[causes],
               sd_outcome_strikes = sdStrike[[o]], sd_outcome_fish = sdFish[[o]],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  totalTable$var_share_multiplier      <- totalTable$total_std^2
  totalTable$var_share_multiplier_fish <- totalTable$total_std_fish^2
  totalTable <- totalTable[order(totalTable$outcome, -abs(totalTable$total_std_fish)), ]
  write.csv(totalTable, file.path(outNow, "total_effects.csv"), row.names = FALSE)
  print(paste0("Total effects: ", sum(totalTable$has_path), " variable-outcome pairs with a path"))


#################################
#  3b. Partial effects by group #
#################################

  # Effect of x on an outcome through group G = total effect - the effect left when every path that enters G is cut
    # (all arrows into G's variables removed: B[G, ] = 0). If x itself is in G, all its effect goes through G
    # "none": the paths that pass through none of the groups (total - the effect through any group)
  # Standardized effects multiply along paths, so the same calculation on the standardized direct effects gives the standardized
    # partial effects, and the bootstrap samples of the direct effects give their intervals
  totalsOf <- function(Bnow) solve(diag(nrow(Bnow)) - Bnow) - diag(nrow(Bnow))
  groupMembers <- lapply(partialGroups, function(pat) nodeOrder[grepl(pat, nodeOrder)])
  groupMembers$any_group <- unique(unlist(groupMembers))
  for(g in names(partialGroups)) print(paste0("Group ", g, ": ", length(groupMembers[[g]]), " variables (", paste(head(groupMembers[[g]], 6), collapse = ", "),
                                              if(length(groupMembers[[g]]) > 6) ", ..." else "", ")"))
  # Partial effects for one direct-effect matrix: list of matrices [outcome, cause], one per group
  partialOf <- function(Bnow) {
    totNow <- totalsOf(Bnow)
    out <- lapply(groupMembers, function(members) {
      Bcut <- Bnow
      Bcut[members, ] <- 0
      throughNow <- totNow - totalsOf(Bcut)
      inGroup <- colnames(Bnow) %in% members
      throughNow[, inGroup] <- totNow[, inGroup]
      throughNow[outcomes, , drop = FALSE]
    })
    out$none <- totNow[outcomes, , drop = FALSE] - out$any_group
    # Rounding: differences that are 0 up to floating-point error are set to 0
    tolNow <- 1e-9 * abs(totNow[outcomes, , drop = FALSE]) + 1e-14
    lapply(out, function(m) { m[abs(m) <= tolNow] <- 0; m })
  }
  partRaw <- partialOf(main$B)
  partStd <- partialOf(main$directStd)

  # Bootstrap: rebuild the standardized direct-effect matrix of each sample from the saved arrows
  bootDS <- res$bootDirectStd
  partBoot <- NULL
  if(!is.null(bootDS) && any(is.finite(bootDS))) {
    edgeKeysNow <- strsplit(colnames(bootDS), " -> ", fixed = TRUE)
    edgeFrom <- sapply(edgeKeysNow, `[`, 1)
    edgeTo   <- sapply(edgeKeysNow, `[`, 2)
    okRows <- which(apply(bootDS, 1, function(z) any(is.finite(z))))
    partBoot <- lapply(okRows, function(k) {
      Bk <- matrix(0, nrow(main$B), ncol(main$B), dimnames = dimnames(main$B))
      valsNow <- bootDS[k, ]
      valsNow[!is.finite(valsNow)] <- 0
      Bk[cbind(edgeTo, edgeFrom)] <- valsNow
      tryCatch(partialOf(Bk), error = function(e) NULL)
    })
    partBoot <- partBoot[!sapply(partBoot, is.null)]
  }

  partialTable <- do.call(rbind, lapply(outcomes, function(o) {
    causes <- setdiff(nodeOrder, o)
    do.call(rbind, lapply(c(names(partialGroups), "none"), function(g) {
      rawNow <- partRaw[[g]][o, causes]
      stdNow <- partStd[[g]][o, causes]
      bootMat <- if(length(partBoot) > 0) sapply(partBoot, function(pb) pb[[g]][o, causes]) else NULL
      data.frame(outcome = o, variable = causes, group = g, role = sapply(causes, roleOf), group_guess = sapply(causes, groupOf),
                 variable_in_group = if(g == "none") FALSE else causes %in% groupMembers[[g]],
                 partial = rawNow, partial_std = stdNow,
                 partial_std_low  = if(!is.null(bootMat)) apply(bootMat, 1, quantile, lowQ, na.rm = TRUE) else NA_real_,
                 partial_std_high = if(!is.null(bootMat)) apply(bootMat, 1, quantile, upQ, na.rm = TRUE) else NA_real_,
                 partial_std_fish = rawNow * sdFish[causes] / sdFish[[o]],
                 share_of_total = ifelse(main$totalE[o, causes] != 0, rawNow / main$totalE[o, causes], NA_real_),
                 stringsAsFactors = FALSE, row.names = NULL)
    }))
  }))
  partialTable$var_share_multiplier_fish <- partialTable$partial_std_fish^2
  partialTable <- partialTable[partialTable$partial != 0 | partialTable$group == "none", ]
  write.csv(partialTable, file.path(outNow, "partial_effects_by_group.csv"), row.names = FALSE)

  # Wide versions (one row per variable, one column per group), fish-level standardized, for the per-structure QTL maps
  for(o in outcomes) {
    wideNow <- data.frame(variable = setdiff(nodeOrder, o), stringsAsFactors = FALSE)
    wideNow$total_std_fish <- totalTable$total_std_fish[totalTable$outcome == o][match(wideNow$variable, totalTable$variable[totalTable$outcome == o])]
    for(g in c(names(partialGroups), "none")) {
      rowsG <- partialTable[partialTable$outcome == o & partialTable$group == g, ]
      wideNow[[g]] <- rowsG$partial_std_fish[match(wideNow$variable, rowsG$variable)]
      wideNow[[g]][is.na(wideNow[[g]])] <- 0
    }
    write.csv(wideNow, file.path(outNow, paste0("partial_effects_wide_", o, ".csv")), row.names = FALSE)
  }
  write.csv(data.frame(group = names(groupMembers), pattern = c(partialGroups, "(union of the groups)"),
                       variables = sapply(groupMembers, paste, collapse = ", ")),
            file.path(outNow, "partial_effect_groups.csv"), row.names = FALSE)
  print(paste0("Partial effects: ", length(partialGroups), " groups + none (", if(length(partBoot) > 0) paste0(length(partBoot), " bootstrap samples") else "no bootstrap", ")"))


#################################
#  Model details and the docs   #
#################################

  rowNow <- comparison[comparison$label == exportModel, , drop = FALSE]
  r2Now <- res$r2
  infoLines <- c(paste0("# SEM export: ", exportModel), "",
                 paste0("Exported ", format(Sys.time(), "%Y-%m-%d %H:%M"), " from `", basename(resultsFile), "` by `6_SEM_exports.R`."), "",
                 "## Model", "",
                 paste0("- Time basis: ", res$variant$timeBasis, "; width model: ", res$variant$widthModel, "; flow version: ", res$variant$flowVersion,
                        "; KT: ", res$variant$ktVersion, "; structure: ", res$variant$structure,
                        if(!is.null(res$variant$pressure)) paste0("; pressure: ", res$variant$pressure, " (\"corrected\" = script 4's pressureVersion)") else ""),
                 paste0("- Rows (strikes): ", nrow(d), "; fish: ", length(fishIDs)),
                 paste0("- Path equations: ", length(pathEqs), "; identities: ", length(activeIdentities), "; exogenous variables: ", length(exogenous)),
                 paste0("- Bootstrap samples: ", if(hasBoot) sum(apply(bootTS, 1, function(z) any(is.finite(z)))) else 0,
                        " (intervals are ", round(100 * ciLevel), "% percentile intervals)"),
                 if(nrow(rowNow) == 1) paste0("- Fisher's C = ", signif(rowNow$fisherC, 4), ", df = ", rowNow$fisherC_df, ", p = ", signif(rowNow$fisherC_p, 3)) else NULL,
                 "", "## Variance explained in the outcomes", "",
                 "| outcome | R2 local | R2 whole model | R2 whole model, cross-validated |", "|---|---|---|---|",
                 sprintf("| %s | %s | %.3f | %.3f |", r2Now$outcome, ifelse(is.na(r2Now$R2_local), "identity", sprintf("%.3f", r2Now$R2_local)),
                         r2Now$R2_whole, r2Now$R2_whole_cv),
                 if(!is.null(res$r2Decomposition)) c("", "## Whole-model R2 by group of starting variables", "",
                   "`unique_share` = R2 lost when the group is held at its mean (variance only that group explains); `R2_alone` = the group by itself. In-sample, from the main fit.", "",
                   "| outcome | group | unique share | R2 alone |", "|---|---|---|---|",
                   sprintf("| %s | %s | %.3f | %.3f |", res$r2Decomposition$outcome, res$r2Decomposition$group,
                           res$r2Decomposition$unique_share, res$r2Decomposition$R2_alone)) else NULL,
                 "", "## Files", "",
                 "- Path: `path_nodes.csv`, `path_edges.csv`, `path_equation_terms.csv`, `path_model.dot` (see `SEM_export_path.md`)",
                 "- QTL traits: `QTL_traits_per_strike.csv`, `QTL_traits_per_fish.csv`, `QTL_variable_list.csv` (see `SEM_export_QTL_traits.md`)",
                 "- Total effects: `total_effects.csv` (see `SEM_export_total_effects.md`)",
                 "- Partial effects by anatomical structure / group: `partial_effects_by_group.csv`, `partial_effects_wide_<outcome>.csv`, `partial_effect_groups.csv` (see `SEM_export_total_effects.md`)")
  writeLines(infoLines, file.path(outNow, "export_info.md"))

  if(!is.na(docsFolder)) {
    for(docNow in c("SEM_export_path.md", "SEM_export_QTL_traits.md", "SEM_export_total_effects.md")) {
      srcNow <- file.path(docsFolder, docNow)
      if(file.exists(srcNow)) file.copy(srcNow, file.path(outNow, docNow), overwrite = TRUE) else print(paste0("Not found (not copied): ", srcNow))
    }
  }
  print(paste0("Saved in: ", normalizePath(outNow)))
