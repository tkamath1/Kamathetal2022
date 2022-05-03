library(rliger)
library(Seurat)

`%notin%` <- Negate(`%in%`)
## Liger extra utilities

AddMito <- function(object, species){
  if (species == 'mouse'){
    pattern.use <- '^mt-'
  }
  else if (species == 'human'){
    pattern.use <- '^MT-'
  }
  if (class(object)[1] == 'rliger'){
    mito.genes <- Reduce(intersect,lapply(object@raw.data,function(q){
      grep(pattern = pattern.use, x = rownames(q), value = TRUE)
    }))
    percent.mito <- unlist(lapply(object@raw.data, function(x) {
      Matrix::colSums(x[mito.genes, ])/Matrix::colSums(x)
    }))
    object@cell.data$percent.mito <- percent.mito
    return(object)
  }
  else
    if (class(object)[1] == 'seurat'){
      mito.genes <- grep(pattern = pattern.use, x = rownames(object@raw.data), value = TRUE)
      percent.mito <-Matrix::colSums(object@raw.data[mito.genes, ])/Matrix::colSums(object@raw.data)
      object@meta.data$percent.mito <- percent.mito
      return(object)
    }
}

## Do MAST

##################
seuratToLiger <- function(objects, combined.seurat = F, names = "use-projects", meta.var = NULL,
                          assays.use = NULL, raw.assay = "RNA", remove.missing = T, renormalize = T,
                          use.seurat.genes = T, num.hvg.info = NULL, use.idents = T, use.tsne = T,
                          cca.to.H = F) {
  
  # Remind to set combined.seurat
  if ((typeof(objects) != "list") & (!combined.seurat)) {
    stop("Please pass a list of objects or set combined.seurat = T")
  }
  # Get Seurat versions
  if (typeof(objects) != "list") {
    version <- package_version(objects@version)$major
  } else {
    version <- sapply(objects, function(x) {
      package_version(x@version)$major
    })
    if (min(version) != max(version)) {
      stop("Please ensure all Seurat objects have the same major version.")
    } else {
      version <- version[1]
    }
  }
  
  # Only a single seurat object expected if combined.seurat
  if (combined.seurat) {
    if ((is.null(meta.var)) & (is.null(assays.use))) {
      stop("Please provide Seurat meta.var or assays.use to use in identifying individual datasets.")
    }
    if (!is.null(meta.var)) {
      # using meta.var column as division split
      if (version > 2) {
        # if integrated assay present, want to make sure to use original raw data
        object.raw <- GetAssayData(objects, assay = raw.assay, slot = "counts")
      } else {
        object.raw <- objects@raw.data
      }
      if (nrow(objects@meta.data) != ncol(object.raw)) {
        cat("Warning: Mismatch between meta.data and raw.data in this Seurat object. \nSome cells",
            "will not be assigned to a raw dataset. \nRepeat Seurat analysis without filters to",
            "allow all cells to be assigned.\n")
      }
      raw.data <- lapply(unique(objects@meta.data[[meta.var]]), function(x) {
        cells <- rownames(objects@meta.data[objects@meta.data[[meta.var]] == x, ])
        object.raw[, cells]
      })
      names(raw.data) <- unique(objects@meta.data[[meta.var]])
    } else {
      # using different assays in v3 object
      raw.data <- lapply(assays.use, function(x) {
        GetAssayData(objects, assay = x, slot = "counts")
      })
      names(raw.data) <- assays.use
    }
    
    if (version > 2) {
      var.genes <- VariableFeatures(objects)
      idents <- Idents(objects)
      if (is.null(objects@reductions$tsne)) {
        cat("Warning: no t-SNE coordinates available for this Seurat object.\n")
        tsne.coords <- NULL
      } else {
        tsne.coords <- objects@reductions$tsne@cell.embeddings
      }
    } else {
      # Get var.genes
      var.genes <- objects@var.genes
      # Get idents/clusters
      idents <- objects@ident
      # Get tsne.coords
      if (is.null(objects@dr$tsne)) {
        cat("Warning: no t-SNE coordinates available for this Seurat object.\n")
        tsne.coords <- NULL
      } else {
        tsne.coords <- objects@dr$tsne@cell.embeddings
      }
    }
  } else {
    # for multiple Seurat objects
    raw.data <- lapply(objects, function(x) {
      if (version > 2) {
        # assuming default assays have been set for each v3 object
        GetAssayData(x, slot = "counts")
      } else {
        x@raw.data
      }
    })
    names(raw.data) <- lapply(seq_along(objects), function(x) {
      if (identical(names, "use-projects")) {
        if (!is.null(meta.var)) {
          cat("Warning: meta.var value is set - set names = 'use-meta' to use meta.var for names.\n")
        }
        objects[[x]]@project.name
      } else if (identical(names, "use-meta")) {
        if (is.null(meta.var)) {
          stop("Please provide meta.var to use in naming individual datasets.")
        }
        objects[[x]]@meta.data[[meta.var]][1]
      } else {
        names[x]
      }
    })
    # tsne coords not very meaningful for separate objects
    tsne.coords <- NULL
    
    if (version > 2) {
      var.genes <- Reduce(union, lapply(objects, function(x) {
        VariableFeatures(x)
      }))
      # Get idents, label by dataset
      idents <- unlist(lapply(seq_along(objects), function(x) {
        idents <- rep("NA", ncol(raw.data[[x]]))
        names(idents) <- colnames(raw.data[[x]])
        idents[names(Idents(objects[[x]]))] <- as.character(Idents(objects[[x]]))
        idents <- paste0(names(raw.data)[x], idents)
      }))
      idents <- factor(idents)
    } else {
      var.genes <- Reduce(union, lapply(objects, function(x) {
        if (!is.null(num.hvg.info)) {
          rownames(head(x@hvg.info, num.hvg.info))
        } else {
          x@var.genes
        }
      }))
      # Get idents, label by dataset
      idents <- unlist(lapply(seq_along(objects), function(x) {
        idents <- rep("NA", ncol(objects[[x]]@raw.data))
        names(idents) <- colnames(objects[[x]]@raw.data)
        idents[names(objects[[x]]@ident)] <- as.character(objects[[x]]@ident)
        idents <- paste0(names(raw.data)[x], idents)
      }))
      idents <- factor(idents)
    }
  }
  new.liger <- createLiger(raw.data = raw.data, remove.missing = remove.missing)
  if (renormalize) {
    new.liger <- rliger::normalize(new.liger)
  }
  if (use.seurat.genes) {
    # Include only genes which appear in all datasets
    for (i in 1:length(new.liger@raw.data)) {
      var.genes <- intersect(var.genes, rownames(new.liger@raw.data[[i]]))
      # Seurat has an extra CheckGenes step which we can include here
      # Remove genes with no expression anywhere
      var.genes <- var.genes[Matrix::rowSums(new.liger@raw.data[[i]][var.genes, ]) > 0]
      var.genes <- var.genes[!is.na(var.genes)]
    }
    
    new.liger@var.genes <- var.genes
  }
  if (use.idents) {
    new.liger@clusters <- idents
  }
  if ((use.tsne) & (!is.null(tsne.coords))) {
    new.liger@tsne.coords <- tsne.coords
  }
  # Get CCA loadings if requested
  if (cca.to.H & combined.seurat) {
    if (version > 2) {
      cat("Warning: no CCA loadings available for Seurat v3 objects.\n")
      return(new.liger)
    }
    if (is.null(objects@dr$inmf)) {
      cat("Warning: no CCA loadings available for this Seurat object.\n")
    } else {
      new.liger@H <- lapply(unique(objects@meta.data[[meta.var]]), function(x) {
        cells <- rownames(objects@meta.data[objects@meta.data[[meta.var]] == x, ])
        objects@dr$inmf@cell.embeddings[cells, ]
      })
      new.liger@H <- lapply(seq_along(new.liger@H), function(x) {
        addMissingCells(new.liger@raw.data[[x]], new.liger@H[[x]])
      })
      names(new.liger@H) <- names(new.liger@raw.data)
    }
    if (is.null(objects@dr$inmf)) {
      cat("Warning: no aligned CCA loadings available for this Seurat object.\n")
    } else {
      new.liger@H.norm <- objects@dr$inmf@cell.embeddings
      new.liger@H.norm <- addMissingCells(Reduce(rbind, new.liger@H), new.liger@H.norm,
                                          transpose = T)
    }
  }
  return(new.liger)
}

addMissingCells <- function(matrix1, matrix.subset, transpose = F) {
  if (transpose) {
    matrix1 <- t(matrix1)
  }
  if (ncol(matrix1) != nrow(matrix.subset)) {
    extra <- matrix(NA, nrow = ncol(matrix1) - nrow(matrix.subset),
                    ncol = ncol(matrix.subset))
    colnames(extra) <- colnames(matrix.subset)
    rownames(extra) <- setdiff(colnames(matrix1), rownames(matrix.subset))
    matrix.subset <- rbind(matrix.subset, extra)
  }
  return(matrix.subset)
}




MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                 verbose = FALSE, save_models = FALSE, save_model_dir = NULL) {
  # Check inputs
  if (is.factor(dataset[[contrast]]) == FALSE) {
    stop("Specified contrast term is not coded as a factor in dataset")
  }
  
  # Generate design matrix from cluster assignments
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  
  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  
  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      # For now, do not allow models without mixed effects terms
      #stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }
  
  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  
  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    test_cluster <- sprintf("`%s`", test_cluster)
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    #Run null and full mixed-effects models
    #null_model <- lme4::glmer(formula = null_fm, data = dataset,
    #                          family = binomial,
    #                          control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial,control = glmerControl(optimizer = 'bobyqa'))
    cluster_models[[i]]$output <- broom.mixed::tidy(full_model)
    #model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    #contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    #contrast_ci <- confint.merMod(full_model, method = "Wald",
    #                              parm = contrast_lvl2)
    # Save model objects to list
    #cluster_models[[i]]$null_model <- null_model
    #cluster_models[[i]]$full_model <- full_model
    #cluster_models[[i]]$model_lrt <- model_lrt
    #cluster_models[[i]]$confint <- contrast_ci
  }
  
  # Organize results into output dataframe
  #output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
  #                     size = colSums(designmat))
  #output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  #output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  #output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  #output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))
  
  # Return MASC results and save models if specified
  if (save_models == TRUE) {
    saveModelObj(cluster_models, save_dir = save_model_dir)
    return(output)
  } else {
    return(cluster_models)
  }
}


## Next few functions all have to do with plotFactor_new (with png option)
get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


plotGene_ALL = function (object, gene, methylation.indices = NULL, pt.size = 0.1, 
                         min.clip = 0, max.clip = 1, points.only = F, low.col = "yellow", 
                         high.col = "red", return.plots = F,option="plasma",zero_color="#F5F5F5",by.dataset=T,
                         reduction.use = 'tsne') 
{
  if (class(object)[[1]] == "liger") {
    gene_vals <- c()
    for (i in 1:length(object@norm.data)) {
      if (i %in% methylation.indices) {
        tmp <- object@norm.data[[i]][gene, ]
        max_v <- quantile(tmp, probs = max.clip, na.rm = T)
        min_v <- quantile(tmp, probs = min.clip, na.rm = T)
        tmp[tmp < min_v & !is.na(tmp)] <- min_v
        tmp[tmp > max_v & !is.na(tmp)] <- max_v
        gene_vals <- c(gene_vals, tmp)
      }
      else {
        if (gene %in% rownames(object@norm.data[[i]])) {
          gene_vals_int <- log2(10000 * object@norm.data[[i]][gene,] + 1)
          gene_vals_int[gene_vals_int==0]=NA
        }
        else {
          gene_vals_int <- rep(list(NA), ncol(object@norm.data[[i]]))
          names(gene_vals_int) <- colnames(object@norm.data[[i]])
        }
        gene_vals <- c(gene_vals, gene_vals_int)
      }
    }
    gene_df <- data.frame(object@tsne.coords)
    rownames(gene_df) <- names(object@clusters)
    gene_df$gene <- as.numeric(gene_vals[rownames(gene_df)])
    colnames(gene_df) <- c("tSNE1", "tSNE2", "gene")
    gene_plots <- list()
    if(by.dataset){
      for (i in 1:length(object@norm.data)) {
        gene_df.sub <- gene_df[rownames(object@scale.data[[i]]),]
        max_v <- max(gene_df.sub["gene"], na.rm = T)
        min_v <- min(gene_df.sub["gene"], na.rm = T)
        midpoint <- (max_v - min_v)/2
        plot_i <- (ggplot(gene_df.sub, aes_string(x = "tSNE1", 
                                                  y = "tSNE2", color = "gene")) + geom_point_rast(size = pt.size) + 
                     scale_color_viridis(option=option,direction=-1,na.value=zero_color) + labs(col = gene) + 
                     ggtitle(names(object@scale.data)[i]))
        gene_plots[[i]] <- plot_i
      }
    }
    else
    {
      plot_i <- ggplot(gene_df, aes_string(x = "tSNE1", 
                                           y = "tSNE2", color = "gene")) + geom_point_rast(size = pt.size) + 
        scale_color_viridis(option=option,direction=-1,na.value=zero_color) + labs(col = gene) + ggtitle("Full dataset")
      gene_plots[[1]] <- plot_i
    }
    if (points.only) {
      for (i in 1:length(gene_plots)) {
        gene_plots[[i]] <- gene_plots[[i]] + theme(axis.line = element_blank(), 
                                                   axis.text.x = element_blank(), axis.text.y = element_blank(), 
                                                   axis.ticks = element_blank(), axis.title.x = element_blank(), 
                                                   axis.title.y = element_blank(), 
                                                   panel.background = element_blank(), panel.border = element_blank(), 
                                                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                   plot.background = element_blank())
        gene_plots[[i]] <- gene_plots[[i]] + ggtitle(as.character(gene))
      }
    }
    if (return.plots) {
      return(gene_plots)
    }
    else {
      for (i in 1:length(gene_plots)) {
        print(gene_plots[[i]])
      }
    }
    
  } else if (class(object)[[1]] == "seurat") {
    
    if (gene %in% rownames(object@data)) {
      gene_vals <- object@data[gene,]
      gene_vals[gene_vals==0]=NA
    }
    else {
      gene_vals <- rep(list(NA), ncol(object@data))
      names(gene_vals) <- colnames(object@data)
    }
    
    gene_df <- data.frame(object@dr[[reduction.use]]@cell.embeddings)
    rownames(gene_df) <- names(object@ident)
    gene_df$gene <- as.numeric(gene_vals[rownames(gene_df)])
    colnames(gene_df) <- c(paste0(reduction.use,"_1"), paste0(reduction.use,"_2"), "gene")
    gene_plots <- list()
    plot_i <- ggplot(gene_df, aes_string(x = colnames(gene_df)[1], 
                                         y = colnames(gene_df)[2], color = "gene")) + geom_point(size = pt.size) + 
      scale_color_viridis(option=option,direction=-1,na.value=zero_color) + labs(col = gene) + ggtitle(gene)
    
    if (points.only) {
      plot_i <- plot_i + theme(axis.line = element_blank(), 
                               axis.text.x = element_blank(), axis.text.y = element_blank(), 
                               axis.ticks = element_blank(), axis.title.x = element_blank(), 
                               axis.title.y = element_blank(), 
                               panel.background = element_blank(), panel.border = element_blank(), 
                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                               plot.background = element_blank())#, plot.title = element_blank())
      plot_i <- plot_i + ggtitle(as.character(gene))
    }
    if (return.plots) {
      return(plot_i)
    } else {
      print(plot_i)
    }
  }  else {
    paste0("Object of class ", class(object)[[1]], " not recognized.")
  }
  
}


# note that setting point size to -1 seems to make it disappear 
# (trick that I can't find in documentation, but got from Seurat)
# also tried Pethukov's ggrastr, but it caused too much distortion 
# also started taking way too long when playing with raster.width > 4

# still working on fixing legend for first page of plots 

plotFactors_new <- function(object, num.genes = 8, cells.highlight = NULL, plot.tsne = F,
                            return.plots = F, pt.size1 = 0.5, pt.size2 = 1, 
                            tsne.colors = c('lemonchiffon', 'red'), use.png = F,
                            width = 7, height = 6, dpi = 100) {
  k <- ncol(object@H.norm)
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  
  W <- t(object@W)
  rownames(W) <- colnames(object@scale.data[[1]])
  Hs_norm <- object@H.norm
  H_raw = do.call(rbind, object@H)
  plot_list = list()
  tsne_list = list()
  
  for (i in 1:k) {
    # creates array of 2 by 1 plots 
    par(mfrow = c(2, 1))
    top_genes.W <- rownames(W)[order(W[, i], decreasing = T)[1:num.genes]]
    top_genes.W.string <- paste0(top_genes.W, collapse = ", ")
    factor_textstring <- paste0("Factor", i)
    
    plot_title1 <- paste(factor_textstring, "\n", top_genes.W.string, "\n")
    
    h_df = data.frame(x = 1:nrow(Hs_norm), h_norm = Hs_norm[, i],
                      h_raw = H_raw[, i], dataset = object@cell.data$dataset,
                      highlight = FALSE)
    
    top <- ggplot(h_df, aes(x = x, y=h_raw, col = dataset)) + geom_point(size = pt.size1) + 
      labs(x = 'Cell', y = 'Raw H Score') + ggtitle(plot_title1) + theme(legend.position = 'none')
    bottom <- ggplot(h_df, aes(x = x, y=h_norm, col = dataset)) + geom_point(size = pt.size1) + 
      labs(x = 'Cell', y = 'H_norm Score') + 
      theme(legend.position = 'top',
            legend.title = element_blank()) +
      guides(colour = guide_legend(override.aes = list(size = 2)))
    
    if (!is.null(cells.highlight)) {
      h_df[cells.highlight, 'highlight'] = TRUE
      top <- top + geom_point(data = subset(h_df, highlight == TRUE),
                              aes(x, h_raw),
                              col = "black",
                              size = pt.size1)
      bottom <- bottom + geom_point(data = subset(h_df, highlight == TRUE),
                                    aes(x, h_norm),
                                    col = "black",
                                    size = pt.size1)
    }
    
    if (use.png){
      top = AugmentPlot(top, width = width, height = height, dpi = dpi)
      top = top + labs(x = 'Cell', y = 'Raw H Score')
      # save legend before losing it -- might have to recreate manually
      bottom_legend = get_legend(bottom) 
      bottom = AugmentPlot(bottom, width = width, height = height, dpi = dpi)
      bottom = bottom + labs(x = 'Cell', y = 'H_norm Score') 
    }
    
    full = plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 1))
    plot_list[[i]] = full
    if (plot.tsne) {
      tsne_df <- data.frame(Hs_norm[, i], object@tsne.coords)
      factorlab <- paste0("Factor", i)
      colnames(tsne_df) <- c(factorlab, "tSNE1", "tSNE2")
      p1 <- ggplot(tsne_df, aes_string(x = "tSNE1", y = "tSNE2", color = factorlab)) + 
        geom_point(size = pt.size2) +
        scale_color_gradient(low = tsne.colors[1], high = tsne.colors[2]) + ggtitle(label = paste('Factor', i)) + 
        theme(legend.position = 'none')
      if (use.png){
        p1 = AugmentPlot(p1, width = width, height = height, dpi = dpi)
      }
      tsne_list[[i]] = p1
    }
    setTxtProgressBar(pb, i)
  }
  if (return.plots) {
    return(list(plot_list, tsne_list))
  } else {
    for (i in 1:k) {
      print(plot_list[[i]])
      print(tsne_list[[i]])
    }
  }
}

library(png)
# AugmentPlot taken from Seurat v3

#' Augments ggplot2-based plot with a PNG image.
#'
#' Creates "vector-friendly" plots. Does this by saving a copy of the plot as a PNG file,
#' then adding the PNG image with \code{\link[ggplot2]{annotation_raster}} to a blank plot
#' of the same dimensions as \code{plot}. Please note: original legends and axes will be lost
#' during augmentation.
#'
#' @param plot A ggplot object
#' @param width,height Width and height of PNG version of plot
#' @param dpi Plot resolution
#'
#' @return A ggplot object
#'
#' @importFrom png readPNG
#' @importFrom ggplot2 ggplot_build ggsave ggplot aes_string geom_blank annotation_raster ggtitle
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot <- DimPlot(object = pbmc_small)
#' AugmentPlot(plot = plot)
#' }
#'
AugmentPlot <- function(plot, width = 10, height = 10, dpi = 100) {
  pbuild.params <- ggplot_build(plot = plot)$layout$panel_params[[1]]
  range.values <- c(
    pbuild.params$x.range,
    pbuild.params$y.range
  )
  xyparams <- GetXYAesthetics(
    plot = plot,
    geom = class(x = plot$layers[[1]]$geom)[1]
  )
  title <- plot$labels$title
  tmpfile <- tempfile(fileext = '.png')
  ggsave(
    filename = tmpfile,
    plot = plot + NoLegend() + NoAxes() + theme(plot.title = element_blank()),
    width = width,
    height = height,
    dpi = dpi
  )
  img <- readPNG(source = tmpfile)
  file.remove(tmpfile)
  blank <- ggplot(
    data = plot$data,
    mapping = aes_string(x = xyparams$x, y = xyparams$y)
  ) + geom_blank()
  blank <- blank + plot$theme + ggtitle(label = title)
  blank <- blank + annotation_raster(
    raster = img,
    xmin = range.values[1],
    xmax = range.values[2],
    ymin = range.values[3],
    ymax = range.values[4]
  )
  return(blank)
}

# other functions which are necessary -- GetXYAesthetics, NoLegend, NoAxes

# Get X and Y aesthetics from a plot for a certain geom
#
# @param plot A ggplot2 object
# @param geom Geom class to filter to
# @param plot.first Use plot-wide X/Y aesthetics before geom-specific aesthetics
#
# @return A named list with values 'x' for the name of the x aesthetic and 'y' for the y aesthetic
#
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list('x' = x, 'y' = y))
}

NoLegend <- function(...) {
  no.legend.theme <- theme(
    # Remove the legend
    legend.position = 'none',
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(no.legend.theme)
}

#' @inheritParams SeuratTheme
#' @param keep.text Keep axis text
#' @param keep.ticks Keep axis ticks
#'
#' @importFrom ggplot2 theme element_blank
#' @export
#'
#' @rdname SeuratTheme
#' @aliases NoAxes
#'
#' @examples
#' # Generate a plot with no axes
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoAxes()
#'
NoAxes <- function(..., keep.text = FALSE, keep.ticks = FALSE) {
  blank <- element_blank()
  no.axes.theme <- theme(
    # Remove the axis elements
    axis.line.x = blank,
    axis.line.y = blank,
    # Validate the theme
    validate = TRUE,
    ...
  )
  if (!keep.text) {
    no.axes.theme <- no.axes.theme + theme(
      axis.text.x = blank,
      axis.text.y = blank,
      axis.title.x = blank,
      axis.title.y = blank,
      validate = TRUE,
      ...
    )
  }
  if (!keep.ticks){
    no.axes.theme <- no.axes.theme + theme(
      axis.ticks.x = blank,
      axis.ticks.y = blank,
      validate = TRUE,
      ...
    )
  }
  return(no.axes.theme)
}

# also some hadley wickham functions 
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}


calculate_celltype_associations <- function(ctd,gwas_sumstats_path,analysis_name="MainRun",upstream_kb=10,downstream_kb=1.5,genome_ref_path,specificity_species="mouse",genesOutCOND=NA,EnrichmentMode="Linear"){
  # Check EnrichmentMode has correct values
  if(!EnrichmentMode %in% c("Linear","Top 10%")){stop("EnrichmentMode argument must be set to either 'Linear' or 'Top 10%")}
  
  gwas_sumstats_path = path.expand(gwas_sumstats_path)
  magmaPaths = get.magma.paths(gwas_sumstats_path,upstream_kb,downstream_kb)
  
  # Check for errors in arguments
  check_inputs_to_magma_celltype_analysis(ctd,gwas_sumstats_path,analysis_name,upstream_kb,downstream_kb,genome_ref_path)
  
  output = list()
  for(annotLevel in 1:length(ctd)){
    sumstatsPrefix2 = sprintf("%s.level%s",magmaPaths$filePathPrefix,annotLevel)
    
    if(EnrichmentMode=="Linear"){
      # First match quantiles to the genes in the genes.out file... then write as the genesCovar file (the input to MAGMA)
      geneCovarFile = create_gene_covar_file(genesOutFile = sprintf("%s.genes.out",magmaPaths$filePathPrefix),ctd,annotLevel,specificity_species=specificity_species,genesOutCOND)
      
      if(is.na(genesOutCOND)){
        #magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --gene-covar '%s' onesided --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,sumstatsPrefix2,analysis_name)
        magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --gene-covar '%s' --model direction=pos --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,sumstatsPrefix2,analysis_name)
      }else{
        magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --gene-covar '%s' --model direction=pos  condition='ZSTAT' --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,sumstatsPrefix2,analysis_name)
      }
    }else if(EnrichmentMode=="Top 10%"){
      # First match quantiles to the genes in the genes.out file... then write as the genesCovar file (the input to MAGMA)
      geneCovarFile = create_top10percent_genesets_file(genesOutFile = sprintf("%s.genes.out",magmaPaths$filePathPrefix),ctd,annotLevel,specificity_species=specificity_species)
      
      if(is.na(genesOutCOND)){
        magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --set-annot '%s' --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,sumstatsPrefix2,analysis_name)
      }else{
        geneCovarFile2 = create_gene_covar_file(genesOutFile = sprintf("%s.genes.out",magmaPaths$filePathPrefix),ctd,annotLevel,specificity_species=specificity_species,genesOutCOND)
        magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --set-annot '%s' twosided --gene-covar '%s' condition-only='ZSTAT' --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,geneCovarFile2,sumstatsPrefix2,analysis_name)
      }
    }
    print(magma_cmd)
    system(magma_cmd)
    
    # Prepare output list
    tmp = list()
    tmp$geneCovarFile = geneCovarFile
    #if(EnrichmentMode=="Linear"){
    path = sprintf("%s.%s.gsa.out",sumstatsPrefix2,analysis_name)
    print(path)
    #}else if(EnrichmentMode=="Top 10%"){
    #    path = sprintf("%s.%s.sets.out",sumstatsPrefix2,analysis_name)
    #}
    tmp$results = load.magma.results.file(path,annotLevel,ctd,genesOutCOND=genesOutCOND,EnrichmentMode=EnrichmentMode,ControlForCT="BASELINE")
    output[[length(output)+1]] = tmp
  }
  
  # Calculate total number of tests performed
  totalTests = 0
  for(annotLevel in 1:length(output)){
    totalTests = totalTests + dim(output[[annotLevel]]$results)[1]
  }
  output$total_baseline_tests_performed = totalTests
  
  output$gwas_sumstats_path = gwas_sumstats_path
  output$analysis_name = analysis_name
  output$upstream_kb = upstream_kb
  output$downstream_kb = downstream_kb
  output$genome_ref_path = genome_ref_path
  
  return(output)
}


aggregate.over.celltypes <- function(rowOfMeans,celltypes,func="mean"){
  if(func=="mean"){
    exp_out = as.matrix(data.frame(aggregate(rowOfMeans,by=list(celltypes),FUN=mean)))
  }else if(func=="median"){
    exp_out = as.matrix(data.frame(aggregate(rowOfMeans,by=list(celltypes),FUN=median)))
  }
  rownames(exp_out) = exp_out[,"Group.1"]
  exp_out = exp_out[,2]
  exp_out2 = as.numeric(exp_out)
  names(exp_out2) = names(exp_out)
  return(exp_out2)
}
calculate.meanexp.for.level <- function(ctd_oneLevel,expMatrix){
  if(dim(expMatrix)[2]==length(unique(ctd_oneLevel$annot))){
    print(dim(expMatrix)[2])
    print(length(ctd_oneLevel$annot))
    if(sum(!colnames(expMatrix)==ctd_oneLevel$annot)!=0){
      stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
    }
    ctd_oneLevel$mean_exp = expMatrix
  }else{
    mean_exp = apply(expMatrix,1,aggregate.over.celltypes,ctd_oneLevel$annot)
    ctd_oneLevel$mean_exp = t(mean_exp)
  }
  return(ctd_oneLevel)
}
calculate.medianexp.for.level <- function(ctd_oneLevel,expMatrix){
  if(dim(expMatrix)[2]==length(unique(ctd_oneLevel$annot))){
    print(dim(expMatrix)[2])
    print(length(ctd_oneLevel$annot))
    if(sum(!colnames(expMatrix)==ctd_oneLevel$annot)!=0){
      stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
    }
    ctd_oneLevel$median_exp = expMatrix
  }else{
    median_exp = apply(expMatrix,1,aggregate.over.celltypes,ctd_oneLevel$annot,func="median")
    ctd_oneLevel$median_exp = t(median_exp)
  }
  return(ctd_oneLevel)
}
calculate.specificity.for.level <- function(ctd_oneLevel){
  normalised_meanExp = t(t(ctd_oneLevel$mean_exp)*(1/colSums(ctd_oneLevel$mean_exp)))
  normalised_medianExp = t(t(ctd_oneLevel$median_exp)*(1/colSums(ctd_oneLevel$mean_exp)))
  ctd_oneLevel$specificity = normalised_meanExp/(apply(normalised_meanExp,1,sum)+0.000000000001)
  ctd_oneLevel$median_specificity = normalised_medianExp/(apply(normalised_meanExp,1,sum)+0.000000000001)
  return(ctd_oneLevel)
}


annotate_mat = function(exp_mat, gset_10x, DAC)
{
  # Intersect dataset with 10x gset
  all_g = data.frame(row.names = rownames(exp_mat), ensembl_id = rep(0,nrow(exp_mat)),
                     MAG = rep(0,nrow(exp_mat)),
                     MIG = rep(0,nrow(exp_mat)))
  # Flag genes with (MAG) and without (MIG) ensembl id 
  all_g[intersect(rownames(all_g),rownames(gset_10x)),'MAG'] = 1
  all_g[intersect(rownames(all_g),rownames(gset_10x)),'ensembl_id'] = gset_10x[as.character(intersect(rownames(all_g),rownames(gset_10x))),'ensembl_id']
  all_g[setdiff(rownames(all_g),rownames(gset_10x)),'MIG'] = 1
  # look through http://mygene.info API to recover ensembl id from MIG
  if(length(rownames(all_g)[which(all_g$MIG == 1)]) > 0){ 
    xl = queryMany(rownames(all_g)[which(all_g$MIG == 1)], scopes="symbol", fields="ensembl.gene", species="mouse")
    matched = as.data.frame(xl[which(!is.na(xl$ensembl.gene)),])
    # Flag MAG recovered genes 
    all_g[matched$query,'ensembl_id'] = matched$ensembl.gene
    all_g[matched$query,'MIG'] = 0
    all_g[matched$query,'MAG'] = 1
    not_found = all_g[which(all_g$MIG == 1),]
  }
  MAG = all_g[which(all_g$MAG==1),]
  # Flag dupes and print the counts for dupes 
  # typically they differ and we throw them out 
  dupes = which(duplicated(MAG$ensembl_id))
  all_g$dupes = 0
  if(length(dupes)>0){
    dupes = MAG[dupes,'ensembl_id']
    gname_dupes = rownames(MAG[which(MAG$ensembl_id %in% dupes),])
    print(t(apply(exp_mat[gname_dupes,],1,sum)))
    all_g[gname_dupes,'dupes'] = 1
  }
  return(all_g)
}


# Spit out marker gene list to excel sheets
WriteExcelOut <- function(markerlist,output){
  a1<-split(markerlist, markerlist$group)
  if (is.numeric(unique(markerlist$group))){
    names(a1) <- as.vector(unlist(lapply(unique(markerlist$group),function(x){paste('Cluster',x,sep = "_" )})))
  }
  else
    names(a1) <- as.vector(unlist(lapply(unique(markerlist$group),function(x){x})))
  
  to.pad <- max(sapply(a1, nrow)) - sapply(a1,nrow)
  for (i in 1:length(a1)){
    a1[[i]][as.numeric(nrow(a1[[i]])+to.pad[i]),] <- NA
  }
  write.xlsx(a1, file = output)
}



## Convert seurat to monocle3
seuratToMonocle3 = function(object) {
  object_genes = data.frame(gene_short_name = row.names(object@raw.data))
  row.names(object_genes)= row.names(object@raw.data)
  
  phendat <- AnnotatedDataFrame(object@meta.data)
  #new_cds <- newCellDataSet(object@raw.data,
  #                             phenoData = phendat,
  #                             featureData = object_genes,
  #                             expressionFamily=negbinomial.size())
  
  new_cds <- new_cell_data_set(object@raw.data, 
                               cell_metadata = object@meta.data, 
                               gene_metadata = object_genes)
  return(new_cds)
}


## Makes ratio plot with median bars when you input a liger object
MakeRatioPlot <- function(object){
  t5 <- table(object@cell.data$dataset, object@clusters)
  t5 <- t5/rowSums(t5)
  t5 <- as.data.frame(t5)
  t5$status <- as.vector(unlist(lapply(t5$Var1, function(x){return(sample.sheet[sample.sheet$Code1 == x,]$Status)})))
  median.df <- aggregate(.~status+Var2,t5, FUN = median)
  median.ctrl <- median.df[grep('Ctrl',median.df$status),]
  
  t5$ctrlfreq <- rep(0,dim(t5)[1])
  ctrl.freq <- as.vector(unlist(lapply(t5$Var2, function(x){
    idx <- which(median.ctrl$Var2 == x)
    return(median.ctrl[idx,]$Freq)
  })))
  t5$ctrlfreq <- ctrl.freq
  t5$foldchange <- t5$Freq/t5$ctrlfreq
  t5$status <- factor(t5$status, levels= c('AbetaTau','Abeta','Ctrl'))
  
  p<-ggplot(t5,aes(x=Var2, y=foldchange, fill = status)) +
    geom_point(aes(x=Var2, y=foldchange,color = status), position=position_dodge(width=0.85),show.legend = F,size =2.5)+
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.75,position = position_dodge(width=0.85), color = 'black') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 8)) + geom_boxplot() + coord_flip() + ylim(0,15)
  return(p)
}


louvainCluster <- function(object, resolution = 1.0, k = 20, prune = 1 / 15, eps = 0.1, nRandomStarts = 10,
                           nIterations = 100, random.seed = 1) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  output_path <- paste0(getwd(), '/edge_test.txt')
  knn <- RANN::nn2(object@H.norm, k = k, eps = eps)
  snn <- Seurat:::ComputeSNN(knn$nn.idx, prune = prune)
  Seurat:::WriteEdgeFile(snn, output_path, display_progress = T)
  clusts <- Seurat:::RunModularityClusteringCpp(snn,
                                                modularityFunction = 1, resolution = resolution, nRandomStarts = nRandomStarts,
                                                nIterations = nIterations, algorithm = 1, randomSeed = random.seed, printOutput = T,
                                                edgefilename = output_path
  )
  object@clusters <- as.factor(clusts)
  names(object@clusters) <- unlist(lapply(object@H,rownames))
  unlink(output_path)
  return(object)
}

MakeFactorPlot <- function(object, path){
  pdf(path)
  plotFactors_new(object, num.genes = 8, plot.tsne = T)
  dev.off()
}




plotbyMetaVar <- function(object, metavar = NULL, title = NULL, pt.size = 0.3,
                          text.size = 3, do.shuffle = T, rand.seed = 1,
                          axis.labels = NULL, do.legend = T, legend.size = 5,
                          return.plots = F, downsample = F) {
  tsne_df <- data.frame(object@tsne.coords)
  if (downsample == T){
    idx.length <- min(table(object@cell.data[[metavar]]))
    cells.use2 <- as.vector(unlist(lapply(unique(object@cell.data[[metavar]]),function(x){
      cells.use <- sample(rownames(tsne_df)[which(object@cell.data[[metavar]] == x)],size = idx.length)
    })))
  }
  else
    cells.use2 <- rownames(tsne_df)
  tsne_df <- tsne_df[cells.use2,]
  colnames(tsne_df) <- c("tsne1", "tsne2")
  
  tsne_df$MetaVar <- object@cell.data[cells.use2,][[metavar]]
  if (do.shuffle) {
    set.seed(rand.seed)
    idx <- sample(1:nrow(tsne_df))
    tsne_df <- tsne_df[idx, ]
  }
  p1 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = MetaVar)) +
    geom_point_rast(size = pt.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size)))
  
  centers <- tsne_df %>% group_by(MetaVar) %>% summarize(
    tsne1 = median(x = tsne1),
    tsne2 = median(x = tsne2)
  )
  
  if (!is.null(title)) {
    p1 <- p1 + ggtitle(title[1])
  }
  if (!is.null(axis.labels)) {
    p1 <- p1 + xlab(axis.labels[1]) + ylab(axis.labels[2])
  }
  if (!do.legend) {
    p1 <- p1 + theme(legend.position = "none")
  }
  if (return.plots) {
    return(list(p1))
  } else {
    print(p1)
  }
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}



plotByDatasetAndCluster2 <- function(object, clusters = NULL, title = NULL, pt.size = 0.3,
                                     text.size = 3, do.shuffle = T, rand.seed = 1,
                                     axis.labels = NULL, do.legend = T, legend.size = 5,
                                     return.plots = F,cells.use = NULL,cells.highlight = NULL) {
  tsne_df <- data.frame(object@tsne.coords)
  colnames(tsne_df) <- c("tsne1", "tsne2")
  tsne_df$Dataset <- unlist(lapply(1:length(object@raw.data),function(x){
    rep(names(object@raw.data[x]),ncol(object@raw.data[[x]]))
  }))
  c_names <- names(object@clusters)
  if (is.null(clusters)) {
    if (length(object@clusters) == 0) {
      clusters <- rep(1, nrow(object@tsne.coords))
      names(clusters) <- c_names <- rownames(object@tsne.coords)
    } else {
      clusters <- object@clusters
      c_names <- names(object@clusters)
    }
  }
  tsne_df$Cluster <- clusters[c_names]
  if (do.shuffle) {
    set.seed(rand.seed)
    idx <- sample(1:nrow(tsne_df))
    tsne_df <- tsne_df[idx, ]
  }
  
  if (!is.null(cells.use)){
    tsne_df <- tsne_df[match(cells.use,rownames(tsne_df)),]
  }
  
  tsne_df$Cluster <- droplevels(tsne_df$Cluster)
  
  p1 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = Dataset)) +
    geom_point(size = pt.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size)))
  
  centers <- tsne_df %>% group_by(Cluster) %>% summarize(
    tsne1 = median(x = tsne1),
    tsne2 = median(x = tsne2)
  )
  centers$Cluster <- levels(tsne_df$Cluster)
  p2 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = Cluster)) + geom_point(size = pt.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size))) + 
    ggrepel::geom_text_repel(data = centers, mapping = aes(label = Cluster),
                             colour = "black", size = text.size)
  
  if (!is.null(cells.highlight)){
    tsne_df$colorval <- rep('black',nrow(tsne_df))
    tsne_df[cells.highlight,]$colorval <- rep('red',length(cells.highlight))
    p2 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = colorval)) + geom_point(size = pt.size) +
      guides(color = guide_legend(override.aes = list(size = legend.size)))
  }
  
  
  if (!is.null(title)) {
    p1 <- p1 + ggtitle(title[1])
    p2 <- p2 + ggtitle(title[2])
  }
  if (!is.null(axis.labels)) {
    p1 <- p1 + xlab(axis.labels[1]) + ylab(axis.labels[2])
    p2 <- p2 + xlab(axis.labels[1]) + ylab(axis.labels[2])
  }
  if (!do.legend) {
    p1 <- p1 + theme(legend.position = "none")
    p2 <- p2 + theme(legend.position = "none")
  }
  if (return.plots) {
    return(list(p1,p2))
  } else {
    print(p1)
    print(p2)
  }
}




pbdac_df <- function(object, clusters = NULL, title = NULL, pt.size = 0.3,
                     text.size = 3, do.shuffle = T, rand.seed = 1,
                     axis.labels = NULL, do.legend = T, legend.size = 5,
                     return.plots = F,cells.use = NULL,cells.highlight = NULL) {
  tsne_df <- object[,c('tsne1','tsne2','dataset','Cluster')]
  if (do.shuffle) {
    set.seed(rand.seed)
    idx <- sample(1:nrow(tsne_df))
    tsne_df <- tsne_df[idx, ]
  }
  
  if (!is.null(cells.use)){
    tsne_df <- tsne_df[match(cells.use,rownames(tsne_df)),]
  }
  
  p1 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = dataset)) +
    geom_point(size = pt.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size)))
  
  centers <- tsne_df %>% dplyr::group_by(Cluster) %>% summarize(
    tsne1 = median(x = tsne1),
    tsne2 = median(x = tsne2)
  )
  
  centers$Cluster <- levels(tsne_df$Cluster)
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(tsne_df$Cluster)))
  print(mycolors)
  p2 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = Cluster)) + geom_point(size = pt.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size))) + 
    ggrepel::geom_label_repel(data = centers, mapping = aes(label = Cluster),
                              colour = "black", size = text.size,fill = alpha(c("white"),0.5)) +
    scale_fill_manual(values = mycolors)
  
  if (!is.null(cells.highlight)){
    tsne_df$colorval <- rep('black',nrow(tsne_df))
    tsne_df[cells.highlight,]$colorval <- rep('red',length(cells.highlight))
    p2 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = colorval)) + geom_point(size = pt.size) +
      guides(color = guide_legend(override.aes = list(size = legend.size)))+ 
      scale_color_manual(values = mycolors)
  }
  
  if (!is.null(title)) {
    p1 <- p1 + ggtitle(title[1])
    p2 <- p2 + ggtitle(title[2])
  }
  if (!is.null(axis.labels)) {
    p1 <- p1 + xlab(axis.labels[1]) + ylab(axis.labels[2])
    p2 <- p2 + xlab(axis.labels[1]) + ylab(axis.labels[2])
  }
  if (!do.legend) {
    p1 <- p1 + theme(legend.position = "none")
    p2 <- p2 + theme(legend.position = "none")
  }
  if (return.plots) {
    return(list(p1,p2))
  } else {
    print(p1)
    print(p2)
  }
}





plotGene2 <- function(object, gene, use.raw = F, use.scaled = F, scale.by = 'dataset', 
                      log2scale = NULL, methylation.indices = NULL, plot.by = 'dataset', 
                      set.dr.lims = F, pt.size = 0.1, min.clip = NULL, max.clip = NULL, 
                      clip.absolute = F, points.only = F, option = 'plasma', cols.use = NULL, 
                      zero.color = '#F5F5F5', axis.labels = NULL, do.legend = T, return.plots = F) {
  if ((plot.by != scale.by) & (use.scaled)) {
    warning("Provided values for plot.by and scale.by do not match; results may not be very
            interpretable.")
  }
  if (use.raw) {
    if (is.null(log2scale)) {
      log2scale <- FALSE
    }
    # drop only outer level names
    gene_vals <- getGeneValues(object@raw.data, gene, log2scale = log2scale)
  } else {
    if (is.null(log2scale)) {
      log2scale <- TRUE
    }
    # rescale in case requested gene not highly variable
    if (use.scaled) {
      # check for feature 
      if (!(scale.by %in% colnames(object@cell.data)) & scale.by != 'none') {
        stop("Please select existing feature in cell.data to scale.by, or add it before calling.")
      }
      gene_vals <- getGeneValues(object@norm.data, gene)
      cellnames <- names(gene_vals)
      # set up dataframe with groups
      gene_df <- data.frame(gene = gene_vals)
      if (scale.by == 'none') {
        gene_df[['scaleby']] = 'none'
      } else {
        gene_df[['scaleby']] = factor(object@cell.data[[scale.by]])
      }
      gene_df1 <- gene_df %>%
        group_by(scaleby) %>%
        # scale by selected feature
        mutate_at(vars(-group_cols()), function(x) { scale(x, center = F)})
      gene_vals <- gene_df1$gene
      names(gene_vals) <- cellnames
      if (log2scale) {
        gene_vals <- log2(10000 * gene_vals + 1)
      }
    } else {
      # using normalized data
      # indicate methylation indices here 
      gene_vals <- getGeneValues(object@norm.data, gene, methylation.indices = methylation.indices,
                                 log2scale = log2scale)
    }
  }
  gene_vals[gene_vals == 0] <- NA
  dr_df <- data.frame(object@tsne.coords)
  rownames(dr_df) <- rownames(object@cell.data)
  dr_df$gene <- as.numeric(gene_vals[rownames(dr_df)])
  colnames(dr_df) <- c("dr1", "dr2", "gene")
  # get dr limits for later
  lim1 <- c(min(dr_df$dr1), max(dr_df$dr1))
  lim2 <- c(min(dr_df$dr2), max(dr_df$dr2))
  
  if (plot.by != 'none') {
    if (!(plot.by %in% colnames(object@cell.data))) {
      stop("Please select existing feature in cell.data to plot.by, or add it before calling.")
    }
    dr_df$plotby <- factor(object@cell.data[[plot.by]])
  } else {
    dr_df$plotby <- factor("none")
  }
  # expand clip values if only single provided
  num_levels <- length(levels(dr_df$plotby))
  if (length(min.clip) == 1) {
    min.clip <- rep(min.clip, num_levels)
    names(min.clip) <- levels(dr_df$plotby)
  }
  if (length(max.clip) == 1) {
    max.clip <- rep(max.clip, num_levels)
    names(max.clip) <- levels(dr_df$plotby)
  }
  if (!is.null(min.clip) & is.null(names(min.clip))) {
    if (num_levels > 1) {
      message("Adding names to min.clip according to levels in plot.by group; order may not be 
              preserved as intended if multiple clip values passed in. Pass in named vector to 
              prevent this.")
    }
    names(min.clip) <- levels(dr_df$plotby)
  }
  if (!is.null(max.clip) & is.null(names(max.clip))) {
    if (num_levels > 1) {
      message("Adding names to max.clip according to levels in plot.by group; order may not be 
              preserved as intended if multiple clip values passed in. Pass in named vector to 
              prevent this.")
    }
    names(max.clip) <- levels(dr_df$plotby)
  }
  p_list <- list()
  for (sub_df in split(dr_df, f = dr_df$plotby)) {
    # maybe do quantile cutoff here
    group_name <- as.character(sub_df$plotby[1])
    if (!clip.absolute) {
      max_v <- quantile(sub_df$gene, probs = max.clip[group_name], na.rm = T)
      min_v <- quantile(sub_df$gene, probs = min.clip[group_name], na.rm = T)
    } else {
      max_v <- max.clip[group_name]
      min_v <- min.clip[group_name]
    }
    sub_df$gene[sub_df$gene < min_v & !is.na(sub_df$gene)] <- min_v
    sub_df$gene[sub_df$gene > max_v & !is.na(sub_df$gene)] <- max_v
    
    ggp <- ggplot(sub_df, aes(x = dr1, y = dr2, color = gene)) + geom_point(size = pt.size) +
      labs(col = gene)
    
    if (!is.null(cols.use)) {
      ggp <- ggp + scale_color_gradientn(colors = cols.use,
                                         na.value = zero.color)
    } else {
      ggp <- ggp + scale_color_viridis_c(option = option,
                                         direction = -1,
                                         na.value = zero.color)
    }
    if (set.dr.lims) {
      ggp <- ggp + xlim(lim1) + ylim(lim2)
    }
    
    if (plot.by != 'none') {
      base <- as.character(sub_df$plotby[1])
    } else {
      base <- ""
    }
    ggp <- ggp + ggtitle(base)
    
    if (!is.null(axis.labels)) {
      ggp <- ggp + xlab(axis.labels[1]) + ylab(axis.labels[2])
    }
    if (!do.legend) {
      ggp <- ggp + theme(legend.position = "none")
    }
    if (points.only) {
      ggp <- ggp + theme(
        axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "none",
        panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank()
      ) + ggtitle(as.character(gene))
    }
    p_list[[as.character(sub_df$plotby[1])]] <- ggp
  }
  if (plot.by == 'dataset') {
    p_list <- p_list[names(object@raw.data)]
  }
  
  if (return.plots){
    if (length(p_list) == 1) {
      return(p_list[[1]])
    } else {
      return(p_list)
    }
  } else {
    for (plot in p_list) {
      print(plot)
    }
  }
}



loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}

source_https <- function(u, unlink.tmp.certs = FALSE) {
  # load package
  require(RCurl)
  
  # read script lines from website using a security certificate
  if(!file.exists("cacert.pem")) download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.pem")
  script <- getURL(u, followlocation = TRUE, cainfo = "cacert.pem")
  if(unlink.tmp.certs) unlink("cacert.pem")
  
  # parase lines and evealuate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}


Stacked_VlnPlot <- function(liger.object, features, 
                            x_lab_rotate = FALSE, plot_margin = unit(c(-0.5, 0, -0.5, 0), "cm"), ...) {
  # Check whether features are present in object
  possible_features <- rownames(liger.object@raw.data)
  if (any(!features %in% possible_features)) {
    bad_features <- features[!features %in% possible_features]
    features <- features[features %in% possible_features]
    if(length(x = features) == 0) {
      stop("No requested features found.")
    }
    warning("The following features were omitted as they were not found",
            ": ", paste(bad_features, collapse = ", "))
  }
  ## remove the x-axis text and tick
  ## plot_margin to adjust the white space between each plot.
  Modify_VlnPlot<- function(liger.object, features, pt.size = 0, 
                            plot_margin = unit(c(-0.5, 0, -0.5, 0), "cm"), ...) {
    VlnPlot(liger.object, features = features, pt.size = pt.size, ...)  +
      xlab("") +
      ylab(features) +
      ggtitle("") +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = rel(1), angle = 0),
            axis.text.y = element_text(size = rel(1)),
            plot.margin = plot_margin)
  }
  
  ## extract the max value of the y axis
  extract_max <- function(p){
    ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
  }
  
  plot_list <- map(features, function(x) Modify_VlnPlot(seurat_object = seurat_object, 
                                                        features = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  # Add ability to rotate the X axis labels to the function call
  if (x_lab_rotate) {
    plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
      theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks.x = element_line())
  }
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value
  ymaxs <- map_dbl(plot_list, extract_max)
  plot_list <- suppressMessages(map2(plot_list, ymaxs, function(x,y) x +
                                       scale_y_continuous(breaks = c(y)) +
                                       expand_limits(y = y)))
  
  wrap_plots(plotlist = plot_list, ncol = 1)
}

