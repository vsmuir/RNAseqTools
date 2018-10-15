#' Plot the median expression value for a set of genes using a counts object and annotation dataframe
#'
#' @param counts Matrix or dataframe of counts or batch corrected residuals.  The dimension corresponding to genes must have dimnames, which will be matched to the gene identifiers in \code{gene_sets}. The dimension corresponding to library or sample IDs must have dimnames, which will be matched to the sample identifiers in \code{annotation}.
#' @param module The module or gene set, as a character vector
#' @param annotation Dataframe where libraries that should be included in plot are in rows and phenotypic/quality characteristics are in columns.
#' @param id Name of the column in \code{annotation} object corresponding to the library IDs used to name samples in the \code{counts} object. 
#' @param group Name of the column in \code{annotation} object corresponding to experimental groups that should be plotted.
#' @param timepoint If NULL, data is not expected to be shown over time (default). Otherwise, the name of the column in \code{annotation} that will contain the existing time information.
#' @param palette The discrete color palette to be used for groups.  Defaults to a 6-color scheme.
#' @param sig_comparison Argument to pass on to geom_signif.  This should be a list of length-2 vectors that correspond to the group names in a pairwise comparison.  These names must be present in the \code{group} column of the \code{annotation} object.
#' @param ... Optional parameters passed to \code{gene_set_median}.
#' @export
#' @usage \code{plot_medians(
#'   counts, module,
#'   annotation, id,
#'   group, timepoint = NULL,
#'   palette, sig_comparison = NULL,
#'   ...)}
#' @return A plot showing the median expression values separated by experimental group, with either one element for each library in \code{counts} (if time series data was not collected) or one element for each group at each time point.

plot_medians <- function(counts, module, annotation, id, group, 
                         timepoint = NULL, sig_comparison = NULL,
                         palette = c("#CDBE6A", "#2F5075", "#A45B37", "#873651", "#308491", "#391335"), ...){
  
  if (!is.data.frame(annotation)) {
    stop("Annotation object must be a data frame.")
  }
  
  if (!is.character(module)) {
    stop("Module object must be a character vector.")
  }  
  
  if (!require(dplyr) & !require(geneSetTools) & !require(ggbeeswarm) & !require(ggsignif)) {
    stop("Sorry, I'm missing some needed packages.\nPlease install the tidyverse, geneSetTools, ggbeeswarm, and ggsignif packages.\n")
  }
  
  
  med_mod_counts = geneSetTools::gene_set_median_count(module, counts, ...) %>%
    as_data_frame() %>%
    rownames_to_column(var = id) %>%
    dplyr::rename(median_count = value)
  
  plot_df = annotation %>% 
    dplyr::select(id, group, timepoint) %>%
    left_join(med_mod_counts)
  
  if(is.null(timepoint)){
    
    plot = ggplot2::ggplot(plot_df, aes(x = get(group), y = median_count, color = get(group))) +
      geom_boxplot(outlier.colour = NA, fill = NA) +
      ggbeeswarm::geom_beeswarm(shape = 1) +
      ggsignif::geom_signif(comparisons = sig_comparison, map_signif_level = T) +
      scale_color_manual(values = palette, guide=FALSE) +
      labs(x = "Group", y = "Median Module\nExpression")
    
  } else {
    
    plot_time_df = plot_df %>%
      group_by_(timepoint, group) %>%
      summarize(mean = mean(median_count), sem = sd(median_count)/sqrt(length(median_count)))
    
    mw_results = plot_df %>%
      subset(get(group) %in% unlist(sig_comparison)) %>%
      group_by_(timepoint) %>%
      summarize(p = wilcox.test(median_count ~ get(group))$p.value) %>%
      mutate(fdr = p.adjust(p))
    
    if(sum(mw_results$fdr < 0.05) < 1){
      cat("The groups you are interested in did not show significant differences at any timepoint.")
    } else {
      cat("The groups you are interested in are significantly different at times: ", unlist(mw_results[mw_results$fdr < 0.05,timepoint]))
    }
    
    plot = ggplot2::ggplot(plot_time_df, aes(x = as.numeric(get(timepoint)), y = mean, color = get(group))) +
      geom_line(size = 1.25) +
      geom_point() +
      geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), width = 2) +
      scale_color_manual(values = palette) +
      labs(x = "Time", y = "Median Module Expression\nAveraged by Group", color = "Group")
    
  }
  
  return(plot)
}
