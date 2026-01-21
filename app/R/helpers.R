empty_bubble_df <- function() {
  data.frame(
    Tissue = character(),
    Time = character(),
    Regulation = character(),
    Count = integer(),
    Genes = character(),
    stringsAsFactors = FALSE
  )
}

parse_gene_input <- function(text) {
  if (is.null(text) || !nzchar(text)) return(character())
  genes <- unlist(strsplit(text, "\\r?\\n"))
  genes <- trimws(genes)
  unique(genes[genes != ""])
}

build_tree <- function(pathway_hierarchy, pathway_list) {
  id2name <- setNames(pathway_list$PathwayName, pathway_list$PathwayID)
  kids <- split(pathway_hierarchy$child, pathway_hierarchy$parent)

  make_node <- function(id) {
    child_ids <- kids[[id]]

    if (is.null(child_ids)) {
      structure(list(), stid = id)
    } else {
      sorted_child_ids <- child_ids[order(id2name[child_ids])]
      sub_nodes <- lapply(sorted_child_ids, make_node)
      names(sub_nodes) <- id2name[sorted_child_ids]
      structure(sub_nodes, stid = id)
    }
  }

  roots <- setdiff(pathway_hierarchy$parent, pathway_hierarchy$child)
  sorted_roots <- roots[order(id2name[roots])]
  forest <- lapply(sorted_roots, make_node)
  names(forest) <- id2name[sorted_roots]

  forest
}

get_bubble_data <- function(selected_pathway_id, protein_map, diff_expr,
                            padj_cutoff, logfc_cutoff,
                            all_tissues, all_times, all_regulations) {
  if (is.null(selected_pathway_id)) return(empty_bubble_df())

  proteins_in_pathway <- protein_map$UniProt[protein_map$PathwayID == selected_pathway_id]
  if (length(proteins_in_pathway) == 0) return(empty_bubble_df())

  bubble_df <- diff_expr %>%
    dplyr::filter(
      Proteins %in% proteins_in_pathway,
      adjPval < padj_cutoff,
      abs(logFC) > logfc_cutoff
    ) %>%
    dplyr::mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
    dplyr::group_by(Tissue, Time, Regulation) %>%
    dplyr::summarise(
      Count = dplyr::n(),
      Genes = {
        genes <- unique(PG.Genes)
        paste0(sapply(seq(1, length(genes), by = 10), function(i) {
          paste(genes[i:min(i + 9, length(genes))], collapse = ", ")
        }), collapse = "<br>")
      },
      .groups = "drop"
    )

  bubble_df <- bubble_df %>%
    tidyr::complete(
      Tissue = all_tissues,
      Time = all_times,
      Regulation = all_regulations,
      fill = list(Count = 0, Genes = "")
    )

  bubble_df
}

make_bubble_plot <- function(bubble_df) {
  if (nrow(bubble_df) == 0) return(NULL)

  bubble_df$Tissue <- factor(bubble_df$Tissue, levels = rev(sort(unique(bubble_df$Tissue))))
  bubble_df$Time <- factor(bubble_df$Time, levels = c("Ob", "STR", "MTR", "LTR"))
  bubble_df$Regulation <- factor(bubble_df$Regulation, levels = c("Up", "Down"))

  plot_df <- dplyr::filter(bubble_df, Count > 0)

  p <- suppressWarnings(
    ggplot2::ggplot(bubble_df, ggplot2::aes(x = Time, y = Tissue)) +
      ggplot2::geom_point(
        data = plot_df,
        ggplot2::aes(
          size = Count,
          color = Regulation,
          group = Regulation,
          text = paste0(
            "<b>Tissue: </b>", Tissue, "\n",
            "<b>Time: </b>", Time, "\n",
            "<b>Regulation: </b>", Regulation, "\n",
            "<b>Count: </b>", Count, "\n",
            "<b>Genes: </b>", Genes, "\n"
          )
        ),
        position = ggplot2::position_dodge(width = 0.7),
        alpha = 0.9,
        shape = 16
      ) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#2C7FB8")) +
      ggplot2::scale_size(range = c(1, 10), breaks = pretty) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = 10),
        legend.position = "none"
      ) +
      ggplot2::labs(size = "Protein Count")
  )

  plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(hoverlabel = list(align = "left")) %>%
    plotly::config(displayModeBar = FALSE)
}

make_gene_plot <- function(plot_df, tissue_colors) {
  if (nrow(plot_df) == 0) return(NULL)

  plot_df$Time <- factor(plot_df$Time, levels = c("Ob", "STR", "MTR", "LTR"))
  plot_df$Tissue <- factor(plot_df$Tissue, levels = names(tissue_colors))

  p <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    ggplot2::geom_line(
      ggplot2::aes(x = Time, y = logFC, group = Tissue, color = Tissue),
      alpha = 1,
      linewidth = 0.7
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "#6b7280", linewidth = 0.4),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 9),
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 10),
      strip.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0,0,0,0)
    ) +
    ggplot2::scale_color_manual(values = tissue_colors) +
    ggplot2::facet_wrap(~PG.Genes, ncol = 2) +
    ggplot2::labs(x = NULL, y = NULL)

  plotly_obj <- plotly::ggplotly(p, tooltip = c("x", "y", "colour"))
  time_levels <- levels(plot_df$Time)

  for (i in seq_along(plotly_obj$x$data)) {
    trace <- plotly_obj$x$data[[i]]
    tissue_name <- ""
    if (!is.null(trace$name)) {
      tissue_name <- gsub("^Tissue\\s*:?\\s*", "", trace$name)
    }
    if (!is.null(trace$x)) {
      if (is.numeric(trace$x)) {
        time_index <- pmax(1, pmin(length(time_levels), round(trace$x)))
        time_labels <- time_levels[time_index]
      } else {
        time_labels <- as.character(trace$x)
      }
      trace$customdata <- cbind(time_labels, rep(tissue_name, length(time_labels)))
    }
    trace$hovertemplate <- paste0(
      "Time: %{customdata[0]}<br>",
      "logFC: %{y:.2f}<br>",
      "Tissue: %{customdata[1]}<extra></extra>"
    )
    trace$mode <- "lines+markers"
    trace$marker <- list(size = 6, opacity = 0)
    trace$hoveron <- "points"
    plotly_obj$x$data[[i]] <- trace
  }

  axis_names <- grep("^xaxis", names(plotly_obj$x$layout), value = TRUE)
  for (axis_name in axis_names) {
    axis <- plotly_obj$x$layout[[axis_name]]
    axis$matches <- NULL
    axis$showticklabels <- TRUE
    axis$showline <- TRUE
    axis$linecolor <- "#cbd5e1"
    axis$linewidth <- 0.6
    axis$ticks <- "outside"
    axis$mirror <- FALSE
    axis$zeroline <- FALSE
    axis$automargin <- TRUE
    plotly_obj$x$layout[[axis_name]] <- axis
  }

  y_axis_names <- grep("^yaxis", names(plotly_obj$x$layout), value = TRUE)
  for (axis_name in y_axis_names) {
    axis <- plotly_obj$x$layout[[axis_name]]
    axis$matches <- NULL
    axis$showticklabels <- TRUE
    axis$showline <- TRUE
    axis$linecolor <- "#cbd5e1"
    axis$linewidth <- 0.6
    axis$ticks <- "outside"
    axis$mirror <- FALSE
    axis$zeroline <- FALSE
    axis$automargin <- TRUE
    plotly_obj$x$layout[[axis_name]] <- axis
  }

  plotly_obj %>%
    plotly::layout(margin = list(t = 50, b = 50, l = 40, r = 20)) %>%
    plotly::config(displayModeBar = FALSE)
}
