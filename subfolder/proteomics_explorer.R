# Load Libraries ----
library(shiny)
library(shinyTree)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(readr)

# Load Your Data ----
protein_map <- read.delim2(file = 'sub_folder/UniProt2Reactome_All_Levels.txt', sep = '\t', header = FALSE)
protein_map <- protein_map[grepl("R-MMU", protein_map$V2), ]
colnames(protein_map) <- c("UniProt", "PathwayID", "PathwayBrowser", "PathwayName", "OrganismCode", "OrganismName")

pathway_list <- read.delim2(file = 'sub_folder/ReactomePathways.txt', sep = '\t', header = FALSE)
pathway_list <- pathway_list[grepl("R-MMU", pathway_list$V1), ]
colnames(pathway_list) <- c("PathwayID", "PathwayName", "OrganismName")

pathway_hierarchy <- read.delim2(file = 'sub_folder/ReactomePathwaysRelation.txt', sep = '\t', header = FALSE)
pathway_hierarchy <- pathway_hierarchy[grepl("R-MMU", pathway_hierarchy$V1) & grepl("R-MMU", pathway_hierarchy$V2), ]
colnames(pathway_hierarchy) <- c("parent", "child")

diff_expr <- read.delim2('sub_folder/combined_all.txt', sep = '\t', header = TRUE)
diff_expr$logFC <- as.numeric(diff_expr$logFC)
diff_expr$adjPval <- as.numeric(diff_expr$adjPval)
diff_expr$Regulation <- ifelse(diff_expr$logFC > 0, "Up", "Down")

all_tissues <- unique(diff_expr$Tissue)
all_times <- c("Ob", "STR", "MTR", "LTR")
all_regulations <- c("Up", "Down")

# Helper: Build Tree Structure ----
build_tree <- function(pathway_hierarchy, pathway_list) {
  
  id2name <- setNames(pathway_list$PathwayName, pathway_list$PathwayID)
  kids    <- split(pathway_hierarchy$child, pathway_hierarchy$parent)
  
  make_node <- function(id) {
    child_ids <- kids[[id]]
    
    if (is.null(child_ids)) {
      # ---- leaf ---------------------------------------------------
      structure(list(), stid = id)
    } else {
      # ---- internal node -----------------------------------------
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

# Helper: Find Selected Pathway ID ----
find_pathway_id <- function(selected_name, pathway_list) {
  pathway_list$PathwayID[pathway_list$PathwayName == selected_name][1]
}

# Helper: Get Bubble Plot Data ----
get_bubble_data <- function(selected_pathway_id, protein_map, diff_expr, padj_cutoff, logfc_cutoff) {
  
  if (is.null(selected_pathway_id)) return(tibble())
  
  proteins_in_pathway <- protein_map$UniProt[
    protein_map$PathwayID == selected_pathway_id]
  
  if (length(proteins_in_pathway) == 0) return(tibble())
  
  bubble_df <- diff_expr %>%
    filter(Proteins %in% proteins_in_pathway,
           adjPval < padj_cutoff,
           abs(logFC) > logfc_cutoff) %>%
    mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
    group_by(Tissue, Time, Regulation) %>%
    summarise(
      Count = n(),
      Genes = {
        genes <- unique(PG.Genes)
        # Insert <br> every 5 genes
        paste0(sapply(seq(1, length(genes), by = 10), function(i) {
          paste(genes[i:min(i+9, length(genes))], collapse = ", ")
        }), collapse = "<br>")
      },
      .groups = "drop"
    )
  
  bubble_df <- bubble_df %>%
    complete(
      Tissue = all_tissues,
      Time = all_times,
      Regulation = all_regulations,
      fill = list(Count = 0,
                  Genes = "")
    )
  
  return(bubble_df)
}

# Helper: Make Bubble Plot ----
make_bubble_plot <- function(bubble_df) {
  if (nrow(bubble_df) == 0) return(NULL)
  
  bubble_df$Tissue <- factor(bubble_df$Tissue, levels = rev(sort(unique(bubble_df$Tissue))))
  bubble_df$Time <- factor(bubble_df$Time, levels = c("Ob", "STR", "MTR", "LTR"))
  bubble_df$Regulation <- factor(bubble_df$Regulation, levels = c("Up", "Down"))
  
  plot_df <- dplyr::filter(bubble_df, Count > 0)
  
  p <- suppressWarnings(ggplot(bubble_df, aes(x = Time, y = Tissue)) +
    geom_point(data = plot_df,
               aes(size = Count,
                   color = Regulation,
                   group = Regulation,
                   text = paste0(
                     "<b>Tissue: </b>", Tissue, "\n",
                     "<b>Time: </b>",   Time,   "\n",
                     "<b>Regulation: </b>", Regulation, "\n",
                     "<b>Count: </b>",  Count, "\n",
                     "<b>Genes: </b>",  Genes, "\n")),
               position = position_dodge(width = 0.7),
               alpha = 0.9, shape = 16) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    scale_color_manual(values = c("Up"   = "#E64B35FF",
                                  "Down" = "#4DBBD5FF")) +
    scale_size(range = c(1, 10), breaks = pretty) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text  = element_text(size = 10),
          legend.position = "none") +
    labs(size = "Protein Count"))
  
  return(ggplotly(p, tooltip = "text") %>% 
           layout(hoverlabel = list(align = "left")) %>%
           config(displayModeBar = FALSE))
}

# UI ----
ui <- fluidPage(
  div(
    style = "background-color: #333333; padding: 20px;",
    h1("Reactome Pathway-Based Proteomics Explorer for Obesity and Its Regression", 
       style = "color: white; text-align: center; margin: 0;")
  ),
  sidebarLayout(
    sidebarPanel(
      style = "border-radius: 0px 0px 10px 10px; overflow-x: auto; white-space: nowrap; font-size: 85%;",
      numericInput("padj_cutoff", "BH-Adjusted P-Value Cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("logfc_cutoff", "Absolute logFC Cutoff:", value = 0, min = 0, step = 0.1),
      p(strong("Select a pathway to plot:")),
      div(
        style = "min-width: 300px; max-width: 100%; overflow-x: auto;",
        shinyTree("pathway_tree", search = FALSE, theme = "proton")
      )
    ),
    mainPanel(
      uiOutput("selected_pathway_name"),
      plotlyOutput("bubble_plot", height = "700px")
    )
  )
)

# Server ----
server <- function(input, output, session) {
  
  tree_data <- reactive({
    build_tree(pathway_hierarchy, pathway_list)
  })
  
  output$pathway_tree <- renderTree({
    tree_data()
  })
  
  selected_pathway_id <- reactive({
    sel <- get_selected(input$pathway_tree, format = "classid")
    if (is.null(sel) || length(sel) == 0) return(NULL)
    
    attr(sel[[1]], "stid")
  })
  
  bubble_data <- reactive({
    pid <- selected_pathway_id()
    if (is.null(pid)) return(NULL)
    
    bd <- get_bubble_data(
      pid,
      protein_map,
      diff_expr,
      input$padj_cutoff,
      input$logfc_cutoff
    )
    
    max_count <- if (all(is.na(bd$Count))) 0 else max(bd$Count, na.rm = TRUE)
    
    list(bubble = bd, max_count = max_count)
  })
  
  output$selected_pathway_name <- renderUI({
    pid <- selected_pathway_id()
    bd_all <- bubble_data()
    
    if (is.null(pid) || is.null(bd_all)) {
      div(
        h3("No pathway selected", style = "text-align: center;")
      )
    } else {
      pname <- pathway_list$PathwayName[pathway_list$PathwayID == pid][1]
      max_count <- bd_all$max_count
      
      div(
        h3(pname, style = "text-align: center;"),
        div(pid, style = "text-align: center; color: grey; font-size: 14px; margin-top: -10px;"),
        div(paste0("Largest proteinset: ", max_count), 
            style = "text-align: center; color: #555555; font-size: 14px; margin-top: 5px;")
      )
    }
  })
  
  output$bubble_plot <- renderPlotly({
    bd_all <- bubble_data()
    if (is.null(bd_all)) return(NULL)
    bd <- bd_all$bubble
    make_bubble_plot(bd)
  })
}
