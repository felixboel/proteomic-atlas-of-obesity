library(shiny)
library(shinyTree)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)

REVIEW_MODE <- Sys.getenv("REVIEW_MODE", "0") == "1"
if (REVIEW_MODE) {
  library(shinymanager)
}

source("R/config.R")
source("R/data_load.R")
source("R/helpers.R")

app_data <- load_processed_data()
protein_map <- app_data$protein_map
pathway_list <- app_data$pathway_list
pathway_hierarchy <- app_data$pathway_hierarchy
diff_expr <- app_data$diff_expr

diff_expr$logFC <- as.numeric(diff_expr$logFC)
diff_expr$adjPval <- as.numeric(diff_expr$adjPval)
diff_expr$Regulation <- ifelse(diff_expr$logFC > 0, "Up", "Down")

all_tissues <- sort(unique(diff_expr$Tissue))
all_times <- c("Ob", "STR", "MTR", "LTR")
all_regulations <- c("Up", "Down")
tissue_palette <- c(
  "Bone" = "#413425",
  "Brain" = "#3E7B38",
  "Duodenum" = "#3BC14A",
  "eWAT" = "#C2DF69",
  "Heart" = "#85FFC7",
  "iBAT" = "#52D1DC",
  "iWAT" = "#3A86FF",
  "Kidney" = "#9899A6",
  "Liver" = "#7F055F",
  "Muscle" = "#E07BE0",
  "Pancreas" = "#EE0077",
  "Serum" = "#B10E35",
  "Small_Intestine" = "#E3170A",
  "Spleen" = "#F16B0B",
  "Thymus" = "#FFBE0B"
)
tissue_colors <- setNames(rep("#9ca3af", length(all_tissues)), all_tissues)
palette_names <- intersect(names(tissue_palette), names(tissue_colors))
tissue_colors[palette_names] <- tissue_palette[palette_names]

ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = "anonymous"),
    tags$link(
      rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Fraunces:wght@600;700&family=IBM+Plex+Sans:wght@300;400;600&display=swap"
    )
  ),
  div(
    class = "hero",
    div(
      class = "hero-inner",
      div(
        class = "hero-text",
        h1(APP_TITLE),
        p("Explore tissue and time-specific proteomic shifts."),
        div(class = "version", paste0("App version ", APP_VERSION))
      )
    )
  ),
  div(
    class = "app-shell",
    sidebarLayout(
      sidebarPanel(
        class = "sidebar",
        div(
          class = "view-toggle-row",
          div(
            class = "view-toggle",
            radioButtons(
              "view_mode",
              NULL,
              choiceNames = list(tags$span("Pathway Viewer"), tags$span("Gene Viewer")),
              choiceValues = c("Pathway", "Gene"),
              selected = "Pathway",
              inline = TRUE
            )
          )
        ),
        conditionalPanel(
          "input.view_mode == 'Pathway'",
          div(
            class = "card",
            h3("Pathway tree"),
            p("Select a pathway to render the bubble plot."),
            div(
              class = "tree-wrapper",
              shinyTree("pathway_tree", search = FALSE, theme = "proton")
            )
          )
        ),
        conditionalPanel(
          "input.view_mode == 'Gene'",
          div(
            class = "card",
            h3("Gene list"),
            p("Enter one gene per line."),
            textAreaInput(
              "gene_input",
              NULL,
              value = "",
              placeholder = "Psma1\nPsma3\nPsmb4",
              rows = 6
            )
          )
        ),
        conditionalPanel(
          "input.view_mode == 'Gene'",
          div(
            class = "card",
            h3("Tissues"),
            uiOutput("gene_tissue_select")
          )
        ),
        conditionalPanel(
          "input.view_mode == 'Pathway'",
          div(
            class = "card",
            h3("Filters"),
            numericInput("padj_cutoff", "BH-adjusted P-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.005),
            numericInput("logfc_cutoff", "Absolute logFC cutoff:", value = 0, min = 0, step = 0.1)
          )
        )
      ),
      mainPanel(
        class = "main-panel",
        conditionalPanel(
          "input.view_mode == 'Pathway'",
          uiOutput("selected_pathway_name"),
          div(
            class = "plot-card",
            plotlyOutput("bubble_plot", height = "700px")
          )
        ),
        conditionalPanel(
          "input.view_mode == 'Gene'",
          uiOutput("gene_status"),
          uiOutput("gene_plot_container")
        )
      )
    ),
    div(class = "section-divider"),
    div(
      class = "footer-row",
      div(
        class = "card cite",
        h3("How to cite"),
        p(CITATION_TEXT),
        tags$a("View publication", href = PUBLICATION_URL, target = "_blank", rel = "noopener")
      ),
      div(
        class = "card data",
        h3("Data access"),
        p(DATA_ACCESS_PRIDE),
        tags$a("View PRIDE", href = PRIDE_URL, target = "_blank", rel = "noopener"),
        p(DATA_ACCESS_GITHUB),
        tags$a("View github", href = GITHUB_URL, target = "_blank", rel = "noopener")
      )
    )
  )
)

server <- function(input, output, session) {
  tree_data <- reactive({
    build_tree(pathway_hierarchy, pathway_list)
  })

  output$pathway_tree <- renderTree({
    tree_data()
  })

  output$gene_tissue_select <- renderUI({
    choice_names <- lapply(all_tissues, function(tissue) {
      tags$span(
        class = "tissue-pill",
        tissue,
        style = paste0("--tissue-color: ", tissue_colors[[tissue]], ";")
      )
    })

    div(
      class = "tissue-toggle",
      checkboxGroupInput(
        "gene_tissues",
        NULL,
        choiceNames = choice_names,
        choiceValues = all_tissues,
        selected = all_tissues
      )
    )
  })

  selected_pathway_id <- reactive({
    sel <- get_selected(input$pathway_tree, format = "classid")
    if (is.null(sel) || length(sel) == 0) return(NULL)
    attr(sel[[1]], "stid")
  })

  bubble_data <- reactive({
    pid <- selected_pathway_id()
    bd <- get_bubble_data(
      pid,
      protein_map,
      diff_expr,
      input$padj_cutoff,
      input$logfc_cutoff,
      all_tissues,
      all_times,
      all_regulations
    )

    max_count <- if (all(is.na(bd$Count))) 0 else max(bd$Count, na.rm = TRUE)
    list(bubble = bd, max_count = max_count)
  })

  output$selected_pathway_name <- renderUI({
    req(input$view_mode == "Pathway")
    pid <- selected_pathway_id()
    bd_all <- bubble_data()

    if (is.null(pid) || is.null(bd_all)) {
      div(h3("No pathway selected", class = "plot-title"))
    } else {
      pname <- pathway_list$PathwayName[pathway_list$PathwayID == pid][1]
      max_count <- bd_all$max_count

      div(
        h3(pname, class = "plot-title"),
        div(pid, class = "plot-subtitle"),
        div(paste0("Largest protein set: ", max_count), class = "plot-meta")
      )
    }
  })

  output$bubble_plot <- renderPlotly({
    req(input$view_mode == "Pathway")
    bd_all <- bubble_data()
    if (is.null(bd_all)) return(NULL)
    make_bubble_plot(bd_all$bubble)
  })

  gene_list <- reactive({
    parse_gene_input(input$gene_input)
  })

  gene_plot_data <- reactive({
    genes <- gene_list()
    tissues <- input$gene_tissues
    if (length(genes) == 0 || is.null(tissues) || length(tissues) == 0) return(NULL)

    plot_df <- diff_expr %>%
      dplyr::filter(PG.Genes %in% genes, Tissue %in% tissues)

    if (nrow(plot_df) == 0) return(NULL)
    plot_df
  })

  gene_plot_gene_count <- reactive({
    plot_df <- gene_plot_data()
    if (is.null(plot_df)) return(0)
    length(unique(plot_df$PG.Genes))
  })

  output$gene_status <- renderUI({
    req(input$view_mode == "Gene")
    genes <- gene_list()
    tissues <- input$gene_tissues

    if (length(genes) == 0) {
      return(div(h3("No genes selected", class = "plot-title")))
    }

    if (is.null(tissues) || length(tissues) == 0) {
      return(div(h3("No tissues selected", class = "plot-title")))
    }

    if (is.null(gene_plot_data())) {
      return(div(h3("No data for selected genes/tissues", class = "plot-title")))
    }

    NULL
  })

  output$gene_plot_container <- renderUI({
    req(input$view_mode == "Gene")
    gene_count <- gene_plot_gene_count()
    rows <- max(1, ceiling(gene_count / 2))
    row_height <- 320
    height <- rows * row_height

    div(
      class = "plot-card",
      plotlyOutput("gene_plot", height = paste0(height, "px"))
    )
  })

  output$gene_plot <- renderPlotly({
    req(input$view_mode == "Gene")
    plot_df <- gene_plot_data()
    if (is.null(plot_df)) return(NULL)
    make_gene_plot(plot_df, tissue_colors)
  })
}

REVIEW_MODE <- Sys.getenv("REVIEW_MODE", "0") == "1"

if (REVIEW_MODE) {
  credentials <- data.frame(
    user = Sys.getenv("REVIEW_USER", "username"),
    password = Sys.getenv("REVIEW_PASS", "password"),
    stringsAsFactors = FALSE
  )
  
  ui_secure <- secure_app(ui)
  
  server_secure <- function(input, output, session) {
    secure_server(check_credentials = check_credentials(credentials))
    server(input, output, session)
  }
  
  shinyApp(ui_secure, server_secure)
  
} else {
  shinyApp(ui, server)
}