# UI
additional_info_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        width = 12,
        title = "Additional Information",
        status = "info",
        solidHeader = TRUE,
        textOutput(ns("reference_status"))
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "Reference Data Overview",
        status = "info",
        solidHeader = TRUE,
        DT::DTOutput(ns("reference_overview"))
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "Study Index",
        status = "info",
        solidHeader = TRUE,
        DT::DTOutput(ns("study_index"))
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "Per-Sequence Additional Information",
        status = "info",
        solidHeader = TRUE,
        DT::DTOutput(ns("sequence_info_table"))
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "Selected Sequence Additional Information",
        status = "info",
        solidHeader = TRUE,
        uiOutput(ns("selected_sequence_details"))
      )
    )
  )
}

# Server
additional_info_server <- function(id, shared_traffic_evidence = NULL) {
  moduleServer(id, function(input, output, session) {
    traffic_evidence <- reactive({
      req(shared_traffic_evidence)
      shared_traffic_evidence()
    })

    output$reference_status <- renderText({
      traffic_evidence()$reference$status
    })

    output$reference_overview <- DT::renderDT({
      DT::datatable(
        traffic_evidence()$reference$overview,
        rownames = FALSE,
        options = list(dom = "t", paging = FALSE, searching = FALSE, info = FALSE, scrollX = TRUE)
      )
    })

    output$study_index <- DT::renderDT({
      DT::datatable(
        traffic_evidence()$reference$study_index,
        rownames = FALSE,
        options = list(dom = "t", pageLength = 8, scrollX = TRUE, autoWidth = TRUE)
      )
    })

    per_sequence_info <- reactive({
      traffic_evidence()$additional_info
    })

    output$sequence_info_table <- DT::renderDT({
      DT::datatable(
        per_sequence_info(),
        rownames = FALSE,
        selection = "single",
        options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE)
      )
    })

    selected_sequence <- reactive({
      selected_row <- input$sequence_info_table_rows_selected
      if (is.null(selected_row) || length(selected_row) == 0) {
        return(NULL)
      }
      selected_id <- per_sequence_info()[selected_row, sequence_id]
      traffic_evidence()$sequence_results[[selected_id]]
    })

    output$selected_sequence_details <- renderUI({
      selected <- selected_sequence()
      if (is.null(selected)) {
        return(tags$p("Click a row in the table above to see additional information for that sequence."))
      }

      known_hit_ui <- if (length(selected$known_hits) == 0) {
        tags$p(tags$b("Matched Study Evidence: "), "None")
      } else {
        tagList(
          tags$p(tags$b("Matched Study Evidence:")),
          lapply(selected$known_hits, function(hit) {
            tags$div(
              style = "margin-bottom: 12px;",
              tags$p(tags$b("Mutation Set: "), hit$mutation_signature),
              tags$p(tags$b("Match Type: "), hit$match_type),
              tags$p(tags$b("Matched Mutations: "), collapse_values(hit$matched_mutations)),
              if (length(hit$unmatched_mutations) > 0) {
                tags$p(tags$b("Unmatched Combination Members: "), collapse_values(hit$unmatched_mutations))
              },
              tags$p(tags$b("Study: "), hit$reference),
              tags$p(tags$b("Study Title: "), hit$study_title),
              tags$p(tags$b("DOI: "), hit$DOI %||% "Unavailable"),
              tags$p(tags$b("Study Type: "), hit$study_type %||% "Unavailable"),
              tags$p(tags$b("Population: "), hit$population %||% "Unavailable"),
              if (!is.null(hit$N) && nzchar(hit$N)) {
                tags$p(tags$b("Observed Count (N): "), hit$N)
              },
              if (!is.null(hit$percentage) && nzchar(hit$percentage)) {
                tags$p(tags$b("Reported Percentage: "), hit$percentage)
              },
              if (!is.null(hit$notes) && nzchar(hit$notes)) {
                tags$p(tags$b("Notes: "), hit$notes)
              }
            )
          })
        )
      }

      tagList(
        tags$p(tags$b("Sequence ID: "), selected$sequence_id),
        tags$p(tags$b("Subtype: "), selected$subtype %||% "Unknown"),
        tags$p(tags$b("Traffic Label: "), selected$traffic_label),
        tags$p(tags$b("Matched Studies: "), collapse_values(selected$matched_studies)),
        tags$p(tags$b("Detected Amino-Acid Mutations: "), collapse_values(selected$detected_aa_mutations)),
        tags$p(tags$b("Novel Amino-Acid Mutations: "), collapse_values(selected$novel_aa_mutations)),
        tags$p(tags$b("Interpretation Note: "), selected$interpretation_note),
        known_hit_ui
      )
    })
  })
}
