# UI
welcome_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        width = 6,
        title = "🧬 Overview",
        status = "info",
        solidHeader = TRUE,
        p("This application helps researchers judge RSV mutation impact by:"),
        tags$ul(
          tags$li("Uploading new RSV sequences"),
          tags$li("Comparing to RSV reference strains"),
          tags$li("Extracting F region sequences"),
          tags$li("Locating F protein amino acid mutations"),
          tags$li("Interpreting and highlighting amino acid mutations based on available knowledge")
        )
      ),
      box(
        width = 6,
        title = "📋 Workflow Guide",
        status = "info",
        solidHeader = TRUE,
        tags$ol(
          tags$li(tags$b("Upload Sequence Data:"), " Use the text/file upload field to upload your RSV nucleotide sequences"),
          tags$li(tags$b("Verify Sequence Integrity:"), " In the 'Sequence Quality' tab view the quality control of your sequences"),
          tags$li(tags$b("View Mutation Interpretation:"), " In the 'Traffic Light' tab view the interpretation of mutations in your uploaded data"),
          tags$li(tags$b("(Optional) Review Additional Information:"), " Navigate to 'Additional Information' tab to get more informations about specific mutations, links to literature, and more details about the analysis")
        )
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "📤 Upload Sequence Data",
        status = "info",
        solidHeader = TRUE,
        textAreaInput(
          ns("sequence_text"),
          "Paste FASTA sequences:",
          value = "",
          placeholder = ">sequence_1\nACGT...",
          rows = 8
        ),
        uiOutput(ns("sequence_text_validation")),
        fileInput(ns("sequence_file"), "or Choose (or drag) FASTA file", accept = c(".fa", ".fasta", ".fna", ".txt")),
        actionButton(ns("start_button"), "Start", class = "btn-success"),
        actionButton(ns("run_example_button"), "Run example", class = "btn-default")
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "⚙️ Technical Details & References",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        tags$ul(
          tags$li(tags$b("R libraries used:"), "BioStrings"),
          tags$li(tags$b("Reference sequence:"), " RSV A (",tags$a("KY883567.1", href = "https://www.ncbi.nlm.nih.gov/nuccore/KY883567.1"),") (", tags$a("https://doi.org/10.1093/ve/vead006", href = "https://doi.org/10.1093/ve/vead006"),")"),
          tags$li(tags$b("Reference sequence:"), " RSV B (",tags$a("KY883569.1", href = "https://www.ncbi.nlm.nih.gov/nuccore/KY883569.1"),") (", tags$a("https://doi.org/10.1093/ve/vead006", href = "https://doi.org/10.1093/ve/vead006"),")")
        )
      )
    )
  )
}

# Server
welcome_server <- function(id, shared_upload = NULL) {
  moduleServer(id, function(input, output, session) {

    read_fasta_with_handling <- function(read_fun, input_source) {
      warning_message <- NULL
      sequences <- tryCatch(
        withCallingHandlers(
          read_fun(),
          warning = function(w) {
            warning_message <<- conditionMessage(w)
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) e
      )

      if (inherits(sequences, "error")) {
        showModal(modalDialog(
          title = "Error",
          paste0("Failed to read the ", input_source, ". Please check the FASTA format and make sure it is not empty."),
          tags$br(),
          HTML("<b>This error was given:</b>"),
          tags$br(),
          sequences[[1]],
          easyClose = TRUE
        ))
        return(NULL)
      }

      if (!is.null(warning_message)) {
        showModal(modalDialog(
          title = "Warning",          
          paste0("Succeeded to read the ", input_source, " with warnings. Please be careful and check the file format and sequences."),
          tags$br(),
          HTML("<b>This warning was given:</b>"),
          tags$br(),
          warning_message,
          easyClose = TRUE
        ))
      }

      if (is.null(sequences) || length(sequences) == 0) {
        showModal(modalDialog(
          title = "Error",
          paste0("Failed to read the ", input_source, ". Please check the FASTA format and make sure it is not empty."),
          easyClose = TRUE,
          footer = NULL
        ))
        return(NULL)
      } else {
        return(sequences)
      }
    }

    load_example_sequences <- function() {
      read_fasta_with_handling(
        function() Biostrings::readDNAStringSet("example_fasta/test_sequences.fasta"),
        "example FASTA file"
      )
    }
    
    parse_fasta_file <- function(file) {
      read_fasta_with_handling(
        function() Biostrings::readDNAStringSet(file$datapath),
        "uploaded FASTA file"
      )
    }

    parse_fasta_text <- function(text_value) {
      read_fasta_with_handling(
        function() {
          fasta_lines <- unlist(strsplit(text_value, "\n", fixed = TRUE), use.names = FALSE)
          fasta_lines <- trimws(fasta_lines)
          header_lines <-  fasta_lines[seq(1, length(fasta_lines), by = 2)]
          sequence_lines <- toupper(fasta_lines[seq(2, length(fasta_lines), by = 2)])
          stringset <- Biostrings::DNAStringSet(sequence_lines)
          names(stringset) <- header_lines
          return(stringset)
        },
        "pasted FASTA text"
      )
    }

    debounced_sequence_text <- debounce(reactive(input$sequence_text), millis = 10)

    output$sequence_text_validation <- renderUI({
      text_value <- debounced_sequence_text()
      if (!nzchar(trimws(text_value))) {
        return(NULL)
      }

      if (is_proper_fasta(text_value)) {
        tags$small(
          style = "display: block; margin-top: 6px; font-size: 0.85em; color: #2e7d32;",
          "Valid FASTA sequence"
        )
      } else {
        tags$small(
          style = "display: block; margin-top: 6px; font-size: 0.85em; color: #c62828;",
          "Invalid FASTA sequence"
        )
      }
    })

    is_proper_fasta <- function(text_value) {

      fasta_lines <- unlist(strsplit(text_value, "\n", fixed = TRUE), use.names = FALSE)
      fasta_lines <- trimws(fasta_lines)

      # Check that number of lines is a multiple of 2 (header + sequence)
      if (length(fasta_lines) %% 2 != 0) {
        return(FALSE)
      }
      
      # Check that headers start with '>' and does not end with '>'
      header_lines <-  fasta_lines[seq(1, length(fasta_lines), by = 2)]
      if (!all(startsWith(header_lines, ">")) || any(endsWith(header_lines, ">"))) {
        return(FALSE)
      }

      # Check that sequence lines are non-empty
      sequence_lines <- toupper(fasta_lines[seq(2, length(fasta_lines), by = 2)])
      if (!all(nzchar(sequence_lines))) {
        return(FALSE)
      }

      # Check that sequence lines contain only valid nucleotide characters (A, C, G, T, N, and IUPAC codes)
      valid_sequence_pattern <- "^[ACGTMRWSYKVHDBN]+$"
      if (!all(grepl(valid_sequence_pattern, sequence_lines))) {
        return(FALSE)
      }

      return(TRUE)
    }

    handle_start_button <- function() {
      text_value <- input$sequence_text

      if (!nzchar(trimws(text_value))) {
        showModal(modalDialog(
          title = "Missing FASTA text",
          "Please paste a valid FASTA sequence in the text field, or use the Run example button.",
          easyClose = TRUE,
          footer = NULL
        ))
        return(invisible(NULL))
      }

      if (!is_proper_fasta(text_value)) {
        showModal(modalDialog(
          title = "Invalid FASTA text",
          "The pasted text does not look like FASTA. Make sure it starts with a header line beginning with '>' and contains sequence lines with valid nucleotide characters.",
          easyClose = TRUE,
          footer = NULL
        ))
        return(invisible(NULL))
      }

      shared_upload(parse_fasta_text(text_value))
      invisible(NULL)
    }
    
    observeEvent(input$start_button, {
      handle_start_button()
    }, ignoreInit = TRUE)

    observeEvent(input$sequence_file, {
      req(input$sequence_file)
      shared_upload(parse_fasta_file(input$sequence_file))
    }, ignoreInit = TRUE)

    observeEvent(input$run_example_button, {
      shared_upload(load_example_sequences())
    }, ignoreInit = TRUE)
  })
}
