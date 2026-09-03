# UI
sequence_quality_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        width = 12,
        title = "🔎 QC Output",
        status = "info",
        solidHeader = TRUE,
        shinycssloaders::withSpinner(
          DT::DTOutput(ns("qc_summary"))
        ),
        br(),
        uiOutput(ns("run_analysis_button")),
        br(),
        # Put download button and translate switch side by side:
        downloadButton(ns("download_F_button"), "Download F Sequences"), HTML("&nbsp;&nbsp;&nbsp;"),
        shinyWidgets::materialSwitch(ns("translate_switch"), "Translate", status = "primary", inline = T, right = T)
      )
    )
  )
}

# Helper functions
check_valid_stringset <- function(stringset) {
  allowed <- colnames(nucleotideSubstitutionMatrix())
  for (n in seq_along(stringset)) {
    present <- unique(str_split_1(as.character(stringset[n]),pattern = ''))
    if (all(present%in%allowed)){
      next
    }
    return(FALSE)
  }
  return(TRUE)
}

check_input <- function(fasta, ref_path, F_gene_path) {
    if (!check_valid_stringset(fasta)){
      showModal(modalDialog(
        title = "Invalid Sequence",
        paste0("The uploaded sequences contain invalid characters not in ", paste(colnames(nucleotideSubstitutionMatrix()), collapse = ", "), ". Please upload valid DNA sequences."),
        easyClose = TRUE,
        footer = NULL
      ))
      stop()
    }
    
    ref_seqs <- readDNAStringSet(ref_path)
    if (!check_valid_stringset(ref_seqs)){
      showModal(modalDialog(
        title = "Invalid Reference Sequence",
        paste0("The reference sequences contain invalid characters not in ", paste(colnames(nucleotideSubstitutionMatrix()), collapse = ", "), "."),
        easyClose = TRUE,
        footer = NULL
      ))
      stop()
    }
    
    F_gene_seqs <- readDNAStringSet(F_gene_path) # (maybe better to use a range in the original reference to extract the F region)
    if (!check_valid_stringset(F_gene_seqs)){
      showModal(modalDialog(
        title = "Invalid F Gene Sequence",
        paste0("The F gene sequences contain invalid characters not in ", paste(colnames(nucleotideSubstitutionMatrix()), collapse = ", "), "."),
        easyClose = TRUE,
        footer = NULL
      ))
      stop()
    }
    return(list(fasta = fasta, ref_seqs = ref_seqs, F_gene_seqs = F_gene_seqs))
}

check_reference_alignments <- function(fasta, ref_seqs) {
  index <- IndexSeqs(ref_seqs, K = 6L)
  seq_names <- names(fasta)
  total_sequences <- length(fasta)

  # Assign default names if they are missing, missing names should ideally not happen and be checked in the upload module. (can only happen if someone uploads a fasta without names, I think)
  # if (is.null(seq_names) || any(seq_names == "")) {
  #   seq_names <- ifelse(
  #     is.null(seq_names) | seq_names == "",
  #     paste0("seq_", seq_along(fasta)),
  #     seq_names
  #   )
  # }

  ref_results <- withProgress(
    message = sprintf("Processing sequences (0/%d)", total_sequences),
    detail = "Preparing reference alignments",
    value = 0,
    {
      lapply(seq_along(fasta), function(i) {
        subject_name <- seq_names[i]
        sequence <- fasta[i]

        setProgress(
          value = (i - 1) / total_sequences,
          message = sprintf("Processing sequence %d/%d", i, total_sequences),
          detail = sprintf("Full reference alignment: %s", subject_name)
        )

        index_search <- SearchIndex(sequence, index, scoreOnly = FALSE)
        aln <- AlignPairs(sequence, ref_seqs, index_search)
        
        # # Maybe useful, to force DECIPHER to perform a global alignment
        # index_search$Position <- mapply(function(x, y, z)
        #   cbind(matrix(0L, 4), x, matrix(c(y, y, z, z), 4)),
        #   index_search$Position,
        #   width(sequence)[index_search$Pattern] + 1L,
        #   width(ref_seqs)[index_search$Subject] + 1L,
        #   SIMPLIFY=FALSE)
        # 
        # global_aln <- AlignPairs(sequence, ref_seqs, index_search, type = "both")
        
        aln_summary <- as.data.table(aln)

        #Rename columns to less confusing names
        setnames(aln_summary, old = c("Pattern", "PatternStart","PatternEnd",
                                      "PatternGapPosition","PatternGapLength",
                                      "Subject", "SubjectStart", "SubjectEnd",
                                      "SubjectGapPosition", "SubjectGapLength"),
                              new = c("SubjectID", "SubjectStart","SubjectEnd",
                                      "SubjectGapPosition","SubjectGapLength",
                                      "ReferenceID", "ReferenceStart", "ReferenceEnd",
                                      "ReferenceGapPosition","ReferenceGapLength"))

        aln_summary[, SubjectID := rep(subject_name, 2)]
        aln_summary[, ReferenceID := names(ref_seqs)]
        aln_summary[, Identity := Matches / width(sequence)]
        aln_summary[, best_reference := ReferenceID[which.max(Identity)]]
        aln_summary[, ReferenceSize := width(ref_seqs[best_reference])]
        aln_summary[, best_identity := max(Identity)]
        aln_summary <- cbind(aln_summary, letterFrequency(sequence, DNA_ALPHABET))

        list(
          aln_summary = aln_summary,
          best_reference = aln_summary[, best_reference][1]
        )
      })
    }
  )
  return(ref_results)
}

check_F_gene_alignments <- function(fasta, F_gene_seqs, ref_results) {
  index_F <- IndexSeqs(F_gene_seqs, K = 6L)
  seq_names <- names(fasta)
  total_sequences <- length(fasta)
  
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3)

  f_gene_results <- withProgress(
    message = sprintf("Processing sequences (0/%d)", total_sequences),
    detail = "Preparing F gene alignments",
    value = 0,
    {
      lapply(seq_along(fasta), function(i) {
        subject_name <- seq_names[i]
        sequence <- fasta[i]

        setProgress(
          value = (i - 1) / total_sequences,
          message = sprintf("Processing sequence %d/%d", i, total_sequences),
          detail = sprintf("F gene alignment: %s", subject_name)
        )
        
        best_hit <- ref_results[[i]]$best_reference
        best_F <- F_gene_seqs[best_hit]
        aln <- pairwiseAlignment(sequence, best_F, type = "overlap",
                                  substitutionMatrix = mat, gapOpening = 16, gapExtension = 1.2)
        aln_trans <- AAStringSet(list(translate_safe(alignedPattern(aln)),translate_safe(alignedSubject(aln))))
        
        index <- IndexSeqs(best_F, K = 6L)
        index_search <- SearchIndex(sequence, index, scoreOnly = FALSE)

        index_search$Position <- mapply(function(x, y, z)
          cbind(matrix(0L, 4), x, matrix(c(y, y, z, z), 4)),
          index_search$Position,
          width(sequence)[index_search$Pattern] + 1L,
          width(best_F)[index_search$Subject] + 1L,
          SIMPLIFY=FALSE)

        global_aln <- AlignPairs(sequence, best_F, index_search, type = "both")
        aln_trans <- AAStringSet(list(translate_safe(global_aln$PatternAligned),translate_safe(global_aln$SubjectAligned)))
        
        # if (names(sequence) == "HRSV-A11 F unstable, last 25 aa removed  (1653 bp)"){
        #   browser()
        # }
        
        freq_table <- letterFrequency(alignedSubject(aln), DNA_ALPHABET)
        identity <- (nmatch(aln) / nchar(aln))
        hit_info_F_row <- data.table(freq_table)
        hit_info_F_row[, SubjectID := subject_name]
        hit_info_F_row[, identity := identity]

        setProgress(
          value = i / total_sequences,
          message = sprintf("Processing sequence %d/%d", i, total_sequences),
          detail = sprintf("Completed sequence: %s", subject_name)
        )

        list(
          hit_info_F_row = hit_info_F_row,
          aln_F = aln,
          aln_F_trans = aln_trans
        )
      })
    }
  )
  return(f_gene_results)
}

# Server
sequence_quality_server <- function(id, shared_sequences = NULL, shared_upload = NULL) {
  moduleServer(id, function(input, output, session) {

    sequence_qc <- reactive({
      req(shared_upload())
      
      # Load and check that the input sequences are valid DNA sequences
      input_data <- check_input(shared_upload(), "RSV_RefSeq.fasta", "RSV_RefSeq_F.fasta")
      fasta <- input_data$fasta
      ref_seqs <- input_data$ref_seqs
      F_gene_seqs <- input_data$F_gene_seqs

      # Perform pairwise alignment of each uploaded sequence against the reference sequences
      ref_results <- check_reference_alignments(fasta, ref_seqs)
      
      # Do the same for the F gene region, using the best reference sequence for each uploaded sequence
      F_gene_results <- check_F_gene_alignments(fasta, F_gene_seqs, ref_results)

      hit_info <- rbindlist(lapply(ref_results, `[[`, "aln_summary"))
      hits <- setNames(lapply(ref_results, `[[`, "aln_summary"), names(fasta))

      hit_info_F <- rbindlist(lapply(F_gene_results, `[[`, "hit_info_F_row"))
      hits_F <- setNames(lapply(F_gene_results, `[[`, "aln_F"), names(fasta))
      hits_F_trans <- setNames(lapply(F_gene_results, `[[`, "aln_F_trans"), names(fasta))

      return(list(hit_info, hits, hit_info_F, hits_F, hits_F_trans))
    })

    # Update the shared_sequences reactive value when the analysis button is pressed, this will also switch the user to the next tab in the UI
    observeEvent(input$run_analysis_button, {
      # Get QC check results
      qc_results <- req(sequence_qc())
      
      AAStringSet(sapply(sequence_qc()[[5]],function(x){x[[2]]}))
      browser()
      
      ## Tell the user if a sequence was ignored because of the quality
      # showNotification("Analysis completed (scaffold mode). Traffic Light and Additional Information are now available.", type = "message")
      
      shared_sequences(sequence_qc())
    })

    ## Not run: ------------------------------------
    # # In server.R:
    # output$downloadData <- downloadHandler(
    #   filename = function() {
    #     paste('data-', Sys.Date(), '.csv', sep='')
    #   },
    #   content = function(con) {
    #     write.csv(data, con)
    #   }
    # )
    # 
    # # In ui.R:
    # downloadLink('downloadData', 'Download')
    ## ---------------------------------------------

    # Add button press event to generate and download F protein multiple alignment
    output$download_F_button <- downloadHandler(
      filename = function() {
        paste0("F_protein_alignment_", Sys.Date(), ".fasta")
      },
      content = function(file) {
        hits_F <- sequence_qc()[[4]]
        fasta_F <- DNAStringSet(sapply(hits_F, function(x) unlist(alignedPattern(x))))
        if (input$translate_switch) {
          fasta_F <- AAStringSet(sapply(fasta_F, translate_safe))
        }
        writeXStringSet(fasta_F, file)
      }
    )

    output$qc_summary <- DT::renderDT({
      req(sequence_qc())
      table_data <- sequence_qc()[[1]]
      #Format percentages
      table_data[, best_identity := sprintf("%.2f%%", best_identity * 100)]
      table_data[, Identity := sprintf("%.2f%%", Identity * 100)]
      
      #Add genome coverage percentage column based on alignment length and number of N positions:
      table_data[, genome_coverage := (AlignmentLength - N)/max(ReferenceSize,AlignmentLength),]
      table_data[, genome_coverage := sprintf("%.2f%%", genome_coverage * 100)]
      
      #Create typing column based on best reference
      table_data[, subtype := subtype_mapping[best_reference]]

      #Reorder some columns
      setcolorder(table_data, c("SubjectID", "subtype", "best_reference", "best_identity", "ReferenceID", "Identity", "genome_coverage"))
      DT::datatable(
        table_data,
        rownames = FALSE,
        extensions = 'Scroller',
        options = list(
          dom = 't',
          scrollY = 500,
          scrollX = TRUE,
          scroller = TRUE,
          columnDefs = list(list(
            targets = c(2,4),
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 11 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 11) + '...</span>' : data;",
              "}")
          ))
        )
      )
    })

    output$run_analysis_button <- renderUI({
      req(sequence_qc())
      actionButton(session$ns("run_analysis_button"), "Continue", class = "btn-success")
    })
    
  })
}
