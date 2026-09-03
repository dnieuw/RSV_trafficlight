# UI
traffic_light_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        width = 12,
        title = "🚥 Per-Sequence Result",
        status = "info",
        solidHeader = TRUE,
        shinycssloaders::withSpinner(
          DT::DTOutput(ns("traffic_table"))
        )
      )
    ),
    fluidRow(
      box(
        width = 12,
        title = "ℹ️ Selected Sequence Details",
        status = "info",
        solidHeader = TRUE,
        shinycssloaders::withSpinner(
          uiOutput(ns("selected_sequence_info"))
        )
      )
    )
  )
}

# Helper functions
translate_safe <- function(dna) {
  normalized_nt <- toupper(gsub("[^ACGTMRWSYKVHDB]", "N", dna))
  AA_sequence <- Biostrings::translate(
    Biostrings::DNAString(normalized_nt),
    if.fuzzy.codon = "solve"
  )
  return(AA_sequence)
}

collapse_values <- function(values, empty_label = "None") {
  if (is.null(values) || length(values) == 0) {
    return(empty_label)
  }

  paste(unique(values), collapse = ", ")
}

collapse_mutation_labels <- function(mutations, empty_label = "None") {
  if (is.null(mutations) || nrow(mutations) == 0) {
    return(empty_label)
  }

  collapse_values(mutations$mutation, empty_label = empty_label)
}

build_sequence_metadata <- function(sequence_bundle) {
  hit_info <- sequence_bundle[[1]]

  metadata <- data.table::as.data.table(hit_info)[, .(
    best_reference = best_reference[1]
  ), by = .(sequence_id = SubjectID)]
 
  if (!all(metadata$best_reference %in% names(subtype_mapping))) {
    modalDialog(
      title = "Reference Picking Error",
      paste0("Some of the best matching reference sequences are not in the subtype mapping, 
      did you add a new reference sequence to the reference file? 
      If so, update traffic_light_module.R"),
      easyClose = TRUE,
      footer = NULL
    )
    stop("Reference mapping error: Some best_reference values are not in the subtype mapping.")
  }

  metadata[, subtype := subtype_mapping[best_reference]]
  metadata
}

match_known_hits <- function(detected_aa_mutations, subtype, reference_index, drug) {
  entries <- reference_index$entries
  selected_drug <- drug
  selected_subtype <- subtype

  if (length(detected_aa_mutations) == 0) {
    return(list(matches = list(), relevant_entries = entries[0]))
  }

  relevant_entries <- entries[drug == selected_drug & subtype == selected_subtype]

  matches <- lapply(seq_len(nrow(relevant_entries)), function(index) {
    entry <- relevant_entries[index]
    entry_mutations <- unique(unlist(entry$mutations, use.names = FALSE))
    matched_mutations <- intersect(entry_mutations, detected_aa_mutations)

    if (length(entry_mutations) == 1 && length(matched_mutations) == 1) {
      match_type <- "single"
    } else if (length(entry_mutations) > 1 && length(matched_mutations) == length(entry_mutations)) {
      match_type <- "combination_full"
    } else if (length(entry_mutations) > 1 && length(matched_mutations) > 0) {
      match_type <- "combination_multiple"
    } else {
      return(NULL)
    }

    list(
      reference = entry$reference[[1]],
      study_title = entry$study_title[[1]],
      study_type = entry$study_type[[1]],
      population = entry$population[[1]],
      notes = entry$notes[[1]],
      DOI = entry$DOI[[1]],
      drug = entry$drug[[1]],
      subtype = entry$subtype[[1]],
      mutation_signature = entry$mutation_signature[[1]],
      mutation_count = entry$mutation_count[[1]],
      matched_mutations = matched_mutations,
      unmatched_mutations = setdiff(entry_mutations, matched_mutations),
      match_type = match_type,
      N = entry$N[[1]],
      percentage = entry$percentage[[1]]
    )
  })

  matches <- Filter(Negate(is.null), matches)
  list(matches = matches, relevant_entries = relevant_entries)
}

build_drug_assessment <- function(detected_aa_mutations, nt_mutations, subtype, reference_index, drug) {
  known_hit_result <- match_known_hits(
    detected_aa_mutations = detected_aa_mutations,
    subtype = subtype,
    reference_index = reference_index,
    drug = drug
  )

  known_hits <- known_hit_result$matches
  relevant_entries <- known_hit_result$relevant_entries
  known_mutation_universe <- if (nrow(relevant_entries) == 0) {
    character()
  } else {
    unique(unlist(relevant_entries$mutations, use.names = FALSE))
  }
  novel_aa_mutations <- setdiff(detected_aa_mutations, known_mutation_universe)
  has_indel <- any(nt_mutations$mutation_type %in% c("insertion", "deletion"))

  traffic_label <- if (length(detected_aa_mutations) == 0) {
    "🟢"
  } else if (length(known_hits) > 0) {
    "🔴"
  } else {
    "🟡"
  }

  interpretation_note <- build_interpretation_note(
    detected_aa_mutations = detected_aa_mutations,
    known_hits = known_hits,
    novel_aa_mutations = novel_aa_mutations,
    has_indel = has_indel,
    subtype = subtype
  )

  matched_studies <- unique(vapply(known_hits, function(hit) hit$reference, character(1)))

  list(
    drug = drug,
    known_hits = known_hits,
    relevant_entries = relevant_entries,
    novel_aa_mutations = novel_aa_mutations,
    traffic_label = traffic_label,
    interpretation_note = interpretation_note,
    known_hit_summary = build_known_hit_summary(known_hits),
    novel_mutation_summary = collapse_values(novel_aa_mutations),
    known_hit_count = length(known_hits),
    matched_study_count = length(matched_studies),
    matched_studies = matched_studies,
    has_indel = has_indel
  )
}

build_known_hit_summary <- function(known_hits) {
  if (length(known_hits) == 0) {
    return("None")
  }
  
  all_matches <- unique(as.vector(sapply(known_hits, function(hit) {
    return(hit$matched_mutations)
  })))

  return(paste(all_matches, collapse = ","))
  
  # labels <- vapply(known_hits, function(hit) {
  #   if (identical(hit$match_type, "single")) {
  #     return(hit$mutation_signature)
  #   }
  # 
  #   prefix <- if (identical(hit$match_type, "combination_full")) {
  #     "Full"
  #   } else {
  #     paste("Partial", collapse_values(hit$matched_mutations))
  #   }
  # 
  #   paste(prefix, hit$mutation_signature)
  # }, character(1))
  # 
  # collapse_values(labels)
}

build_interpretation_note <- function(detected_aa_mutations, known_hits, novel_aa_mutations, has_indel, subtype) {
  note <- if (length(detected_aa_mutations) == 0) {
    paste("No amino-acid differences detected for", subtype %||% "the selected subtype", "against the nirsevimab evidence set.")
  } else if (length(known_hits) > 0) {
    paste(length(known_hits), "known study hit(s) detected for", subtype %||% "the selected subtype", ".")
  } else {
    paste(length(novel_aa_mutations), "novel amino-acid mutation(s) detected for", subtype %||% "the selected subtype", ".")
  }

  if (has_indel) {
    note <- paste(note, "Nucleotide insertion/deletion also detected; review alignment manually.")
  }

  note
}

build_sequence_evidence <- function(sequence_id, mutation_bundle, sequence_metadata, reference_index, drug_panel) {
  ###This is still a bit AI sloppy:
  selected_sequence_id <- sequence_id
  metadata_row <- sequence_metadata[sequence_id == selected_sequence_id]
  subtype <- if (nrow(metadata_row) == 0) NA_character_ else metadata_row$subtype[[1]]
  best_reference <- if (nrow(metadata_row) == 0) NA_character_ else metadata_row$best_reference[[1]]
  nt_mutations <- mutation_bundle$nt_mutations
  aa_mutations <- mutation_bundle$aa_mutations
  detected_aa_mutations <- unique(aa_mutations$mutation)
  ###

  drug_assessments <- lapply(drug_panel, function(drug_name) {
    build_drug_assessment(
      detected_aa_mutations = detected_aa_mutations,
      nt_mutations = nt_mutations,
      subtype = subtype,
      reference_index = reference_index,
      drug = drug_name
    )
  })
  names(drug_assessments) <- drug_panel
  primary_drug_assessment <- drug_assessments[["nirsevimab"]]
  combined_known_mutation_universe <- unique(unlist(lapply(drug_assessments, function(assessment) {
    if (nrow(assessment$relevant_entries) == 0) {
      return(character())
    }

    unique(unlist(assessment$relevant_entries$mutations, use.names = FALSE))
  }), use.names = FALSE))
  combined_novel_aa_mutations <- setdiff(detected_aa_mutations, combined_known_mutation_universe)
  combined_novel_mutation_summary <- collapse_values(combined_novel_aa_mutations)

  summary_row <- data.table::data.table(
    sequence_id = sequence_id,
    subtype = subtype,
    nirsevimab = drug_assessments[["nirsevimab"]]$traffic_label,
    clesrovimab = drug_assessments[["clesrovimab"]]$traffic_label,
    palivizumab = drug_assessments[["palivizumab"]]$traffic_label,
    nirsevimab_mutation = drug_assessments[["nirsevimab"]]$known_hit_summary,
    clesrovimab_mutation = drug_assessments[["clesrovimab"]]$known_hit_summary,
    palivizumab_mutation = drug_assessments[["palivizumab"]]$known_hit_summary,
    novel_mutations = combined_novel_mutation_summary
  )

  list(
    sequence_id = sequence_id,
    subtype = subtype,
    best_reference = best_reference,
    drug = "nirsevimab",
    drug_assessments = drug_assessments,
    nt_mutations = nt_mutations,
    aa_mutations = aa_mutations,
    all_mutations = mutation_bundle$all_mutations,
    detected_nt_mutations = unique(nt_mutations$mutation),
    detected_aa_mutations = detected_aa_mutations,
    known_hits = primary_drug_assessment$known_hits,
    novel_aa_mutations = combined_novel_aa_mutations,
    traffic_label = primary_drug_assessment$traffic_label,
    interpretation_note = primary_drug_assessment$interpretation_note,
    known_hit_summary = primary_drug_assessment$known_hit_summary,
    novel_mutation_summary = combined_novel_mutation_summary,
    known_hit_count = primary_drug_assessment$known_hit_count,
    matched_study_count = primary_drug_assessment$matched_study_count,
    matched_studies = primary_drug_assessment$matched_studies,
    has_indel = primary_drug_assessment$has_indel,
    summary_row = summary_row,
    additional_info_row = data.table::data.table(
      sequence_id = sequence_id,
      subtype = subtype,
      traffic_label = primary_drug_assessment$traffic_label,
      known_hit_count = primary_drug_assessment$known_hit_count,
      novel_mutation_count = length(combined_novel_aa_mutations),
      matched_study_count = primary_drug_assessment$matched_study_count,
      note = primary_drug_assessment$interpretation_note
    )
  )
}

extract_nt_mutations <- function(alignment, sequence_id) {
  subject_aligned <- alignedSubject(alignment)
  pattern_aligned <- alignedPattern(alignment)
  total_positions <- width(subject_aligned)
  
  check_position <- function(index) {
    ref_base <- as.character(subseq(subject_aligned, start = index, width = 1))
    query_base <- as.character(subseq(pattern_aligned, start = index, width = 1))
    
    if (ref_base != "-" && query_base != "-") {
      if (ref_base != query_base) {
        # Substitution
        return(data.table(
          sequence_id = sequence_id,
          level = "nt",
          mutation_type = "substitution",
          position = index,
          reference = ref_base,
          observed = query_base,
          mutation = paste0(ref_base, index, query_base)
        ))
      } else {
        return(NULL)
      }
    } else if (ref_base == "-") {
      # Insertion
      return(data.table(
        sequence_id = sequence_id,
        level = "nt",
        mutation_type = "insertion",
        position = index,
        reference = "-",
        observed = query_base,
        mutation = paste0(ref_base, index, query_base)
      ))
    } else if (query_base == "-") {
      # Deletion
      return(data.table(
        sequence_id = sequence_id,
        level = "nt",
        mutation_type = "deletion",
        position = index,
        reference = ref_base,
        observed = "-",
        mutation = paste0(ref_base, index, "-")
      ))
    } else {
      stop(paste0("Unexpected case encountered in mutation extraction. Reference base: ", ref_base, ", Query base: ", query_base))
    }
  }
  
  mutation_table <- rbindlist(lapply(seq(total_positions), check_position))
  
  return(mutation_table)
}

extract_aa_mutations <- function(alignment, nt_mutations, sequence_id) {
  #Return empty dt if there are no mutations
  if (nrow(nt_mutations) == 0){
    return(data.table())
  }
  
  subject_aligned <- alignedSubject(alignment)
  pattern_aligned <- alignedPattern(alignment)
  
  # Check which codons have nucleotide mutations
  mutated_codon_indices <- unique((nt_mutations$position - 1) %/% 3 + 1)
  mutated_codon_starts <- (mutated_codon_indices - 1) * 3 + 1
  
  subject_translated <- translate_safe(subject_aligned)
  pattern_translated <- translate_safe(pattern_aligned)
  
  #FIX: Return empty dt if there are extra stopcodons in the sequence, should just ignore a faulty F protein sequence
  subject_count_trimmed <- Biostrings::countPattern("*", Biostrings::trimLRPatterns("","*",subject_translated))
  pattern_count_trimmed <- Biostrings::countPattern("*", Biostrings::trimLRPatterns("","*",pattern_translated))
  if (subject_count_trimmed > 0 | pattern_count_trimmed > 0){
    return(data.table())
  }
  
  check_position <- function(index) {
    
    ref_codon <- subseq(subject_aligned, start = index, width = 3)
    query_codon <- subseq(pattern_aligned, start = index, width = 3)
    
    aa_position <- index %/% 3 + 1
    
    ref_AA <- as.character(subseq(subject_translated, aa_position, width = 1))
    query_AA <- as.character(subseq(pattern_translated, aa_position, width = 1))
    
    #Ignore synonymous or indel or stop codon mutations
    if (ref_AA == query_AA | ref_AA == "X" | query_AA == "X" | ref_AA == "*" | query_AA == "*") {
      return(NULL)
    }
    
    return(data.table(
      sequence_id = sequence_id,
      level = "aa",
      mutation_type = "substitution",
      position = aa_position,
      reference = ref_AA,
      observed = query_AA,
      mutation = paste0(ref_AA, aa_position, query_AA),
      reference_codon = as.character(ref_codon),
      observed_codon = as.character(query_codon)
    ))
  }

  mutation_table <- rbindlist(lapply(mutated_codon_starts, check_position))

  return(mutation_table)
}

compare_pairwise_alignment_mutations <- function(alignment, sequence_id) {
  nt_mutations <- extract_nt_mutations(
    alignment,
    sequence_id = sequence_id
  )

  aa_mutations <- extract_aa_mutations(
    alignment, nt_mutations,
    sequence_id = sequence_id
  )

  all_mutations <- data.table::rbindlist(
    list(nt_mutations, aa_mutations),
    fill = TRUE
  )

  list(
    nt_mutations = nt_mutations,
    aa_mutations = aa_mutations,
    all_mutations = all_mutations
  )
}

# Server
traffic_light_server <- function(id, shared_sequences = NULL) {
  moduleServer(id, function(input, output, session) {

    traffic_evidence <- reactive({
      sequence_bundle <- req(shared_sequences())
      alignments <- sequence_bundle[[4]]

      sequence_metadata <- build_sequence_metadata(sequence_bundle)
      # This fuction is loaded in global from the R/resistance_reference_helpers.R file
      reference_data <- load_resistance_reference_index(rsv_nirsevimab_resistance_table_path)

      sequence_ids <- names(alignments)
      total_sequences <- length(sequence_ids)

      per_sequence_results <- withProgress(
        message = sprintf("Analyzing mutations (0/%d)", total_sequences),
        detail = "Preparing sequence analysis",
        value = 0,
        {
          lapply(seq_along(sequence_ids), function(i) {
            sequence_id <- sequence_ids[[i]]
            sequence_progress_end <- i / max(total_sequences, 1)

            mutation_bundle <- compare_pairwise_alignment_mutations(
              alignments[[sequence_id]],
              sequence_id = sequence_id
            )

            sequence_evidence <- build_sequence_evidence(
              sequence_id = sequence_id,
              mutation_bundle = mutation_bundle,
              sequence_metadata = sequence_metadata,
              reference_index = reference_data,
              drug_panel = c("nirsevimab", "clesrovimab", "palivizumab")
            )

            setProgress(
              value = sequence_progress_end,
              message = sprintf("Analyzing mutations (%d/%d)", i, total_sequences),
              detail = sprintf("Completed sequence: %s", sequence_id)
            )

            sequence_evidence
          })
        }
      )

      names(per_sequence_results) <- names(alignments)

      summary_table <- data.table::rbindlist(
        lapply(per_sequence_results, `[[`, "summary_row"),
        fill = TRUE
      )
      additional_info_table <- data.table::rbindlist(
        lapply(per_sequence_results, `[[`, "additional_info_row"),
        fill = TRUE
      )
      mutations_table <- data.table::rbindlist(
        lapply(per_sequence_results, function(result) result$all_mutations),
        fill = TRUE
      )

      list(
        summary = summary_table,
        additional_info = additional_info_table,
        mutations = mutations_table,
        sequence_results = per_sequence_results,
        reference = reference_data
      )
    })

    selected_sequence_result <- reactive({
      selected_row <- input$traffic_table_rows_selected
      if (is.null(selected_row) || length(selected_row) == 0) {
        return(NULL)
      }

      selected_id <- traffic_evidence()$summary[selected_row, sequence_id]
      traffic_evidence()$sequence_results[[selected_id]]
    })

    output$traffic_table <- DT::renderDT({
      traffic_dt <- DT::datatable(
        traffic_evidence()$summary,
        rownames = FALSE,
        selection = "single",
        extensions = 'Scroller',
        options = list(
          dom = 't',
          scrollY = 500,
          scrollX = TRUE,
          scroller = TRUE
        )
      )
    })

    output$selected_sequence_info <- renderUI({
      selected <- selected_sequence_result()
      if (is.null(selected)) {
        return(tags$p("Click a row in the table above to see details for that sequence."))
      }

      known_hit_ui <- if (length(selected$known_hits) == 0) {
        tags$p(tags$b("Known Study Hits: "), "None")
      } else {
        tagList(
          tags$p(tags$b("Known Study Hits:")),
          tags$ul(lapply(selected$known_hits, function(hit) {
            label <- if (identical(hit$match_type, "combination_partial")) {
              paste0(
                hit$mutation_signature,
                " (partial: ",
                collapse_values(hit$matched_mutations),
                ")"
              )
            } else {
              hit$mutation_signature
            }

            tags$li(
              tags$b(label),
              tags$span(paste0(" | ", hit$reference)),
              if (!is.na(hit$DOI) && nzchar(hit$DOI)) {
                tags$span(paste0(" | ", hit$DOI))
              }
            )
          }))
        )
      }

      tagList(
        tags$p(tags$b("Sequence ID: "), selected$sequence_id),
        tags$p(tags$b("Subtype: "), selected$subtype %||% "Unknown"),
        tags$p(tags$b("Best Reference: "), selected$best_reference %||% "Unknown"),
        tags$p(tags$b("Traffic Label: "), selected$traffic_label),
        tags$p(tags$b("Detected Nucleotide Mutations: "), collapse_values(selected$detected_nt_mutations)),
        tags$p(tags$b("Detected Amino-Acid Mutations: "), collapse_values(selected$detected_aa_mutations)),
        tags$p(tags$b("Novel Amino-Acid Mutations: "), collapse_values(selected$novel_aa_mutations)),
        tags$p(tags$b("Interpretation Note: "), selected$interpretation_note),
        known_hit_ui
      )
    })

    traffic_evidence
  })
}
