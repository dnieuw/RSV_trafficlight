empty_reference_entries <- function() {
  data.table::data.table(
    reference = character(),
    study_title = character(),
    study_type = character(),
    population = character(),
    notes = character(),
    DOI = character(),
    drug = character(),
    subtype = character(),
    mutation_signature = character(),
    mutation_count = integer(),
    N = character(),
    percentage = character(),
    mutations = I(list())
  )
}

empty_reference_overview <- function() {
  data.table::data.table(
    drug = character(),
    subtype = character(),
    study_count = integer(),
    mutation_entry_count = integer(),
    unique_mutation_count = integer()
  )
}

empty_study_index <- function() {
  data.table::data.table(
    reference = character(),
    study_title = character(),
    drug = character(),
    subtype = character(),
    mutation_entry_count = integer(),
    unique_mutation_count = integer(),
    DOI = character()
  )
}

load_resistance_reference_index <- function(reference_path) {
  result <- list(
    loaded = FALSE,
    status = NULL,
    reference_path = reference_path,
    entries = empty_reference_entries(),
    overview = empty_reference_overview(),
    study_index = empty_study_index(),
    raw_data = NULL
  )

  raw_data <- jsonlite::fromJSON(reference_path, simplifyVector = FALSE)

  studies <- raw_data$studies
  metadata_fields <- c("reference", "study_title", "study_type", "population", "notes", "DOI", "study_setup")

  entry_rows <- list()
  index_rows <- list()
  entry_counter <- 1L

  for (study in studies) {
    drug_fields <- setdiff(names(study), metadata_fields) #Any other fields are drug entries (allows for adding additional drugs)

    for (drug_name in drug_fields) {
      drug_payload <- study[[drug_name]]
      subtype_names <- names(drug_payload)

      for (subtype_name in subtype_names) {
        subtype_entries <- drug_payload[[subtype_name]]
        unique_mutations <- unique(unlist(lapply(subtype_entries, function(entry) entry$mutations), use.names = FALSE))

        index_rows[[length(index_rows) + 1L]] <- data.table::data.table(
          reference = study$reference %||% NA_character_,
          study_title = study$study_title %||% NA_character_,
          drug = drug_name,
          subtype = subtype_name,
          mutation_entry_count = length(subtype_entries),
          unique_mutation_count = length(unique_mutations),
          DOI = study$DOI %||% NA_character_
        )

        for (entry in subtype_entries) {
          mutations <- unique(entry$mutations %||% character())
          entry_rows[[entry_counter]] <- data.table::data.table(
            reference = study$reference %||% NA_character_,
            study_title = study$study_title %||% NA_character_,
            study_type = study$study_type %||% NA_character_,
            population = study$population %||% NA_character_,
            notes = study$notes %||% NA_character_,
            DOI = study$DOI %||% NA_character_,
            drug = drug_name,
            subtype = subtype_name,
            mutation_signature = paste(mutations, collapse = " + "),
            mutation_count = length(mutations),
            N = as.character(entry$N %||% ""),
            percentage = as.character(entry$percentage %||% ""),
            mutations = list(mutations)
          )
          entry_counter <- entry_counter + 1L
        }
      }
    }
  }

  entries <- rbindlist(entry_rows, fill = TRUE)

  study_index <- unique(rbindlist(index_rows, fill = TRUE))

  overview <- entries[, .(
      study_count = data.table::uniqueN(reference),
      mutation_entry_count = .N,
      unique_mutation_count = data.table::uniqueN(unlist(mutations))
    ), by = .(drug, subtype)][order(drug, subtype)]

  result$loaded <- TRUE
  result$status <- paste(
    "Reference data loaded from",
    reference_path,
    sprintf("(%d study records, %d mutation entries).", length(studies), nrow(entries))
  )
  result$entries <- entries
  result$overview <- overview
  result$study_index <- study_index
  result$raw_data <- raw_data
  result
}
