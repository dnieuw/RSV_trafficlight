function(input, output, session) {
  shared_upload <- reactiveVal(NULL)
  shared_sequences <- reactiveVal(NULL)
  shared_traffic_evidence <- reactive(NULL)
  
  #### Render menu items conditional on outputs
  output$welcome_tab <- renderMenu({
    menuItem("Welcome", tabName = "welcome_tab", selected = T)
  })
  output$sequence_quality_tab <- renderMenu({
    req(isTruthy(shared_upload()))
    menuItem("Sequence Quality", tabName = "sequence_quality_tab")
  })
  output$traffic_light_tab <- renderMenu({
    req(isTruthy(shared_sequences()))
    menuItem("Traffic Light", tabName = "traffic_light_tab")
  })
  output$additional_info_tab <- renderMenu({
    req(isTruthy(shared_sequences()))
    menuItem("Additional Information", tabName = "additional_info_tab")
  })
  isolate({updateTabItems(session, "sidebar", "welcome_tab")}) #Somehow it doesn't automatically load the "selected" tab..
  
  #### Observers to switch tabs
  observe({
    req(isTruthy(shared_upload()))
    updateTabItems(session, "sidebar", "sequence_quality_tab")
  })
  
  observe({
    req(isTruthy(shared_sequences()))
    updateTabItems(session, "sidebar", "traffic_light_tab")
  })
  
  #### Module calls
  welcome_server("welcome", shared_upload = shared_upload)
  sequence_quality_server(
    "sequence_quality",
    shared_sequences = shared_sequences,
    shared_upload = shared_upload
  )
  shared_traffic_evidence <- traffic_light_server(
    "traffic_light",
    shared_sequences = shared_sequences
  )
  additional_info_server(
    "additional_info",
    shared_traffic_evidence = shared_traffic_evidence
  )
}
