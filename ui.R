library(shiny)

shinyUI(
  pageWithSidebar(
    headerPanel(""),
    sidebarPanel(
      uiOutput("cut_off"),
      uiOutput("select_cut")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "Semantic Plot", 
          plotOutput("plot")
        ),
        tabPanel(
          title = "Table",
          h2('FPKM'),
          tableOutput('fpkm_table'),
          h2('Probes'),
          tableOutput('probe_table')
        ),
        tabPanel(
          title = "Odds Ratio",
          plotOutput("or_plot")
        )
        
      )
    )
  )
)