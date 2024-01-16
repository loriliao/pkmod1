# Options for Spinner
require(shinycssloaders)
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)

# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  
  # Application title.
  headerPanel(h3("Regimen simulator"),windowTitle = "Regimen simulator"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view. The helpText function is also used to 
  # include clarifying text. Most notably, the inclusion of a 
  # submitButton defers the rendering of output until the user 
  # explicitly clicks the button (rather than doing it immediately
  # when inputs change). This is useful if the computations required
  # to render output are inordinately time-consuming.
  sidebarPanel(
    numericInput("Dose1", "Day 1, 1st Dose (mg):", 
                 300, min = 1, max = 500),
    numericInput("Dose2", "Day 1, 2nd Dose and Day 2 bid Dose (mg):", 
                 100, min = 1, max = 500),
    numericInput("Dose3", "Day 3-10 bid Dose (mg):", 
                 30, min = 1, max = 500),
    
    # fluidRow(
    #   column(9, sliderInput("ncyc1", "# Cycles of Maintenance:", 
    #                         min=10, max=50, value=14, step=1)),
    #   column(3, selectInput("reg1", "Regimen:", 
    #                         choices = c("QD", "BID","TID"),
    #                         selected="QD"))
    # ),
    
    h4("Plotting options:"),
    fluidRow(
      column(6, radioButtons("semilog", "Linear or semi-log:", 
                             choices=c("linear", "semi-log"),
                             selected = "linear")),
    ),
    actionButton(
      inputId = "submit_mod",
      label = "Run Model")
  ), 
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations. Note the use of the h4 function to provide
  # an additional header above each output section.
  mainPanel(
    h4("Median, 5th and 95th pctl concentration-time profile"),
    withSpinner(plotOutput("CpPlot"), type = 1),
    fluidRow(
      column(4, tableOutput("day3params")),
      column(4, tableOutput("ssparams"))
    )
  )
))

