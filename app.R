# Load packages ----
library(shiny)

# Source helper functions -----
source("SILR.R")

# User interface ----
ui <- fluidPage(
  titlePanel("SILR"),
  
  sidebarLayout(
    sidebarPanel(width = 7, 
      helpText("SILR is abbreviated for single-individual likelihood ratios, which is a new exact test for demonstrating that an effect exists in binary trials."),
      
      helpText("References: "),
      
      helpText("Each participant in your study has undergone t binary trials, succeeding or failing on each trial. List t+1 non-negative integers showing the number of participants who achieved each number of successes from 0 to t. List them in that order, separated by commas."),
      
      textInput("hitfreq", label = "", value = "6, 10, 8, 5, 4, 2"), 
      
      helpText(HTML("Under your null hypothesis, what is the probability of success (a hit) on each trial? This one value (call it <i>pnull</i>), must apply to all trials for all participants.")),
      
      numericInput("pnull", label = "", value = 0.2),
      
      helpText(HTML("If you reject the null hypothesis that <i>pnull</i> is the true success rate for all participants, which of the following conclusions do you want to reach in this test? (You may run this program more than once in order to test for several different conclusions).")),
      
      helpText(HTML("upper-tail: The true long-term hit rate exceeds <i>pnull</i> for at least some participants.")),
      
      helpText(HTML("lower-tail: The true long-term hit rate falls below <i>pnull</i> for at least some participants.")),
      
      helpText(HTML("some-tail: The true hit rate deviates from <i>pnull</i> in some direction for at least some participants.")),
      
      selectInput("type", 
                  label= "",
                  choices = c("upper-tail", "lower-tail", "some-tail"),
                  selected = "upper-tail"),
      
      helpText(HTML("In normal use, this program finds a lower confidence limit on the number of participants in your study whose hit rate appears to differ from")),
      
      numericInput("alpha", label = "", value = 0.05),
      
      helpText(HTML("The program uses a 'rounding interval' <i>ri</i>. By default, <i>ri</i> is set to 0.01. Using a smaller value of <i>ri</i> might slightly increase the method's power, but using a larger value will speed the program up. Enter the value of <i>ri</i> you want to use.")),
      
      numericInput("ri", label = "", value = 0.01),
      
      actionButton("go", "Calculate")
    ),
    
    mainPanel(textOutput("selected_type"),
              textOutput("entered_vec"),
              htmlOutput("result"))
  )
)

# Server logic ----
server <- function(input, output) {
  # output$selected_type <- renderText({ 
  #   paste("You have selected", input$type)
  # })
  
  Vals <- eventReactive(input$go, {
    as.numeric(unlist(strsplit(input$hitfreq,",")))
  })
  
  # output$entered_vec <- renderText({
  #   paste("You have entered", Vals())
  # })

  output$result <- renderUI({ 
    hnum <- switch(input$type, 
                  "upper-tail" = 1,
                  "lower-tail" = 2,
                  "some-tail" = 3)
    
    alpha <- input$alpha
    
    pnull <- input$pnull
    
    ri <- input$ri
    
    # hitfreq <- as.numeric(unlist(strsplit(input$hitfreq,",")))
    
    ans = SILR(hnum = hnum, alpha = alpha, pnull = pnull, ri = ri, hitfreq = Vals())
    
    str1 <- paste("maxnsn = ", ans$maxnsn, ", ", "pstoohigh = ",  ans$pstoohigh)
    str2 <- paste("ntoolarge = ", ans$ntoolarge, ", ", "pstoolow = ", ans$pstoolow)
    str3 <- paste("withtrait = ", ans$withtrait)
    
    HTML(paste(str1, str2, str3, sep = '<br/>'))
    
  })  
  
}

# Run app ----
shinyApp(ui, server)