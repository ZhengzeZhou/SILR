# Load packages ----
library(shiny)

# Source helper functions -----
source("SILR.R")

# User interface ----
ui <- fluidPage(
  titlePanel("SILR"),
  
  sidebarLayout(
    sidebarPanel(width = 7, 
      helpText(HTML("SILR stands for 'single-individual likelihood ratios'. The SILR test can be used when each of several individuals has undergone <i>t</i> binary trials, each yielding a 'hit' or a 'miss', and you want to test whether some of those individuals have long-term hit rates that differ from a specified value <i>pnull</i>.")),
      
      helpText("References: "),
    
      helpText(HTML("List <i>t</i>+1 non-negative integers, separated by commas, showing, in order, the number of individuals who got each possible number of hits from 0 to <i>t</i>.")),
      
      textInput("hitfreq", label = "", value = "6, 10, 8, 5, 4, 2"), 
      
      helpText(HTML("Under your null hypothesis, what is the probability of success (a hit) on each trial? This one value (call it <i>pnull</i>), must apply to all trials for all participants.")),
      
      numericInput("pnull", label = "", value = 0.2, min = 0, max = 1.0, step = 0.01),
      
      helpText(HTML("If you reject the null hypothesis that <i>pnull</i> is the true success rate for all participants, which of the following conclusions do you want to reach in this test? (You may run this program more than once in order to test for several different conclusions).")),
      
      helpText(HTML("upper-tail: The true long-term hit rate exceeds <i>pnull</i> for at least some participants.")),
      
      helpText(HTML("lower-tail: The true long-term hit rate falls below <i>pnull</i> for at least some participants.")),
      
      helpText(HTML("some-tail: The true hit rate deviates from <i>pnull</i> in some direction for at least some participants.")),
      
      selectInput("type", 
                  label= "",
                  choices = c("upper-tail", "lower-tail", "some-tail"),
                  selected = "upper-tail"),
      
      helpText(HTML("Enter the value of <i>alpha</i> you want to use in finding confidence limits.")),
      
      numericInput("alpha", label = "", value = 0.05),
      
      helpText(HTML("The program uses a 'rounding interval' <i>ri</i>. By default, <i>ri</i> is set to 0.01. Using a smaller value of <i>ri</i> might slightly increase the method's power, but using a larger value will speed the program up. Enter the value of <i>ri</i> you want to use.")),
      
      numericInput("ri", label = "", value = 0.01),
      
      helpText(HTML("Do you want to find a confidence limit on <i>pmax</i> or <i>pmin</i>?")),
      
      selectInput("pmax", 
                  label= "",
                  choices = c("yes", "no"),
                  selected = "yes"),
      
      actionButton("go", "Calculate")
    ),
    
    mainPanel(textOutput("selected_type"),
              textOutput("entered_vec"),
              h4("Input: "),
              htmlOutput("input"),
              h4("Results: "),
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
  
  output$input <- renderUI({ 
    hnum <- switch(input$type, 
                   "upper-tail" = 1,
                   "lower-tail" = 2,
                   "some-tail" = 3)
    
    str_N <- paste("number of participants <i>N</i> = ", sum(Vals()))
    str_type <- paste("type of test = ", input$type)
    str_t <- paste("number of binary trials <i>t</i> for each participant = ", length(Vals()) - 1)
    str_pnull <- paste("<i>pnull</i> = ", input$pnull)
    str_ri <- paste("rounding interval <i>ri</i> = ", input$ri)
    str_alpha <- paste("significance level <i>alpha</i> = ", input$alpha)
    str_hitfreq <- paste("number of participants with each possible number of hits = ", input$hitfreq)
    
    
    HTML(paste(str_type, str_N, str_t, str_pnull, str_ri, str_alpha, str_hitfreq, sep = '<br/>'))
    
  })  

  output$result <- renderUI({ 
    hnum <- switch(input$type, 
                  "upper-tail" = 1,
                  "lower-tail" = 2,
                  "some-tail" = 3)
    
    alpha <- input$alpha
    
    pnull <- input$pnull
    
    ri <- input$ri
    
    pmax_yn <- switch(input$pmax, 
                   "yes" = 1,
                   "no" = 0)
    
    # hitfreq <- as.numeric(unlist(strsplit(input$hitfreq,",")))
    
    
    full <- SILR(hnum = hnum, alpha = 0, pnull = pnull, ri = ri, hitfreq = Vals())
    
    if (hnum == 3 & pmax_yn == 1){
      if (full$pstoohigh > alpha) {
        str_full <- paste("full sample ps = ", signif(full$pstoohigh, digits = 6))
        str <- "Confidence limits on pmax and pmin are not meaningful for a some-tail test."
        
        HTML(paste(str_full, str, sep = '<br/>')) 
      } else {
        str_full <- paste("full sample ps = ", signif(full$pstoohigh, digits = 6))
        str <- "Confidence limits on pmax and pmin are not meaningful for a some-tail test."
        
        ans = SILR(hnum = hnum, alpha = alpha, pnull = pnull, ri = ri, hitfreq = Vals())
        str3 <- paste("Lower confidence limit on the number of experimental participants with the trait of interest = ", ans$withtrait)
        
        HTML(paste(str_full, str, str3, sep = '<br/>')) 
      }
    
    } else if (full$pstoohigh > alpha) {
      str_full <- paste("full sample ps = ", signif(full$pstoohigh, digits = 6))
      str <- "No confidence limits were printed because <i>ps</i> > <i>alpha</i>."
      
      HTML(paste(str_full, str, sep = '<br/>'))          
                  
    } else{
      
      if (pmax_yn == 0) {
        str_full <- paste("full sample ps = ", signif(full$pstoohigh, digits = 6))
        
        ans = SILR(hnum = hnum, alpha = alpha, pnull = pnull, ri = ri, hitfreq = Vals())
        
        # str1 <- paste("maxnsn = ", ans$maxnsn, ", ", "pstoohigh = ",  ans$pstoohigh)
        # str2 <- paste("ntoolarge = ", ans$ntoolarge, ", ", "pstoolow = ", ans$pstoolow)
        str3 <- paste("Lower confidence limit on the number of experimental participants with the trait of interest = ", ans$withtrait)
        
        # HTML(paste(str_full, str1, str2, str3, sep = '<br/>'))
        HTML(paste(str_full, str3, sep = '<br/>'))
      } else {
        str_full <- paste("full sample ps = ", signif(full$pstoohigh, digits = 6))
        
        ans = SILR(hnum = hnum, alpha = alpha, pnull = pnull, ri = ri, hitfreq = Vals())
        
        # str1 <- paste("maxnsn = ", ans$maxnsn, ", ", "pstoohigh = ",  ans$pstoohigh)
        # str2 <- paste("ntoolarge = ", ans$ntoolarge, ", ", "pstoolow = ", ans$pstoolow)
        str3 <- paste("lower confidence limit on the number of experimental participants with the trait of interest = ", ans$withtrait)
        
        ans = findPmaxlog(hnum, alpha, pnull = pnull, ri = ri, printallps = 0, hitfreq = Vals(), tol = 1e-2)
        
        if (hnum == 1) {
          pmax = ans$xv[length(ans$xv)]
          
          str4 <- paste("lower confidence limit on pmax = ", round(pmax, digits = 6))
        } else {
          pmin = ans$xv[length(ans$xv)]
          
          str4 <- paste("upper confidence limit on pmin = ", round(pmin, digits = 6))
        }
        
        # HTML(paste(str_full, str1, str2, str3, str4, sep = '<br/>'))
        HTML(paste(str_full, str3, str4, sep = '<br/>'))
      }

    }
    
  })  
  
}

# Run app ----
shinyApp(ui, server)