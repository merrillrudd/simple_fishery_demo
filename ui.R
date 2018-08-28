### general population dynamics
### teaching tool
### Merrill Rudd
### August 2018


library(shiny)

shinyUI(fluidPage(
  
  #  Application title
  titlePanel("Selectivity, maturity, and lengths in the catch"),

             sidebarLayout(
               sidebarPanel(
                 numericInput("amax", "Age classes:", value=20, step=1), 
                 numericInput("Linf", "Linf:", value=60, step=1),
                 numericInput("k", "k:", value = 0.20, step=0.01),
                 numericInput("M", "Natural mortality:", value=0.25, step=0.01),
                 sliderInput("amat", "Age at 50% maturity:", value=2, min=1, max=25),
                 sliderInput("aselex", "Age at 50% selectivity:", value = 3, min=1, max=25),
                 sliderInput("dome", "Level of dome-shaped selectivity", value=0, min=0, max=3, step=1),
                 sliderInput("fish", "Fishing mortality rate", min=0, max=1, step=0.01, value=0),
                 numericInput("lsamp", "Sample size of length measurements:", value=1000)
                 ),
                mainPanel(
                  column(4, plotOutput("SelexMature")),
                  column(4, plotOutput("NumbersAtAge")),
                  column(4, plotOutput("plotLegend")),
                  column(6, plotOutput("VBGFplot")),
                  column(6, plotOutput("CatchAtLength"))
                )
               )
             
  
))
