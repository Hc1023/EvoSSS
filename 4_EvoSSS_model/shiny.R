library(shiny)
library(ggplot2)
ui <- fluidPage(
  sliderInput("r1", "Select r1 value", min = 0.1, max = 0.5, value = 0.1, step = 0.1),
  sliderInput("r2", "Select r2 value", min = 0.1, max = 0.5, value = 0.1, step = 0.1),
  plotOutput("plot")
)

server <- function(input, output) {
  output$plot <- renderPlot({
    data <- expand.grid(r1 = input$r1, r2 = input$r2, strain = 1:2, x = seq(0, 10, length.out = 100))
    data$y <- with(data, r1*x + r2*x^2)
    ggplot(data, aes(x, y)) + geom_line() + facet_wrap(~ strain)
  })
}

shinyApp(ui, server)
