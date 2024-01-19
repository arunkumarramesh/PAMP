library(ggplot2)
library(reshape2)
library(shiny)

cpm_log2 <- readRDS("cpm_log.RDS")

server <- function(input, output) {
  gene <- eventReactive(input$do, {
    unlist(strsplit(as.character(input$gene), ',', fixed=TRUE))
  }, ignoreNULL= T)
  waspexp <- reactive({
    cbind(cpm_log2["group"],cpm_log2[colnames(cpm_log2) %in% c(gene())])
  })
  waspexp2 <- reactive({
    melt(waspexp(),id.vars = "group")
  })
  p <- reactive({
    ggplot(data=waspexp2(), aes_string(x=waspexp2()$group, y=waspexp2()$value)) +
      geom_boxplot() +
      facet_wrap(~variable,ncol=3)+
      ylab("log(CPM)")+
      xlab("")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  })
  observeEvent(input$do, {
    output$plot1 <- renderPlot({
      print(p())
    })
  })
}

ui <- fluidPage(
  headerPanel('Wasp-induced gene expression browser'),
  sidebarPanel(
    textAreaInput('gene', 'Input gene names separated by comma:', value = "", width = NULL, placeholder = 'e.g. FBgn0016075,FBgn0000299'),
    actionButton("do", "Evaluate!")
  ),
  mainPanel(
    p("Gene expression in two immune tissues in Drosophila melanogaster larvae: Fat body and Hemocytes (blood cells)"),
    p("Control: Untreated"),
    p("Oil: Immune induction through injection of clear oil droplet into Drosophila"),
    p("Wasp extract: Immune induction through injection of oil droplet containing parasitic wasp (Leptopilina boulardi) tissue into Drosophila"),
    p("CPM: Counts per million, sequencing depth-normalized read counts"),
    p("Full data available at https://catalogue.ceh.ac.uk/documents/06ea87f3-476d-40fd-acce-e6923e786d48"),
    p("Please post questions at https://github.com/arunkumarramesh/PAMP/"),
    plotOutput("plot1")
  )
)

shinyApp(ui = ui, server = server)
