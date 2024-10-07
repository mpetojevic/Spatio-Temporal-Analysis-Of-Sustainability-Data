library(dplyr)
library(tidyr)
library(stringr)
library(rgeoboundaries)
library(class)

getwd()  # Check the current working directory
#setwd("~/Uni/2024S/Bachelor/Spatio-Temporal-Analysis-Of-Sustainability-Data/Cluster-Analysis") 
getwd()

# load data
data <- read.csv("sustainability_data_central_europe.csv")

combined_data <- read.csv("combined_treecover_loss_data.csv") %>%
  mutate(country = ifelse(country == "Czech Republic", "Czechia", country))

net_cover_data <- read.csv("combined_net_cover_change_data.csv")  %>%
  mutate(country = ifelse(country == "Czech Republic", "Czechia", country))

fire_loss_data <- read.csv("combined_fire_loss_data.csv") %>%
  mutate(country = ifelse(country == "Czech Republic", "Czechia", country))


# total loss data by year for each country
total_loss_summary <- combined_data %>%
  group_by(country, umd_tree_cover_loss__year) %>%
  summarise(total_loss = sum(umd_tree_cover_loss__ha, na.rm = TRUE))

# detailed loss data by year and driver for each country
detailed_loss_summary <- combined_data %>%
  group_by(country, umd_tree_cover_loss__year, tsc_tree_cover_loss_drivers__driver) %>%
  summarise(total_loss = sum(umd_tree_cover_loss__ha, na.rm = TRUE))

# vector of all unique years
years <- unique(total_loss_summary$umd_tree_cover_loss__year)


# pivot longer
time_series_data <- data %>%
  select(shapeName, starts_with("water_percentage"), starts_with("tree_percentage"),
         starts_with("flooded_vegetation_percentage"), starts_with("crops_percentage"),
         starts_with("built_area_percentage"), starts_with("bare_ground_percentage"),
         starts_with("snow_ice_percentage"), starts_with("clouds_percentage"),
         starts_with("rangeland_percentage"), starts_with("burned_area"),
         starts_with("CO2_total"), starts_with("PM25_total"),
         starts_with("TPC_total"), starts_with("NMHC_total"),
         starts_with("OC_total"), starts_with("CH4_total"),
         starts_with("SO2_total"), starts_with("BC_total")) %>%
  pivot_longer(cols = -shapeName,
               names_to = "parameter",
               values_to = "value") %>%
  mutate(year = as.integer(str_extract(parameter, "\\d{4}")),
         parameter = str_remove(parameter, "_\\d{4}")) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  group_by(shapeName, year) %>%
  summarize(across(everything(), mean, na.rm = TRUE), .groups = 'drop')


time_series_data <- data %>%
  select(shapeName, shapeGroup, total_area_km2) %>%
  left_join(time_series_data, by = "shapeName")

time_series_data <- time_series_data %>% mutate(country = case_when(
  shapeGroup == "AUT" ~ "Austria",
  shapeGroup == "CZE" ~ "Czechia",
  shapeGroup == "HUN" ~ "Hungary",
  shapeGroup == "DEU" ~ "Germany",
  shapeGroup == "POL" ~ "Poland",
  shapeGroup == "SVK" ~ "Slovakia",
  shapeGroup == "CHE" ~ "Switzerland",
  shapeGroup == "SVN" ~ "Slovenia",
  TRUE ~ NA_character_
)) 


# state boarders for plotting
central_europe <- rbind(gb_adm1("Austria"), gb_adm1("Germany"), 
                        gb_adm1("Czech Republic"), gb_adm1("Poland"), 
                        gb_adm1("Slovakia"), gb_adm1("Hungary"), 
                        gb_adm1("Switzerland"), gb_adm1("Slovenia"))

data_with_state_boarders <- central_europe %>% 
  select(shapeName, geometry) %>%
  left_join(data, by = "shapeName")

data_with_state_boarders <- data_with_state_boarders %>% mutate(country = case_when(
  shapeGroup == "AUT" ~ "Austria",
  shapeGroup == "CZE" ~ "Czechia",
  shapeGroup == "HUN" ~ "Hungary",
  shapeGroup == "DEU" ~ "Germany",
  shapeGroup == "POL" ~ "Poland",
  shapeGroup == "SVK" ~ "Slovakia",
  shapeGroup == "CHE" ~ "Switzerland",
  shapeGroup == "SVN" ~ "Slovenia",
  TRUE ~ NA_character_
)) 

get_color_palette <- function(variable, domain) {
  switch(variable,
         "tree_percentage_change" = colorBin(palette = c("#ff0000", "#006f00"), domain = domain, bins = 5),  # Green to Yellow
         "built_area_percentage_change" = colorBin(palette = c("#f7de7c", "#ff0000"), domain = domain, bins = 5),  # Shades of Red
         "crops_percentage_change" = colorBin(palette = c("#f7de7c", "#ff7f00"), domain = domain, bins = 5),  # Shades of Orange
         "rangeland_percentage_change" = colorBin(palette = c("#f7de7c", "#8B4513"), domain = domain, bins = 5),  # Shades of Brown
         colorBin(palette = "YlOrRd", domain = domain, bins = 10)  
  )
}




library(shiny)
library(shinyjs)
library(shinyalert)
library(shinyBS)
library(plotly)
library(dplyr)
library(tools)  
library(leaflet)
library(rnaturalearth)
library(sf)
library(reshape2)
library(dtw)
library(dbscan)
library(RColorBrewer)
library(mclust)
library(cluster)
library(ggplot2)
library(reshape2)

# Function to format column names
format_param_names <- function(names) {
  names %>%
    gsub("_", " ", .) %>%  
    toTitleCase()   
}

# Create a mapping of formatted names to original names
param_mapping <- setNames(
  names(time_series_data)[!names(time_series_data) %in% c("shapeName", "shapeGroup", "country", "year", "tree_percentage", "total_area_km2")],
  format_param_names(names(time_series_data)[!names(time_series_data) %in% c("shapeName", "year",  "shapeGroup", "country", "tree_percentage", "total_area_km2")])
)

param_mapping_with_tree_percentage <- setNames(
  names(time_series_data)[!names(time_series_data) %in% c("shapeName", "year", "total_area_km2", "shapeGroup", "country")],
  format_param_names(names(time_series_data)[!names(time_series_data) %in% c("shapeName", "year", "shapeGroup", "country", "total_area_km2")])
)

# Define UI
ui <- fluidPage(
  h2("Central Europe Environmental Data Analysis", align = "center"),
  
  fluidRow(
    column(12, 
           downloadButton("downloadAllData", "Download Data", class = "btn-primary"),
           actionButton("help_general", "Help"),
           hr() 
    )
  ),
  
  tabsetPanel(
    # Time Series Analysis tab
    tabPanel("Time Series Analysis",
             h3("State Statistics and Time Series Analysis"),
             p("Select a state by clicking on the map to analyze the time series data of environmental variables."),
             actionButton("help_time_series", "Help"),
             
             div(style = "position: relative; margin-bottom: 30px;", 
                 leafletOutput("statesMap", height = "700px"),
                 absolutePanel(
                   top = 10, right = 10,
                   style = "background: white; padding: 10px;"
                 )
             ),
             
             div(style = "margin-bottom: 30px;", uiOutput("dynamicSidebarLayout")),
             h3("Trend Line Plot"),
             div(style = "margin-bottom: 30px;", plotlyOutput("linePlotWithTrend")),
             h3("Scatter Plot"),
             div(style = "margin-bottom: 30px;", uiOutput("dynamicScatterPlot")),
             
             fluidRow(
               column(6, h3("Correlation Matrix"), 
                      div(style = "margin-bottom: 30px;", 
                          plotlyOutput("correlationMatrix", height = "600px"))),
               column(6, h3("Heatmap of Variables"),
                      div(style = "margin-bottom: 30px;",
                          plotlyOutput("heatmapPlot", height = "600px")))
             )
    ),
    
    tabPanel("Country Analysis",
             h3("Tree Loss Across Central Europe"),
             p("Select a country to view tree cover loss statistics over time, categorized by different drivers."),
             actionButton("help_country_analysis", "Help"),
             
             div(style = "position: relative; margin-bottom: 30px;",
                 leafletOutput("treeLossMap", height = "700px"),
                 absolutePanel(
                   top = 10, right = 10,
                   style = "background: white; padding: 10px;",
                   selectInput("year", "Choose a year:", 
                               sort(unique(combined_data$umd_tree_cover_loss__year)),
                               selected = 2023)
                 )
             ),
             
             h3("Total Tree Loss by Year"),
             fluidRow(
               div(style = "margin-bottom: 30px;", plotlyOutput("totalLossPlot", height = "600px"))
             ),
             h3("Tree Loss by Driver"),
             fluidRow(
               div(style = "margin-bottom: 30px;", plotlyOutput("lossByDriverPlot", height = "600px"))
             ),
             h3("Fire-Related Tree Loss"),
             fluidRow(
               div(style = "margin-bottom: 30px;", plotlyOutput("lossByFirePlot", height = "600px")),
               div(style = "margin-bottom: 30px;", htmlOutput("fireLossSummary"))
             ),
             h3("Net Tree Cover Change"),
             fluidRow(
               div(style = "margin-bottom: 30px;", plotlyOutput("donutChart", height = "600px"))
             ),
             h3("Additional Statistics"),
             div(style = "margin-bottom: 30px;", uiOutput("reactiveHistogram")),
             div(style = "margin-bottom: 30px;", uiOutput("compareChange")),
             
             h3("Radar Plot of Environmental Factors"),
             div(style = "position: relative; margin-bottom: 30px;",
                 plotlyOutput("radarPlot", height = "700px"),
                 absolutePanel(
                   top = 10, left = 10,
                   style = "background: transparent; padding: 10px;",
                   selectInput("year", "Choose a year:", 
                               2018:2023,
                               selected = 2023)
                 ))
    ),
    
    tabPanel("Central Europe Statistics",
             h3("Area-Based Statistics for Central Europe"),
             p("Explore various environmental statistics over the Central European region. Use clustering to group regions based on selected variables."),
             actionButton("help_statistics", "Help"),
             
             div(style = "margin-bottom: 30px;", plotlyOutput("areaPlotCentralEurope", height = "600px")),
             h3("Bubble Chart"),
             div(style = "margin-bottom: 30px;", uiOutput("dynamicBubbleChart")),
             h3("Tree Loss Treemap"),
             div(style = "margin-bottom: 30px;", plotlyOutput("treeLossTreemap", height = "600px")),
             h3("Feature Change Map"),
             div(style = "margin-bottom: 30px;", uiOutput("featureChangeMap")),
             
             
             h3("Clustering Plots"),
             div(style = "margin-bottom: 30px;", uiOutput("dynamicClusters"))
             
             
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  observeEvent(input$help_general, {
    shinyalert(
      title = "General Help",
      text = "This app provides an interactive way to analyze spatio-temporal sustainability data across Central Europe. You can switch between tabs for time series analysis, country analysis, and detailed statistics.",
      type = "info"
    )
  })
  
  observeEvent(input$help_time_series, {
    shinyalert(
      title = "Time Series Analysis Help",
      text = "In this section, you can select a state by clicking on the map and analyze the time series data for different parameters. You can also generate scatter plots, correlation matrices, and more.",
      type = "info"
    )
  })
  
  observeEvent(input$help_country_analysis, {
    shinyalert(
      title = "Country Analysis Help",
      text = "Here, you can explore tree cover loss statistics across countries by year. Click on a country to show statistics such as total loss, loss by driver, and loss from fires.",
      type = "info"
    )
  })
  
  observeEvent(input$help_statistics, {
    shinyalert(
      title = "Central Europe Statistics Help",
      text = "This section provides overall statistics and visualizations across the Central European region. You can explore area changes, clustering plots, and compare state-wise changes over time.",
      type = "info"
    )
  })
  
  
  output$downloadAllData <- downloadHandler(
    filename = function() {
      paste("sustainability-datasets-", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      temp_dir <- tempdir()
      
      write.csv(data, file.path(temp_dir, "sustainability_data_central_europe.csv"), row.names = FALSE)
      write.csv(combined_data, file.path(temp_dir, "combined_treecover_loss_data.csv"), row.names = FALSE)
      write.csv(net_cover_data, file.path(temp_dir, "combined_net_cover_change_data.csv"), row.names = FALSE)
      write.csv(fire_loss_data, file.path(temp_dir, "combined_fire_loss_data.csv"), row.names = FALSE)
      write.csv(time_series_data, file.path(temp_dir, "time_series_converted_data_set.csv"), row.names = FALSE)
      
      zip::zipr(file, files = list.files(temp_dir, full.names = TRUE))
    },
    contentType = "application/zip"
  )
  
  clicked_state <- reactiveVal(NULL)
  
  output$areaPlotCentralEurope <- renderPlotly({
    cumulative_data <- combined_data %>%
      group_by(umd_tree_cover_loss__year) %>%
      summarise(cumulative_loss = cumsum(sum(umd_tree_cover_loss__ha, na.rm = TRUE)))
    
    plot_ly(cumulative_data, x = ~umd_tree_cover_loss__year, y = ~cumulative_loss, type = 'scatter', mode = 'lines',
            fill = 'tozeroy', name = 'Cumulative Tree Loss',
            hovertemplate = paste(
              "<b>Year:</b> %{x}<br>",
              "<b>Cumulative Loss:</b> %{y:.2f} ha",
              "<extra></extra>"  
            )) %>%
      layout(title = "Cumulative Tree Cover Loss Over Time",
             xaxis = list(title = "Year"),
             yaxis = list(title = "Cumulative Loss (ha)"))
  })
  
  output$dynamicBubbleChart <- renderUI({
    fluidRow(
      column(2,
             selectInput("xVar", "Select X Variable:",
                         choices = param_mapping,
                         selected = "CO2_total"),
             selectInput("yVar", "Select Y Variable:",
                         choices = param_mapping,
                         selected = "PM25_total")
      ),
      column(10,  plotlyOutput("bubbleChart"))
    )
  })
  
  output$bubbleChart <- renderPlotly({
    req(input$xVar, input$yVar)
    
    x_var_name <- names(param_mapping)[which(param_mapping == input$xVar)]
    y_var_name <- names(param_mapping)[which(param_mapping == input$yVar)]
    
    plot_ly(time_series_data, x = ~get(input$xVar), y = ~get(input$yVar), type = 'scatter', mode = 'markers',
            size = ~tree_percentage, color = ~country, 
            text = ~paste("Country:", country, "<br>",
                          x_var_name, ": ", get(input$xVar), "<br>",
                          y_var_name, ": ", get(input$yVar), "<br>",
                          "Tree Percentage: ", round(tree_percentage, 2), "%"),  
            hovertemplate = paste(
              "<b>%{text}</b><extra></extra>"  
            )) %>%
      layout(title = paste0(x_var_name, " vs. ", y_var_name, " with Tree Percentage Bubble Size"),
             xaxis = list(title = x_var_name),
             yaxis = list(title = y_var_name))
  })
  
  output$treeLossTreemap <- renderPlotly({
    plot_ly(data = detailed_loss_summary, labels = ~tsc_tree_cover_loss_drivers__driver, parents = "", values = ~total_loss, type = 'treemap',
            textinfo = "label+value+percent entry") %>%
      layout(title = "Tree Cover Loss by Driver")
  })
  
  output$radarPlot <- renderPlotly({
    req(clicked_country())
    
    radar_data <- filter(time_series_data, year == input$year & country == clicked_country())
    
    radar_aggregated <- radar_data %>%
      summarise(tree_percentage = mean(tree_percentage, na.rm = TRUE),
                CO2_total = mean(CO2_total, na.rm = TRUE),
                PM25_total = mean(PM25_total, na.rm = TRUE),
                flooded_vegetation_percentage = mean(flooded_vegetation_percentage, na.rm = TRUE),
                built_area_percentage = mean(built_area_percentage, na.rm = TRUE),
                crops_percentage = mean(crops_percentage, na.rm = TRUE))
    
    plot_ly(type = 'scatterpolar', 
            r = c(radar_aggregated$tree_percentage, radar_aggregated$CO2_total, radar_aggregated$PM25_total,
                  radar_aggregated$flooded_vegetation_percentage, radar_aggregated$built_area_percentage,
                  radar_aggregated$crops_percentage),
            theta = c('Tree Cover Percentage', 'CO2 Emissions', 'PM2.5 Levels', 'Flooded Vegetation Percentage',
                      'Built Area Percentage', 'Crops Percentage'), 
            fill = 'toself') %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, max(radar_aggregated, na.rm = TRUE)))),
             title = paste("Environmental Factor Comparison for", clicked_country(), "in", input$year))
  })
  
  output$statesMap <- renderLeaflet({
    
    colors <- brewer.pal(nrow(data_with_state_boarders), "Set3")
    palette <- colorFactor(palette = colors, domain = data_with_state_boarders$shapeName)
    
    leaflet(data_with_state_boarders) %>%
      addTiles() %>%
      addPolygons(
        fillColor = ~palette(shapeName),
        fillOpacity = 0.7,
        color = "black",
        weight = 1,
        highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.7, bringToFront = TRUE),
        label = ~paste("Show statistics for ", shapeName),
        layerId = ~shapeName  
      ) %>%
      setView(lng = 15, lat = 48, zoom = 4) %>%
      addControl("<strong>Central Europe Areas Map</strong>", 
                 position = "bottomleft", 
                 className = "map-title")
  })
  
  observeEvent(input$help_plot1, {
    shinyalert("How to Use States Map", "This plot shows map of states in Central europe. To display statistiscs for a specific state, click on it.")
  })
  
  observeEvent(input$statesMap_shape_click, {
    clicked_state(input$statesMap_shape_click$id)
  })
  
  output$dynamicClusters <- renderUI({
    sidebarLayout(
      sidebarPanel(
        selectInput("year", "Choose year to plot clusters:", 
                    choices = 2018:2023, 
                    selected = 2023),
        selectInput("clusteringMethod", "Select Clustering Method:", 
                    choices = c("KMeans", "GMM", "ST-DBSCAN", "DTW"))
      ),
      mainPanel(
        plotlyOutput("clusterPlot", height = "600px")
      )
    )
  })
  
  
  selected_year_data <- reactive({
    time_series_data %>% filter(year == input$year)
  })
  
  cluster_data <- reactive({
    tree_data <- selected_year_data() %>%
      select(tree_percentage) %>%
      as.matrix()
    
    dtw_dist <- dist(tree_data, method = "DTW")
    
    # Hierarchical clustering based on DTW distances
    hc_dtw <- hclust(as.dist(dtw_dist), method = "average")
    
    selected_year_data() %>%
      mutate(cluster_dtw = cutree(hc_dtw, k = 4))
  })
  
  
  
  output$clusterPlot <- renderPlotly({
    selected_data <- selected_year_data()
    
    method <- input$clusteringMethod
    if (method == "DTW") {
      # DTW Clustering
      tree_data <- selected_data %>%
        select(tree_percentage) %>%
        as.matrix()
      
      dtw_dist <- dist(tree_data, method = "DTW")
      
      # Hierarchical clustering based on DTW distances
      hc_dtw <- hclust(as.dist(dtw_dist), method = "average")
      
      # Add DTW clusters to the data
      cluster_df <- selected_data %>%
        mutate(cluster = as.factor(cutree(hc_dtw, k = 4)))
      
      plot_ly(cluster_df, 
              x = ~total_area_km2, 
              y = ~tree_percentage, 
              type = 'scatter', 
              mode = 'markers', 
              color = ~cluster, 
              marker = list(size = 10), 
              text = ~paste("State:", shapeName, "<br>",
                            "Cluster:", cluster, "<br>",
                            "Tree Percentage:", round(tree_percentage, 2), "%<br>",
                            "Total Area:", round(total_area_km2, 2), "km²"), 
              hoverinfo = "text") %>%
        layout(title = paste("DTW Clustering of Tree Percentage in", input$year),
               xaxis = list(title = "Total Area (km²)"),
               yaxis = list(title = "Tree Percentage (%)"))
    } else {
      tree_data <- selected_data %>%
        select(tree_percentage) %>%
        as.matrix()
      
      cluster_result <- switch(method,
                               "KMeans" = kmeans(tree_data, centers = 4),
                               "GMM" = Mclust(tree_data, G = 4),
                               "ST-DBSCAN" = dbscan(tree_data, eps = 2, minPts = 5)
      )
      
      cluster_assignment <- switch(method,
                                   "KMeans" = cluster_result$cluster,
                                   "GMM" = cluster_result$classification,
                                   "ST-DBSCAN" = cluster_result$cluster)
      
      cluster_df <- selected_data %>%
        mutate(cluster = as.factor(cluster_assignment))
      
      plot_ly(cluster_df, 
              x = ~total_area_km2, 
              y = ~tree_percentage, 
              type = 'scatter', 
              mode = 'markers', 
              color = ~cluster, 
              marker = list(size = 10), 
              text = ~paste("State:", shapeName, "<br>",
                            "Cluster:", cluster, "<br>",
                            "Tree Percentage:", round(tree_percentage, 2), "%<br>",
                            "Total Area:", round(total_area_km2, 2), "km²"), 
              hoverinfo = "text") %>%
        layout(title = paste(method, "Clustering of Tree Percentage in", input$year),
               xaxis = list(title = "Total Area (km²)"),
               yaxis = list(title = "Tree Percentage (%)"))
    }
  })
  
 
  output$dynamicSidebarLayout <- renderUI({
    req(clicked_state())  
    
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("parameters", "Choose parameters to plot:", 
                           choices = names(param_mapping), 
                           selected = c())
      ),
      mainPanel(
        plotlyOutput("timeSeriesPlot", height = "600px")  
      )
    )
  })
  
  output$dynamicScatterPlot <- renderUI({
    req(clicked_state())
    fluidRow(
      column(2,
             selectInput("xVar", "Select X Variable:",
                         choices = param_mapping_with_tree_percentage,
                         selected = "tree_percentage"),
             selectInput("yVar", "Select Y Variable:",
                         choices = param_mapping_with_tree_percentage,
                         selected = "crops_percentage")
      ),
      column(10,  plotlyOutput("scatterMatrixPlot"))
    )
  })
  
  clicked_country <- reactiveVal(NULL)
  
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  
  
  filtered_data <- reactive({
    combined_data %>%
      filter(umd_tree_cover_loss__year == input$year) %>%
      group_by(country) %>%
      summarise(tree_loss = sum(umd_tree_cover_loss__ha, na.rm = TRUE))
  })
  
  central_europe_countries <- c("Austria", "Czechia", "Germany", "Hungary", 
                                "Poland", "Slovakia", "Slovenia", "Switzerland")
  
  
  geo_data <- reactive({
    countries %>%
      left_join(filtered_data(), by = c("name" = "country")) %>%
      filter(name %in% central_europe_countries) 
  })
  
  
  output$treeLossMap <- renderLeaflet({
    leaflet(data = geo_data()) %>%
      addTiles() %>%
      addPolygons(
        fillColor = ~colorBin("YlOrRd", geo_data()$tree_loss)(tree_loss),
        fillOpacity = 0.7,
        color = "black",
        weight = 1,
        highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.7, bringToFront = TRUE),
        label = ~paste(name, ", Tree Loss:", round(tree_loss, 2), "ha"),
        layerId = ~name 
      ) %>%
      setView(lng = 15, lat = 48, zoom = 4) %>%
      addControl("<strong>Tree Cover Loss Across Central Europe by Year</strong>", 
                 position = "bottomleft", 
                 className = "map-title")
  })
  
  observeEvent(input$treeLossMap_shape_click, {
    clicked_country(input$treeLossMap_shape_click$id)  
  })
  
  output$totalLossPlot <- renderPlotly({
    req(clicked_country())  
    
    country_total_loss <- filter(total_loss_summary, country == clicked_country())
    
    plot_ly(country_total_loss, x = ~umd_tree_cover_loss__year, y = ~total_loss, type = 'bar',
            color = "pink", 
            text = ~paste("Year:", umd_tree_cover_loss__year, "<br>Loss (ha):", round(total_loss, 2)),
            hoverinfo = 'text') %>%
      layout(title = paste("Total Tree Cover Loss by Year in", clicked_country()),
             xaxis = list(title = "Year", 
                          tickvals = years, 
                          tickmode = 'array', 
                          tickangle = -45), 
             yaxis = list(title = "Tree Cover Loss (ha)"),
             plot_bgcolor = 'white',  
             margin = list(l = 60, r = 60, t = 60, b = 60))
  })
  
  output$lossByDriverPlot <- renderPlotly({
    req(clicked_country()) 
    
    country_detailed_loss <- filter(detailed_loss_summary, country == clicked_country())
    
    plot_ly(country_detailed_loss, x = ~umd_tree_cover_loss__year, y = ~total_loss, type = 'bar',
            color = ~tsc_tree_cover_loss_drivers__driver,  
            text = ~paste("Driver:", tsc_tree_cover_loss_drivers__driver, "<br>Year:", umd_tree_cover_loss__year, "<br>Loss (ha):", round(total_loss, 2)),
            hoverinfo = 'text') %>%
      layout(barmode = 'stack',  
             title = paste("Tree Cover Loss by Year and Driver in", clicked_country()),
             xaxis = list(title = "Year", 
                          tickvals = years,  
                          tickmode = 'array', 
                          tickangle = -45),  
             yaxis = list(title = "Tree Cover Loss (ha)"),
             plot_bgcolor = 'white', 
             margin = list(l = 60, r = 60, t = 60, b = 60))  
  })
  
  output$lossByFirePlot <- renderPlotly({
    req(clicked_country())  
    
    country_data <- filter(fire_loss_data, country == clicked_country())
    
    plot_ly(country_data, x = ~umd_tree_cover_loss__year, y = ~umd_tree_cover_loss_from_fires__ha, type = 'bar',
            marker = list(color = '#FF6F6F'),
            text = ~paste("Year:", umd_tree_cover_loss__year, "<br>Loss from Fires (ha):", round(umd_tree_cover_loss_from_fires__ha, 2)),
            hoverinfo = 'text') %>%
      layout(title = paste("Tree Cover Loss from Fires by Year in", clicked_country()),
             xaxis = list(title = "Year", 
                          tickvals = years, 
                          tickmode = 'array', 
                          tickangle = -45), 
             yaxis = list(title = "Tree Cover Loss from Fires (ha)"),
             plot_bgcolor = 'white',  
             margin = list(l = 60, r = 60, t = 60, b = 60))
  })
  
  
  output$fireLossSummary <- renderText({
    req(clicked_country())
    
    country_data <- filter(fire_loss_data, country == clicked_country())
    
    total_fire_loss <- sum(country_data$umd_tree_cover_loss_from_fires__ha, na.rm = TRUE)
    
    total_loss_all_drivers <- sum(combined_data$umd_tree_cover_loss__ha[combined_data$country == clicked_country()], na.rm = TRUE)
    
    max_fire_loss_year <- country_data %>%
      filter(umd_tree_cover_loss_from_fires__ha == max(umd_tree_cover_loss_from_fires__ha, na.rm = TRUE)) %>%
      pull(umd_tree_cover_loss__year)
    
    max_fire_loss_amount <- max(country_data$umd_tree_cover_loss_from_fires__ha, na.rm = TRUE)
    max_fire_loss_percentage <- (max_fire_loss_amount / total_loss_all_drivers) * 100
    
    paste0("<div style='text-align: center; font-size: 18px;'>",
           "From 2001 to 2023, ", clicked_country(),
           " lost <b>", round(total_fire_loss, 2), " ha</b> of tree cover from fires and <b>", 
           round(total_loss_all_drivers / 1000, 1), " kha </b> from all other drivers of loss. ",
           "The year with the most tree cover loss due to fires during this period was <b>", max_fire_loss_year,
           "</b> with <b>", round(max_fire_loss_amount, 2), " ha </b> lost to fires — <b>", 
           round(max_fire_loss_percentage, 1), "%</b> of all tree cover loss for that year.",
           "</div>")
  })
  
  
  output$donutChart <- renderPlotly({
    req(clicked_country())  
    
    country_net_cover <- filter(net_cover_data, country == clicked_country())
    
    if (nrow(country_net_cover) == 0) {
      return(NULL)  
    }
    
    donut_data <- country_net_cover %>%
      select(stable, loss, gain, disturb) %>%
      pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
      mutate(
        metric = case_when(
          metric == "stable" ~ "Stable forest",
          metric == "loss" ~ "Loss",
          metric == "gain" ~ "Gain",
          metric == "disturb" ~ "Disturbed",
          TRUE ~ metric
        )
      )
    
    net_change <- country_net_cover$net / 1000  
    change_percentage <- country_net_cover$change
    
    custom_colors <- c(
      "Stable forest" = "#2ca02c",   
      "Loss" = "#d62728",     
      "Gain" = "#1f77b4",   
      "Disturbed" = "#ff7f0e"     
    )
    
    plot_ly(donut_data, labels = ~metric, values = ~value, type = 'pie', hole = 0.4,
            textinfo = 'label+percent',
            marker = list(colors = custom_colors[donut_data$metric])) %>%
      layout(
        title = paste("From 2000 to 2020, <b>", clicked_country(), 
                      "</b> experienced a net change of <b>", round(net_change, 2), "</b>kha (<b>", round(change_percentage, 2), "%</b>) in tree cover."),
        showlegend = TRUE,
        legend = list(
          title = "Metrics",
          orientation = "v", 
          xanchor = "center",
          yanchor = "bottom"
        )
      )
  })
  
  output$reactiveHistogram <- renderUI({
    req(clicked_country())
    sidebarLayout(
      sidebarPanel(
        selectInput("year", "Choose a year:", 
                    sort(unique(time_series_data$year)),
                    selected = 2023),
        selectInput("variable", "Select Variable:", 
                    choices = param_mapping_with_tree_percentage,
                    selected = "tree_percentage")
      ),
      mainPanel(
        plotlyOutput("histogramPlot", height = "600px")  
      )
    )
  })
  
  output$compareChange <- renderUI({
    req(clicked_country())
    years <- 2018:2023
    sidebarLayout(
      sidebarPanel(
        selectInput("compare_year", "Select Year to Compare with 2018:", choices = years[-1], selected = 2023),
        selectInput("variableStateChange", "Select Change to Display:", choices = 
                      c("tree_percentage_change", "built_area_percentage_change",
                        "crops_percentage_change", "rangeland_percentage_change" ))
      ),
      mainPanel(
        h4(paste0("Display Changes in Areas of ", clicked_country())),
        h5("Click state of interest to display concrete statistics"),
        leafletOutput("stateMapPlot", height = "800px")  
      )
    )
  })
  
  data_for_anual_state_change <- reactive({
    req(clicked_country())
    data_with_state_boarders %>%
      filter(country == clicked_country())
  })
  
  annual_change_df <- reactive({
    country_df <- data_for_anual_state_change()
    
    get_column_name <- sub("_[^_]*$", "", input$variableStateChange)
    
    col_current <- paste0(get_column_name, "_", input$compare_year)
    col_base <- paste0(get_column_name, "_2018")
    col_change <- paste0(get_column_name, "_change_2018_", input$compare_year)
    
    if (!all(c(col_current, col_base) %in% colnames(country_df))) {
      stop(paste("Selected variable", get_column_name, "is not available for the chosen years."))
    }
    
    country_df %>%
      mutate(
        !!col_change := (get(col_current) - get(col_base))
      )
  })
  
  output$stateMapPlot <- renderLeaflet({
    plot_data <- annual_change_df()
    get_column_name <- sub("_[^_]*$", "", input$variableStateChange)
    col_change <- paste0(get_column_name, "_change_2018_", input$compare_year)
    
    if (!col_change %in% colnames(plot_data)) {
      stop(paste("Change data for", input$variableStateChange, "is not available."))
    }
    
    plot_data_sf <- st_as_sf(plot_data)
    
    domain <- range(plot_data_sf[[col_change]], na.rm = TRUE)
    
    color_bins <- get_color_palette(input$variableStateChange, domain)
    
    leaflet(data = plot_data_sf) %>%
      addTiles() %>%
      addPolygons(
        fillColor = ~color_bins(get(col_change)),
        fillOpacity = 0.7,
        color = "black",
        weight = 1,
        highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.7, bringToFront = TRUE),
        popup = ~paste0(
          "<div style='font-size: 14px;'><b>State:</b> ", shapeName, "<br>",
          "<b>", input$variableStateChange, " Change:</b> ", round(get(col_change), 2), "%</div>"
        ),
        layerId = ~shapeName
      ) %>%
      addLegend(
        pal = color_bins,
        values = ~get(col_change),
        title = paste(input$variableStateChange, "Change from 2018 to", input$compare_year),
        position = "bottomright",
        labFormat = labelFormat(suffix = "%")  
      ) %>%
      setView(lng = 15, lat = 48, zoom = 4) %>%
      addControl("<strong>Central Europe Areas Map</strong>", 
                 position = "bottomleft", 
                 className = "map-title")
  })
  
  
  
  filtered_year_data <- reactive({
    req(input$year, clicked_country())
    filter(time_series_data, country == clicked_country(), year == input$year)
  })
  
  output$histogramPlot <- renderPlotly({
    country_data <- filtered_year_data()
    
    if (nrow(country_data) == 0 || !input$variable %in% names(country_data)) {
      return(NULL)
    }
    
    var_name <- names(param_mapping_with_tree_percentage)[which(param_mapping_with_tree_percentage == input$variable)]
    
    plot_ly(country_data, 
            x = ~get(input$variable), 
            type = 'histogram', 
            marker = list(color = '#a4ebe0', line = list(color = '#a39d9d', width = 2)),
            hovertemplate = paste(
              'Bin range: %{x}', 
              '<br>Count: %{y}', 
              '<extra></extra>'   
            )) %>%
      layout(title = list(text = paste("Histogram of", var_name, "for", clicked_country(), "in Year", input$year),
                          font = list(size = 18, family = 'Arial Black')),
             xaxis = list(title = list(text = var_name, font = list(size = 14, family = 'Arial')),
                          tickfont = list(size = 12, family = 'Arial')),
             yaxis = list(title = list(text = 'Frequency', font = list(size = 14, family = 'Arial')),
                          tickfont = list(size = 12, family = 'Arial')),
             showlegend = FALSE,
             hovermode = 'closest')
    
  })
  
  output$timeSeriesPlot <- renderPlotly({
    req(clicked_state())
    
    state_data <- filter(time_series_data, shapeName == clicked_state())
    
    p <- plot_ly(state_data, x = ~year, y = ~tree_percentage, type = 'scatter', mode = 'lines+markers',
                 name = 'Tree Percentage', 
                 text = paste('Year:', state_data$year, '<br>Tree Percentage:', round(state_data$tree_percentage, 2)),
                 hoverinfo = 'text') 
    
    if (length(input$parameters) > 0) {
      for (param in input$parameters) {
        original_param <- param_mapping[param]
        
        p <- p %>% add_trace(y = state_data[[original_param]], mode = 'lines+markers',
                             name = param, 
                             text = paste('Year:', state_data$year, '<br>', param, ':', round(state_data[[original_param]], 2)),
                             hoverinfo = 'text')
      }
    }
    
    p <- add_annotations(p,
                         x = max(state_data$year),
                         y = max(state_data$tree_percentage),
                         text = paste("Total Area of ", clicked_state(), ": ", round(state_data$total_area_km2[1], 2), "km²")
    )
    
    p <- layout(p, title = paste("Time Series for", clicked_state()),
                xaxis = list(title = 'Year'),
                yaxis = list(title = 'Value'))
    
    p
  })
  
  output$correlationMatrix <- renderPlotly({
    req(clicked_state())
    
    state_data <- filter(time_series_data, shapeName == clicked_state())
    
    numeric_data <- state_data %>%
      select(tree_percentage, crops_percentage, burned_area, water_percentage,
             flooded_vegetation_percentage, built_area_percentage,
             rangeland_percentage, 
             CO2_total, PM25_total, TPC_total, NMHC_total,
             OC_total, CH4_total, SO2_total, BC_total ) %>%
      na.omit()
    
    correlation_matrix <- cor(numeric_data)
    
    melted_correlation_matrix <- melt(correlation_matrix)
    
    plot_ly(melted_correlation_matrix, x = ~Var1, y = ~Var2, z = ~value, type = 'heatmap', 
            colorscale = 'Viridis') %>%
      layout(title = paste("Correlation Matrix for", clicked_state()),
             xaxis = list(title = 'Parameters'),
             yaxis = list(title = 'Parameters'))
  })
  
  output$heatmapPlot <- renderPlotly({
    req(clicked_state())
    
    state_data <- filter(time_series_data, shapeName == clicked_state())
    
    heatmap_data <- state_data %>%
      select(year, tree_percentage) %>%
      spread(key = year, value = tree_percentage) %>%
      as.matrix()
    
    plot_ly(z = heatmap_data, type = "heatmap", colorscale = "Viridis") %>%
      layout(title = paste("Heatmap of Tree Percentage by Year for", clicked_state()),
             xaxis = list(title = 'Years'),
             yaxis = list(title = 'Tree Percentage'))
  })
  
  output$linePlotWithTrend <- renderPlotly({
    req(clicked_state())
    
    state_data <- filter(time_series_data, shapeName == clicked_state())
    
    p <- plot_ly(state_data, 
                 x = ~year, 
                 y = ~tree_percentage, 
                 type = 'scatter', 
                 mode = 'lines+markers',
                 name = 'Tree Percentage', 
                 line = list(color = 'blue'),
                 text = ~paste("Year: ", year, "<br>",
                               "Tree Percentage: ", round(tree_percentage, 2), "%<br>",
                               "State: ", clicked_state()),
                 hoverinfo = "text")
    
    p <- p %>% add_trace(y = ~fitted(lm(tree_percentage ~ year, data = state_data)),
                         type = 'scatter', 
                         mode = 'lines',
                         name = 'Trend Line', 
                         line = list(color = 'red', dash = 'dash'),
                         text = ~paste("Year: ", year, "<br>",
                                       "Trend Line Value: ", round(fitted(lm(tree_percentage ~ year, data = state_data)), 2), "%<br>",
                                       "State: ", clicked_state()),
                         hoverinfo = "text")
    
    p <- layout(p, 
                title = paste("Tree Percentage Trend for", clicked_state()),
                xaxis = list(title = 'Year'),
                yaxis = list(title = 'Tree Percentage (%)'))
    
    return(p)
  })
  
  output$scatterMatrixPlot <- renderPlotly({
    req(clicked_state())
    
    state_data <- filter(time_series_data, shapeName == clicked_state())
    
    x_var_name <- names(param_mapping_with_tree_percentage)[which(param_mapping_with_tree_percentage == input$xVar)]
    y_var_name <- names(param_mapping_with_tree_percentage)[which(param_mapping_with_tree_percentage == input$yVar)]
    
    plot_ly(data = state_data, x = ~get(input$xVar), y = ~get(input$yVar),
            type = 'scatter', mode = 'markers',
            marker = list(color = 'blue'),
            text = ~paste("Year:", year, "<br>", x_var_name, ":", round(get(input$xVar), 2), "<br>", 
                          y_var_name, ":", round(get(input$yVar), 2)),
            hoverinfo = 'text') %>%
      layout(title = paste("Scatter Plot Matrix for", clicked_state()),
             xaxis = list(title = x_var_name),
             yaxis = list(title = y_var_name))
  })
  
 
  output$featureChangeMap <-renderUI({
    years <- 2018:2023
    sidebarLayout(
      sidebarPanel(
        selectInput("compare_year", "Select Year to Compare with 2018:", choices = years[-1], selected = 2023),
        selectInput("variableStateChange", "Select Change to Display:", choices = 
                      c("tree_percentage_change", "built_area_percentage_change",
                        "crops_percentage_change", "rangeland_percentage_change" ))
      ),
      mainPanel(
        h4("Display Changes in Areas of Central Europe"),
        h5("Click state of interest to display concrete statistics"),
        leafletOutput("centralEuropChangePlot", height = "800px")  
      )
    )
  })
  

    central_europe_anual_change <- reactive({
    df <- data_with_state_boarders
    
    get_column_name <- sub("_[^_]*$", "", input$variableStateChange)
    
    col_current <- paste0(get_column_name, "_", input$compare_year)
    col_base <- paste0(get_column_name, "_2018")
    col_change <- paste0(get_column_name, "_change_2018_", input$compare_year)
    
    df %>%
      mutate(
        !!col_change := (get(col_current) - get(col_base))
      )
  })
  
  output$centralEuropChangePlot <- renderLeaflet({
    plot_data <- central_europe_anual_change()
    get_column_name <- sub("_[^_]*$", "", input$variableStateChange)
    col_change <- paste0(get_column_name, "_change_2018_", input$compare_year)
    
    if (!col_change %in% colnames(plot_data)) {
      stop(paste("Change data for", input$variableStateChange, "is not available."))
    }
    
    plot_data_sf <- st_as_sf(plot_data)
    
    domain <- range(plot_data_sf[[col_change]], na.rm = TRUE)
    
    color_bins <- get_color_palette(input$variableStateChange, domain)
    
    leaflet(data = plot_data_sf) %>%
      addTiles() %>%
      addPolygons(
        fillColor = ~color_bins(get(col_change)),
        fillOpacity = 0.7,
        color = "black",
        weight = 1,
        highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.7, bringToFront = TRUE),
        popup = ~paste0(
          "<div style='font-size: 14px;'><b>State:</b> ", shapeName, "<br>",
          "<b>", input$variableStateChange, " Change:</b> ", round(get(col_change), 2), "%</div>"
        ),
        layerId = ~shapeName
      ) %>%
      addLegend(
        pal = color_bins,
        values = ~get(col_change),
        title = paste(input$variableStateChange, "Change from 2018 to", input$compare_year),
        position = "bottomright",
        labFormat = labelFormat(suffix = "%")  
      ) %>%
      setView(lng = 15, lat = 48, zoom = 5) %>%
      addControl("<strong>Central Europe Areas Map</strong>", 
                 position = "bottomleft", 
                 className = "map-title")
  })

}

shinyApp(ui, server)
