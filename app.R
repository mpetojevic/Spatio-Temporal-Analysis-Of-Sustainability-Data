library(dplyr)
library(tidyr)
library(stringr)

getwd()  # Check the current working directory
setwd("~/Uni/2024S/Bachelor/Spatio-Temporal-Analysis-Of-Sustainability-Data/Cluster-Analysis") 
getwd()

# Load your data
data <- read.csv("sustainability_data_central_europe.csv")


combined_data <- read.csv("combined_treecover_loss_data.csv") %>%
  mutate(country = ifelse(country == "Czech Republic", "Czechia", country))

net_cover_data <- read.csv("combined_net_cover_change_data.csv")  %>%
  mutate(country = ifelse(country == "Czech Republic", "Czechia", country))

fire_loss_data <- read.csv("combined_fire_loss_data.csv") %>%
  mutate(country = ifelse(country == "Czech Republic", "Czechia", country))


# Prepare total loss data by year for each country
total_loss_summary <- combined_data %>%
  group_by(country, umd_tree_cover_loss__year) %>%
  summarise(total_loss = sum(umd_tree_cover_loss__ha, na.rm = TRUE))

# Prepare detailed loss data by year and driver for each country
detailed_loss_summary <- combined_data %>%
  group_by(country, umd_tree_cover_loss__year, tsc_tree_cover_loss_drivers__driver) %>%
  summarise(total_loss = sum(umd_tree_cover_loss__ha, na.rm = TRUE))

# Create a vector of all unique years
years <- unique(total_loss_summary$umd_tree_cover_loss__year)


# Prepare time series data
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

library(rgeoboundaries)
central_europe <- rbind(gb_adm1("Austria"), gb_adm1("Germany"), 
                        gb_adm1("Czech Republic"), gb_adm1("Poland"), 
                        gb_adm1("Slovakia"), gb_adm1("Hungary"), 
                        gb_adm1("Switzerland"), gb_adm1("Slovenia"))

time_series_data <- central_europe %>%
  select(shapeName, geometry) %>%
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

# Check the result
head(time_series_data)

library(shiny)
library(plotly)
library(dplyr)
library(tools)  
library(leaflet)
library(rnaturalearth)
library(sf)
library(bslib)

# Function to format column names
format_param_names <- function(names) {
  names %>%
    gsub("_", " ", .) %>%  # Replace underscores with spaces
    toTitleCase()   # Capitalize each word
}

# Create a mapping of formatted names to original names
param_mapping <- setNames(
  names(time_series_data)[!names(time_series_data) %in% c("shapeName", "year", "tree_percentage", "country", "shapeGroup", "total_area_km2" )],
  format_param_names(names(time_series_data)[!names(time_series_data) %in% c("shapeName", "year", "tree_percentage", "country", "shapeGroup", "total_area_km2")])
)

# Define UI
ui <- fluidPage(
  titlePanel("Interactive Time Series and Country Analysis"),
  tabsetPanel(
    tabPanel("Country Analysis",
             h3("Tree Loss across Central Europe by Year"),
             div(style = "position: relative;",
                 leafletOutput("treeLossMap", height = "700px"),
                 absolutePanel(
                   top = 10, right = 10,  # Adjust position as needed
                   style = "background: white; padding: 10px;",
                   selectInput("year", "Choose a year:", 
                               sort(unique(combined_data$umd_tree_cover_loss__year)),
                               selected = 2023)
                 )
             ),
             fluidRow(
               plotlyOutput("totalLossPlot", height = "600px")
             ),
             fluidRow(
               plotlyOutput("lossByDriverPlot", height = "600px")
             ),
             fluidRow(
               plotlyOutput("lossByFirePlot", height = "600px"),
               htmlOutput("fireLossSummary")  
             ),
             fluidRow(plotlyOutput("donutChart", height = "600px")),
             h3("Time Series Comparison Of Tree Cover Change And Other Factors"),
             uiOutput("stateSelectorUI"),
             uiOutput("timeSeriesPlotUI")
    )
  )
)


server <- function(input, output, session) {
  
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
      filter(name %in% central_europe_countries)  # Only Central European countries are kept
  })
  
  
  # Leaflet map for tree loss
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
        layerId = ~name  # Use name as the layer ID to capture clicks
      ) %>%
      setView(lng = 15, lat = 48, zoom = 4) %>%
      addControl("<strong>Tree Cover Loss Across Central Europe by Year</strong>", 
                 position = "bottomleft", 
                 className = "map-title")
  })
  
  observeEvent(input$treeLossMap_shape_click, {
    clicked_country(input$treeLossMap_shape_click$id)  
    updateSelectInput(session, "state", selected = unique(filter(country_states(), country == clicked_country())$shapeName)[1])
  })
  
  output$totalLossPlot <- renderPlotly({
    req(clicked_country())  # Only proceed if a country is clicked
    
    # Filter data for selected country
    country_total_loss <- filter(total_loss_summary, country == clicked_country())
    
    # Total Loss Plot
    plot_ly(country_total_loss, x = ~umd_tree_cover_loss__year, y = ~total_loss, type = 'bar',
            color = "pink",  # Change color to pink
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
    req(clicked_country())  # Only proceed if a country is clicked
    
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
    req(clicked_country())  # Only proceed if a country is clicked
    
    # Filter data for selected country
    country_data <- filter(fire_loss_data, country == clicked_country())
    
    # Fire Loss Plot
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
    
    # Filter data for selected country
    country_data <- filter(fire_loss_data, country == clicked_country())
    
    # Calculate total loss from fires
    total_fire_loss <- sum(country_data$umd_tree_cover_loss_from_fires__ha, na.rm = TRUE)
    
    # Calculate total loss from all drivers
    total_loss_all_drivers <- sum(combined_data$umd_tree_cover_loss__ha[combined_data$country == clicked_country()], na.rm = TRUE)
    
    # Find the year with the most loss from fires
    max_fire_loss_year <- country_data %>%
      filter(umd_tree_cover_loss_from_fires__ha == max(umd_tree_cover_loss_from_fires__ha, na.rm = TRUE)) %>%
      pull(umd_tree_cover_loss__year)
    
    max_fire_loss_amount <- max(country_data$umd_tree_cover_loss_from_fires__ha, na.rm = TRUE)
    max_fire_loss_percentage <- (max_fire_loss_amount / total_loss_all_drivers) * 100
    
    # Text summary
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
    req(clicked_country())  # Only proceed if a country is clicked
    
    country_net_cover <- filter(net_cover_data, country == clicked_country())
    
    if (nrow(country_net_cover) == 0) {
      return(NULL)  # No data for the selected country
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
    
    net_change <- country_net_cover$net / 1000  # Convert to kha
    change_percentage <- country_net_cover$change
    
    # Custom colors for the donut chart slices
    custom_colors <- c(
      "Stable forest" = "#2ca02c",   
      "Loss" = "#d62728",     
      "Gain" = "#1f77b4",   
      "Disturbed" = "#ff7f0e"     
    )
    
    # Donut chart
    plot_ly(donut_data, labels = ~metric, values = ~value, type = 'pie', hole = 0.4,
            textinfo = 'label+percent',
            marker = list(colors = custom_colors[donut_data$metric])) %>%
      layout(
        title = paste("From 2000 to 2020, <b>", clicked_country(), 
                      "</b> experienced a net change of <b>", round(net_change, 2), "</b>kha (<b>", round(change_percentage, 2), "%</b>) in tree cover."),
        showlegend = TRUE,
        legend = list(
          title = "Metrics",
          orientation = "v",  # Horizontal legend
          xanchor = "center",
          yanchor = "bottom"
        )
      )
  })
  
  country_states <- reactive({
    req(clicked_country())  
    filter(time_series_data, country == clicked_country())
  })
  
  
  output$stateSelectorUI <- renderUI({
    req(country_states())
    states <- unique(country_states()$shapeName)
    print("states:")
    print(states)
    print("First State in the list:")
    print(states[1])
    
    fluidRow(
      column(2,
             selectInput("state", "Choose a state:", choices = states, selected = states[1]),  
             checkboxGroupInput("parameters", "Choose parameters to plot:", 
                                choices = names(param_mapping),
                                selected = c())
      ),
      column(10, plotlyOutput("timeSeriesPlot", height = "600px"))
    )
  })
  
output$timeSeriesPlot <- renderPlotly({
  req(input$state)
  
  # Filter data for selected state
  state_data <- filter(time_series_data, shapeName == input$state)
  
  # Check if state_data is empty
  if (nrow(state_data) == 0) {
    return(NULL)
  }
  
  # Create the base plot with tree percentage
  p <- plot_ly(state_data, x = ~year, y = ~tree_percentage, type = 'scatter', mode = 'lines+markers',
               name = 'Tree Percentage', 
               text = rep(paste('Year:', state_data$year, '<br>Tree Percentage:', round(state_data$tree_percentage, 2)), each = 1),
               hoverinfo = 'text', 
               line = list(color = 'blue'))
  
  # Add selected parameters
  if (length(input$parameters) > 0) {
    for (param in input$parameters) {
      original_param <- param_mapping[param]
      p <- p %>% add_trace(y = state_data[[original_param]], mode = 'lines+markers',
                           name = param, 
                           text = rep(paste('Year:', state_data$year, '<br>', param, ':', round(state_data[[original_param]], 2)), each = 1),
                           hoverinfo = 'text')
    }
  }
  
  # Add annotations if necessary
  if (nrow(state_data) > 0) {
    p <- add_annotations(p,
                         x = max(state_data$year),
                         y = max(state_data$tree_percentage),
                         text = paste("Total Area of ", input$state, ": ", round(state_data$total_area_km2[1], 2), "km²")
    )
  }
  
  p <- layout(p, title = paste("Time Series for", input$state),
              xaxis = list(title = 'Year'),
              yaxis = list(title = 'Value'))
  p
})

  
  
  
  
  
  
}

shinyApp(ui, server)