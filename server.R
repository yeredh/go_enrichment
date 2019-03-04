library(shiny)
library( scales )
library(ggrepel)
library(colorspace)
library(xtable)
library(ggplot2)
library(gridExtra)
library(rafalib)
# ==== Functions ====
REVIGOplot = function(one.data){
  p1 <- ggplot( data = one.data );
  p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
  # p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
  
  p1 = p1 +  scale_color_gradientn(colours = rev(heat_hcl(7)), limits = c( min(one.data$log10_p_value), max(one.data$log10_p_value)) ) + labs(colour="-log10(p)",size="Term Frequency")
  
  p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.45) )) + scale_size_area();
  p1 <- p1 + scale_size( range=c(2, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  ex <- one.data #[ one.data$dispensability < quantile(one.data$dispensability,0.9), ];
  # p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 2.5 );
  p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 2.5 );
  p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
  p1 <- p1 + theme(legend.key = element_blank()) ;
  one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
  one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
  p1 <- p1 + xlim(min(one.data$plot_X)-2*one.x_range/10,max(one.data$plot_X)+2*one.x_range/10);
  p1 <- p1 + ylim(min(one.data$plot_Y)-2*one.y_range/10,max(one.data$plot_Y)+2*one.y_range/10);
  
  (p1 = p1 + labs(colour="-log10(p)",size="Term Frequency"))
  return(p1)
}

REVIGOplot2 = function(one.data){
  p1 <- ggplot( data = one.data );
  p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = OddsRatio, size = plot_size), alpha = I(0.6) ) + scale_size_area();
  # p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$OddsRatio), 0) );
  
  p1 = p1 +  scale_color_gradientn(colours = rev(heat_hcl(7)), limits = c( min(one.data$OddsRatio), max(one.data$OddsRatio)) ) + labs(colour="-log10(p)",size="Term Frequency")
  
  p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.45) )) + scale_size_area();
  p1 <- p1 + scale_size( range=c(2, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  ex <- one.data #[ one.data$dispensability < quantile(one.data$dispensability,0.9), ]; 
  # p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 2.5 );
  p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 2.5 );
  p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
  p1 <- p1 + theme(legend.key = element_blank()) ;
  one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
  one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
  p1 <- p1 + xlim(min(one.data$plot_X)-2*one.x_range/10,max(one.data$plot_X)+2*one.x_range/10);
  p1 <- p1 + ylim(min(one.data$plot_Y)-2*one.y_range/10,max(one.data$plot_Y)+2*one.y_range/10);
  
  (p1 = p1 + labs(colour="Odds Ratio",size="Term Frequency"))
  return(p1)
}

DisplayPlot = function(one_data){
  if(nrow(one_data) > 0){
    p = REVIGOplot2(one_data)
  }else{
    # empty ggplot
    p = ggplot(one_data, aes(plot_X, plot_Y)) + labs (y = "semantic space x", x = "semantic space y") + theme_bw()
  }
  return(p) 
}

# ==== Data ====
# cut-offs for GO enrichment
p_cut = c(0.0001,0.0005,0.001,0.005,0.01,0.05)
top_n = c(100,250,500,750,1000,1500,2000)

# load results
revigo_probe_pcut = lapply(p_cut,function(x){
  readRDS(paste0("data/revigo_probe_pcut_",format(x),".RDS"))
})
names(revigo_probe_pcut) = as.character(p_cut)

table_probe_pcut = lapply(p_cut,function(x){
  readRDS(paste0("data/toptable_probe_pcut_",format(x),".RDS"))
})
names(table_probe_pcut) = as.character(p_cut)

revigo_fpkm_pcut = lapply(p_cut,function(x){
  readRDS(paste0("data/revigo_fpkm_pcut_",format(x),".RDS"))
})
names(revigo_fpkm_pcut) = as.character(p_cut)

table_fpkm_pcut = lapply(p_cut,function(x){
  readRDS(paste0("data/toptable_fpkm_pcut_",format(x),".RDS"))
})
names(table_fpkm_pcut) = as.character(p_cut)

revigo_probe_topn = lapply(top_n,function(x){
  readRDS(paste0("data/revigo_probe_top_",format(x),".RDS"))
})
names(revigo_probe_topn) = as.character(top_n)

table_probe_topn = lapply(top_n,function(x){
  readRDS(paste0("data/toptable_probe_top_",format(x),".RDS"))
})
names(table_probe_topn) = as.character(top_n)

revigo_fpkm_topn = lapply(top_n,function(x){
  readRDS(paste0("data/revigo_fpkm_top_",format(x),".RDS"))
})
names(revigo_fpkm_topn) = as.character(top_n)

table_fpkm_topn = lapply(top_n,function(x){
  readRDS(paste0("data/toptable_fpkm_top_",format(x),".RDS"))
})
names(table_fpkm_topn) = as.character(top_n)



shinyServer(
  function(input, output) {
    output$cut_off = renderUI({
      selectInput(
        inputId = "pick",
        label = "Selection",
        choices = list(
          "q-value < k" = "p_cut",
          "Top n genes"="top_n"
        ),
        selected = "p_cut"
      )
    })
    output$select_cut = renderUI({
      if(is.null(input$pick)){
        return()
      }
      if(input$pick == "p_cut"){
        selectInput(
          inputId = "p_ind",
          label = "k",
          choices = as.list(p_cut),
          selected = "5e-04"
        )
      }else{
        selectInput(
          inputId = "n_ind",
          label = "n",
          choices = as.list(top_n),
          selected = "2000"
        )
      }
    })
    
    output$plot = renderPlot({
      # q-value cut-off
      if(input$pick == "p_cut"){
        p_lst = list(
          p_fpkm  = DisplayPlot(revigo_fpkm_pcut[[as.character(input$p_ind)]]) + ggtitle("FPKM"),
          p_probe = DisplayPlot(revigo_probe_pcut[[as.character(input$p_ind)]]) + ggtitle("Probes")
        )
        grid.arrange(grobs=p_lst, ncol=length(p_lst))
      }
      # top genes
      if(input$pick == "top_n"){
        p_lst = list(
          p_fpkm  = DisplayPlot(revigo_fpkm_topn[[as.character(input$n_ind)]]) + ggtitle("FPKM"),
          p_probe = DisplayPlot(revigo_probe_topn[[as.character(input$n_ind)]]) + ggtitle("Probes")
        )
        grid.arrange(grobs=p_lst, ncol=length(p_lst))
      }
    })
    
    output$probe_table <- renderTable({
      # q-value cut-off
      if(input$pick == "p_cut"){
        my_tab = table_probe_pcut[[as.character(input$p_ind)]][,c("GO.ID","Term","Annotated","Significant","qval","OddsRatio")]
        my_tab$OddsRatio = as.character(round(my_tab$OddsRatio,3))
        my_tab$qval = as.character(format(my_tab$qval,scientific = T, digits = 3))
        return(my_tab)
      }
      # top genes
      if(input$pick == "top_n"){
        my_tab = table_probe_topn[[as.character(input$n_ind)]][,c("GO.ID","Term","Annotated","Significant","qval","OddsRatio")]
        my_tab$OddsRatio = as.character(round(my_tab$OddsRatio,3))
        my_tab$qval = as.character(format(my_tab$qval,scientific = T, digits = 3))
        return(my_tab)
      }
    },digits=0,rownames=T)
    
    output$fpkm_table <- renderTable({
      # q-value cut-off
      if(input$pick == "p_cut"){
        my_tab = table_fpkm_pcut[[as.character(input$p_ind)]][,c("GO.ID","Term","Annotated","Significant","qval","OddsRatio")]
        my_tab$OddsRatio = as.character(round(my_tab$OddsRatio,3))
        my_tab$qval = as.character(format(my_tab$qval,scientific = T, digits = 3))
        
        return(my_tab)
      }
      # top genes
      if(input$pick == "top_n"){
        my_tab = table_fpkm_topn[[as.character(input$n_ind)]][,c("GO.ID","Term","Annotated","Significant","qval","OddsRatio")]
        my_tab$OddsRatio = as.character(round(my_tab$OddsRatio,3))
        my_tab$qval = as.character(format(my_tab$qval,scientific = T, digits = 3))
        return(my_tab)
      }
    },digits=0,rownames=T)
    
    output$or_plot = renderPlot({
      # q-value cut-off
      if(input$pick == "p_cut"){
        fpkm_data = table_fpkm_pcut[[as.character(input$p_ind)]]
        probe_data = table_probe_pcut[[as.character(input$p_ind)]]
      }
      # top genes
      if(input$pick == "top_n"){
        fpkm_data = table_fpkm_topn[[as.character(input$n_ind)]]
        probe_data = table_probe_topn[[as.character(input$n_ind)]]
      }
      
      if(nrow(fpkm_data) == 0 & nrow(probe_data) == 0){
        plot.new()
      }
      if(nrow(probe_data) == 0 & nrow(fpkm_data) > 0){
        mypar(1,1)
        boxplot(list("FPKM" = fpkm_data$OddsRatio),main="FPKM",ylab="Odds Ratio")
      }
      if(nrow(fpkm_data) == 0 & nrow(probe_data) > 0){
        mypar(1,1)
        boxplot(list("Probes" = probe_data$OddsRatio),main="Probes",ylab="Odds Ratio")
      }
      if(nrow(fpkm_data) > 0 & nrow(probe_data) > 0){
        mypar(1,2)
        boxplot(list("FPKM" = fpkm_data$OddsRatio,"Probes" = probe_data$OddsRatio))
        plot(
          y = fpkm_data$OddsRatio[match(probe_data$GO.ID,fpkm_data$GO.ID)],
          x = probe_data$OddsRatio,
          ylab="FPKM", xlab = "Probes"
        )
        abline(c(0,1))
      }
    })
    
  }
)