library(tidyverse)
library(reshape2)
library(plotly)

options(warn=-1)

#### get_ploty_depth_plot
# Breif: Plots the depth and found targets over region sequenced
# Input:
# detected_var_file, (str) path to file containing all targets found
# target_table_file, (str) path to file containing all targets within region
# bait_table_file, (str) path to file detailing the targeted regions
# gdf_file, (str) path to genomic deth file of target
# exon_table_file, (str) path to exon table for target
# gene, (str) Name of the gene
# Output:
# plotly interactive plot With exons, introns, detected variants, target variants and read depth
#

get_plotly_depth_plot <- function(detected_var_file, target_table_file, bait_table_file, gdf_file, exon_table_file, gene){
detected_variants <- read.csv(
  detected_var_file, sep = "\t", stringsAsFactors = F, row.names = NULL)

exon_table             <- read.csv(exon_table_file, sep = ",", as.is=T, check.names = F) %>% filter(Start != "")
exon_table$Start       <- as.numeric(sapply(exon_table$Start, function(x) str_remove_all(x, ",")))
exon_table$End         <- as.numeric(sapply(exon_table$End, function(x) str_remove_all(x, ",")))
exon_table$Exon_Intron <- exon_table[, "Exon / Intron"]

bait_table <- read.table(bait_table_file, sep = "\t", as.is=T, check.names = F,
                         col.names = c("Chr", "Start", "End", "Name")) %>%
  separate(Name, c("Target", "Gene"), "_") %>%
  filter(Gene == gene)

row.names(bait_table) <- bait_table$Target
bait_table            <- bait_table[order(bait_table$Start), ]
bait_table$Chr        <- as.character(bait_table$Chr)


target_table <- read.table(target_table_file, sep = "\t",
                           col.names = c("Chr", "Start", "End", "Rsid", "Gene"))

gdf <- read.csv(gdf_file, sep = "\t") %>%
  separate(Locus, c("Chr", "Pos"), ":")
gdf$Pos <- as.numeric(gdf$Pos)
gdf$Chr <- sapply(gdf$Chr, function(x) gsub("chr", "", x))


bait_table  <- bait_table[order(bait_table$Start), ]
pos_bait_table <- data.frame(Chr=numeric(0), Pos=numeric(0), Target=character(0), Gene=character(0))
for (i in 1:nrow(bait_table)) {
  start_pos <- bait_table[i, "Start"] - 1
  end_pos   <- bait_table[i, "End"] + 1
  padding   <- 200

  pos_range   <- bait_table[i, "Start"]:bait_table[i, "End"]
  len_pos     <- length(pos_range)

  target_pos  <- data.frame(Chr=rep(bait_table[i, "Chr"], len_pos),
                            Pos=pos_range,
                            Target=rep(bait_table[i, "Target"], len_pos),
                            Gene=rep(bait_table[i, "Gene"], len_pos))


  if (!min(bait_table$Start) == start_pos + 1){
    padding_min <- data.frame(Chr=rep(bait_table[i, "Chr"], padding),
                              Pos=(start_pos-padding + 1):start_pos,
                              Target=rep(paste0("intron", i, "min"), padding),
                              Gene=rep(bait_table[i, "Gene"], padding))
    tmp <- rbind(padding_min, target_pos)
  } else {
    tmp <- target_pos
  }

  if (!max(bait_table$End) == end_pos - 1){
    padding_max <- data.frame(Chr=rep(bait_table[i, "Chr"], padding),
                            Pos=end_pos:(end_pos + padding - 1),
                            Target=rep(paste0("intron", i, "max"), padding),
                            Gene=rep(bait_table[i, "Gene"], padding))
    tmp <- rbind(tmp, padding_max)
  }



  pos_bait_table <- rbind(pos_bait_table, tmp)
}

# Plotting
plot_table <- merge(pos_bait_table,gdf,by = "Pos", all.x = T)
min_cov    <- min(plot_table$Average_Depth_sample, na.rm = T)
plot_table$plot_range <- 1:nrow(plot_table)
plot_table_exon   <- plot_table %>% filter(grepl("exon", Target)) %>%  group_by(Target)
plot_table_intron <- plot_table %>% filter(!grepl("exon", Target)) %>%  group_by(Target)
plot_table_rsids  <- target_table %>%
  filter(Gene == gene) %>%
  melt(id=c("Chr", "Rsid", "Gene"), value.name = "Pos") %>%
  {merge(., plot_table, all.x=T)} %>% group_by(Rsid)
plot_table_detected_variants <- detected_variants %>%
  filter(GENE==gene) %>%  {merge(., plot_table, by.x="POS", by.y="Pos", all.x=T)}

tic_labels     <- sort(unlist(bait_table[, c("Start", "End")]))
pos_range_tics <- plot_table %>% filter(Pos %in% tic_labels)

color_pallete <- c("0275d8", "f0ad4e", "292b2c", 'rgba(92,128,42,1)', 'rgba(213,30,37,1)')
background_color <- "F7F7F7"
titlefont <- list(
  size = 22
)

fig <- plot_ly(plot_table_exon, x=~plot_range, y=min_cov, text =~Target,
            type="scatter", mode="lines", size=0.5, name="Exoner", line=list(color=color_pallete[2])) %>%
  add_trace(data=plot_table_intron, x=~plot_range, y=min_cov, text ="Intron",
            type="scatter", mode="lines",size=0.4 ,name="Introner", line=list(color=color_pallete[3]), inherit=F) %>%
  add_trace(data=plot_table, x=~plot_range, y=~Average_Depth_sample,
            type="scatter", mode="lines", name="Läsdjup", line=list(color=color_pallete[1]), inherit=F) %>%
  add_trace(data=plot_table_rsids, x=~plot_range, y=min_cov-20, text=~Rsid,
            name="Variant targets", size=0.6, type="scatter", mode="markers", inherit = F,
            marker = list(color=color_pallete[4], line=list(color=color_pallete[4]))) %>%
  layout(
    title=list(
      text=paste("Läsdjup över", gene),
      font=list(size=30)),
    yaxis= list(
      title="Läsdjup",
      titlefont=titlefont
      ),
    xaxis = list(
      title="Position",
      tickmode  = "array",
      ticktext = pos_range_tics$Pos,
      tickvals  = pos_range_tics$plot_range,
      titlefont=titlefont
    ),
    autosize = T,
    margin = list(l=50, r=50, b=50, t=100, pad=3),
    legend = list(font=list(size=18)),
    plot_bgcolor=background_color
  )

if (nrow(plot_table_detected_variants) !=0){
  fig <- fig %>%  add_trace(data=plot_table_detected_variants, x=~plot_range, y=min_cov+20, text=~ID,
                            name="Hittad variant", size=0.9, type="scatter", mode="markers", inherit = F,
                            marker=list(color=color_pallete[5], line=list(color=color_pallete[5])))
}

return(fig)
}
