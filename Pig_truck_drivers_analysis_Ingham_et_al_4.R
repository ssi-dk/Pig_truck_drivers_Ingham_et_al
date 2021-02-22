################################################ #
#### Description                              ####
################################################ #
# Prepares the tree for ploting, adds metadata and plots the tree. 

# Note that legend and text elements in the published figures were adjusted 
# using Adobe Illustrator 

# written by Raphael N. Sieber, 2020, 
# Statens Serum Institute, Denmark

################################################ #
## Requirements                               ####
################################################ #
# - Newick-tree file
# - metadata file

################################################ #
#### Load packages, files & functions         ####
################################################ #
library(ggtree)
library(ggplot2)
library(phangorn)
library(openxlsx)
library(pals)

my_color_ramp_creator <- function(data, sort=TRUE) {
  library(pals)
  if(sort) d_unique <- unique(sort(data)) else d_unique <- unique(data)
  l <- length(d_unique)
  if(l<3) cols = c("#C0C0C0", "#ca0020") else
    if(l<9) cols = brewer.set1(l) else
      if(l<26) cols = cols25(l) else
        if(l<37) cols = polychrome(l) else
          cols = rep(polychrome(), ceiling(l/36))[1:l]
  names(cols) <- d_unique
  cols[is.na(names(cols))] <- "#C0C0C0"
  return(cols)
}

my_break_branch <- function(tree_view, node=NULL, labels=NULL, dist=T, factor=10, gap=NULL, seg_width = 0.7, text_size=2, line_size=1) {
  require(ggtree)
  df <- tree_view$data
  if(is.null(node) & is.null(labels)) stop("Error: no node or label given") else
    if(is.null(node) & !is.null(labels)) {
      if(length(labels)>1) node = MRCA(tree_view, labels) else 
        node = which(df$label == labels)
    } 
  br.row.no <- which(df$node == node)
  br.row <- df[br.row.no,]
  x.fact <- (br.row$x - df[br.row$parent,]$x) / br.row$branch.length
  br.row.xred <- x.fact * br.row$branch.length * (1 - 1/factor)
  br.row$branch.length <- br.row$branch.length/factor
  br.row.x0 <- br.row$x
  br.row$x <- br.row$x - br.row.xred
  tree_view$data[br.row.no,] <- br.row

  if(!br.row$isTip){
    sp.df <- tidytree:::offspring(df, node)
    sp <- sp.df$node
    sp.df$x <- sp.df$x - br.row.xred
    tree_view$data[sp,] <- sp.df
  }
  
  if (is.null(gap)) seg_d <- max(df$x)/200 else
    seg_d <- gap
  seg_x <- df[br.row$parent,]$x + (br.row$x - df[br.row$parent,]$x)/2
  seg_width <- seg_width
  seg_skew <- tan(35)*seg_width/max(df$y)*max(df$x)
  if (class(tree_view$coordinates)[1]=="CoordPolar") 
    seg_angle <- ifelse(abs(br.row$angle) < 90 | abs(br.row$angle) >= 270, br.row$angle, br.row$angle-180) else
    seg_angle <- 0
  tree_view <- tree_view +
    do.call(geom_polygon, list(x = c(seg_x-seg_skew-seg_d, seg_x-seg_skew+seg_d, seg_x+seg_skew+seg_d, seg_x+seg_skew-seg_d, rep(NA, nrow(df)-4)),
                               y = c(br.row$y-seg_width, br.row$y-seg_width, br.row$y+seg_width, br.row$y+seg_width, rep(NA, nrow(df)-4)),
                               fill = c(rep("white", 4), rep(NA, nrow(df)-4)), color="white", lwd=NA)) +
    do.call(geom_segment, c(list(x = c(seg_x-seg_skew-seg_d, rep(NA, nrow(df)-1)), xend = seg_x+seg_skew-seg_d, color="#000000"),
                            (l <- list(y = br.row$y-seg_width, yend = br.row$y+seg_width, alpha=1, size=line_size, na.rm = T)))) +
    do.call(geom_segment, c(list(x = c(seg_x-seg_skew+seg_d, rep(NA, nrow(df)-1)), xend = seg_x+seg_skew+seg_d, color="#000000"), l)) 
  if (dist) {
    gap_dist <- br.row.x0 - br.row$x
    tree_view <- tree_view +
      geom_text2(x = c(seg_x, rep(NA, nrow(df)-1)), y = br.row$y, label=paste0(signif(gap_dist, 2)), size = text_size, parse=F, na.rm = T, angle=seg_angle)
  } else {
    tree_view <- tree_view +
      geom_text2(x = c(seg_x, rep(NA, nrow(df)-1)), y = br.row$y, label=paste0("x", factor), size = text_size, parse=F, na.rm = T, angle=seg_angle)
    
  }
  return(tree_view)
}



################################################ #
#### Load, root & adjust tree                 ####
################################################ #
tree <- read.tree(file = "./trees/PigTrDr_27102020.nwk")
tree <- ladderize(reorder(tree))
outgr <- c("H_1953", "H_2046", "H_LY19990171", "H_ST20071083", 
           "H_ST20082015", "H_ST20090121", "H_ST20091155", 
           "H_ST20091526", "H_ST20091826", "H_ST20100011")
all(outgr%in%tree$tip.label)
rooted_tree <- root(unroot(tree), outgroup = outgr, resolve.root = T, edgelabel = T)

## make the root having equal edge lengths on both sides.
root_edges <- which(rooted_tree$edge[,1]==getRoot(rooted_tree))
rooted_tree$edge.length[root_edges] <- mean(rooted_tree$edge.length[root_edges]); rm(root_edges)

tree <- rooted_tree

################################################ #
#### Import metadata                          ####
################################################ #
meta <- read.xlsx("./R_DATA/PigTrDr_isolate_metadata.xlsx", 1)
meta <- meta[!is.na(meta$Lineage),]
meta$MRSA_carrier_status <- gsub(" ", "", meta$MRSA_carrier_status)

meta_extra <- cbind(tree$tip.label[!tree$tip.label%in%meta$Isolate_ID], NA,NA,NA,NA,NA,NA,NA,NA,NA)
colnames(meta_extra) <- colnames(meta)

meta <- rbind(meta, meta_extra)

meta$study <- ifelse(grepl("D0|S0|T0", meta$Isolate_ID), "this_study", ifelse(grepl("._", meta$Isolate_ID), "mbio12", "mbio18"))

identical(sort(tree$tip.label), sort(meta$Isolate_ID)) # TRUE if tip labels match metadata samples

################################################ #
#### Make dataframe for Lineages and extract information for all samples ####
################################################ #
Lin_df <- data.frame(node = c(MRCA(tree, c("T08_2_2_isolate1", "P_34-M-B-1_11", "T02_1_1_isolate3")),
                              MRCA(tree, c("D09_8_0_isolate1", "D03_1_1_isolate1", "T09_2_2_isolate1", "T07_3_1_isolate4")),
                              MRCA(tree, c("P_7413532-2", "T04_1_1_isolate1", "T04_3_1_isolate1", "D06_1_2_isolate1", "D08_3_2_isolate1", "T09_3_1_isolate1", "D02_4_1_isolate1")) ), 
                     Lin_col = c(L1="#386cb0", L2="#7fc97f", L3="#ef3b2c"),
                     col = c("#C3D2E7", "#D8EED8", "#FAC4BF"), 
                     extend_base = 0.000117,
                     l_cex = 4)

Lin_df$exto <- Lin_df$extend_base + c(0, 0, 0)
Lin_df$lab_x <- Lin_df$exto - c(0.00007, 0.00001, 0.00001)

################################################ #
#### Prepare plotting data                    ####
################################################ #

discrete_columns <- c("Time_point_within_day_long", "Sample_type", "MRSA_carrier_status", "Household")
hm_df <- meta[,discrete_columns]
rownames(hm_df) <- meta$Isolate_ID

cols <- c(
  A_morning = "#b2e7ea",
  B_after_first_unloading = "#008c95",
  C_after_second_unloading = "#002325",
  B_Intermittent_MRSA_carrier = "#307fed",
  C_MRSA_carrier = "#b14441",
  D_Truck_of_intermittent_MRSA_carrier = "#abccf7",
  E_Truck_of_MRSA_carrier = "#d8a1a0",
  driver = "#00AFBB",
  spouse = "#FC4E07",
  truck = "#e0c600"
)

hh_cols <- my_color_ramp_creator(sort = F, rev(sort(unique(meta$Household[!is.na(meta$Household)]))))
study_cols <- c(mbio12="grey50", mbio18="grey20", this_study="indianred")
cols <- c(cols, hh_cols)

################################################ #
#### Make plot with ggtree                    ####
################################################ #
p <- ggtree(tree, lwd=0.3, layout = "fan", right = T, open.angle = 0) %<+% meta # Basic tree and link metadata
p <- p + geom_tippoint(aes(color=study), cex=0.4, show.legend=F) # Points at tips

# Add lineage highlights
for (i in rownames(Lin_df)){ 
  p <- p + do.call(geom_hilight, c(list(Lin_df[i, "node"], fill=Lin_df[i, "col"]) , extendto = Lin_df[i, "exto"])) + 
    annotate(geom = "label", x=Lin_df[i ,"lab_x"], y=mean(as.numeric(ggtree:::get_cladelabel_position_(p[["data"]], Lin_df[i, "node"])[1,2:3])),
             label=i, color="black", fill=Lin_df[i, "col"], cex=Lin_df[i, "l_cex"], label.padding = unit(0.25, "lines"),
             label.r = unit(0, "lines"), label.size = 0.5, na.rm = FALSE)
} 
# move Highlight-layers to bottom
hl <- sapply(p$layers, function(x) class(x$stat)[1])=="StatHilight"
p$layers <- c(p$layers[hl], p$layers[!hl]); rm(hl) 

# Break branches
p <- my_break_branch(tree_view=p, labels="P_089B_2", factor=3, text_size = 1, line_size = 0.2, gap = 0.000005) %>%# collaps branch
  my_break_branch(labels=c("H_P23-13_HF-724402", "H_P23-12_HF-80520", "P_P23-14_SD4.1"), factor=3, text_size = 1, line_size = 0.2, gap = 0.000005) %>% # collaps branch
  my_break_branch(labels=c("P_M2009_10003479"), factor=1.5, text_size = 1, line_size = 0.2, gap = 0.000005) # collaps branch

# p
# Add scales, theme parameters, etc.
q <- gheatmap(p, hm_df, width = 0.2, colnames_angle = 90, font.size = 2, colnames_position = "bottom", hjust = 1, color=NA,offset = -0.00001) +  
  scale_color_manual("Color legend", values = study_cols) +
  scale_fill_manual("Color legend", values = cols) +
  theme(panel.background = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = NA,colour = NA),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key.size = unit(9,units = "pt")
  ) +
  xlim(-0.00001,NA) +
  geom_treescale(y=Ntip(tree)/2, x = 0.00001, offset=6, width = 0.00005, fontsize = 2) + geom_rootedge(rootedge = 0.00001)
q
# This throws some errors which can be ignored. 

ggsave(q, filename = "Figures/Phylogeny_Main.pdf", device=cairo_pdf, width = 8, height = 9, bg=NA)

sessionInfo()