library(data.table)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)


usage <- function() {
  cat("
description: plot LD heatmap and markers in LD to a significant GWAS hit

usage: get_markers_high-ld_sig-hit.R [ld_file] [sig_marker] [outfolder] [...]

positional arguments:
  ld_file             plink LD file
  sig_marker          name of significant GWAS hit
  outfolder           folder to save plots

optional argument:
  --help                  show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 3) stop(usage(), "missing positional argument(s)")

# get arguments
ld_file <- args[1]
sig_marker <- args[2]
outfolder <- args[3]
if (!dir.exists(outfolder)) dir.create(outfolder)
# ld_file <- "analysis/sig_marker_yield/test_ld_chr1.ld"
# # sig_marker <- "snp.1.306948638"
# sig_marker <- "snp.1.306703117"
# outfolder <- "analysis/sig_marker_yield"




#### load data ----

# load file
ld <- fread(ld_file, header = TRUE, data.table = FALSE)

# get only the last 2 mb of chr
ld_chr_end <- subset(ld, BP_A > 305000000 & BP_B > 305000000)

# get only markers in ld to significant marker
ld_sig_marker <- subset(ld, SNP_A == sig_marker | SNP_B == sig_marker)

# remove big file
rm(ld)




#### ld triangle ----

# order markers by position
markers_to_order <- c(ld_chr_end$BP_A, ld_chr_end$BP_B)
names(markers_to_order) <- c(ld_chr_end$SNP_A, ld_chr_end$SNP_B)
ordered_markers <- sort(markers_to_order[!duplicated(markers_to_order)])

# create a empty matrix with all marker positions
ld_matrix <- matrix(NA, nrow = length(ordered_markers), ncol = length(ordered_markers))
rownames(ld_matrix) <- names(ordered_markers)
colnames(ld_matrix) <- names(ordered_markers)

# transform to 3-column format
ld_matrix <- reshape2::melt(ld_matrix)[reshape2::melt(upper.tri(ld_matrix))$value, ]
colnames(ld_matrix) <- c("SNP_A", "SNP_B", "R2")

# merge R2 values
ld_matrix <- merge(x = ld_matrix, y = ld_chr_end, by = c("SNP_A", "SNP_B"), all.x = TRUE)
ld_matrix <- ld_matrix[, c("SNP_A", "SNP_B", "R2.y")]
colnames(ld_matrix)[3] <- "R2"

# adjust arrow/label coordinates based on where the sig marker is located
# this will avoid any arrow/label get out of bounds
sig_marker_pos <- ordered_markers[which(names(ordered_markers) == sig_marker)]
sig_marker_xend <- sig_marker
sig_marker_yend <- sig_marker
if (sig_marker_pos < quantile(ordered_markers)[2]) {
  # adjust arrow
  sig_marker_x <- sig_marker
  sig_marker_y <- names(ordered_markers)[which(names(ordered_markers) == sig_marker) + 100]
  # adjust label
  sig_marker_vjust <- -0.3
  sig_marker_hjust <- 0.1
}
if (sig_marker_pos >= quantile(ordered_markers)[2] & sig_marker_pos <= quantile(ordered_markers)[4]) {
  # adjust arrow
  sig_marker_x <- names(ordered_markers)[which(names(ordered_markers) == sig_marker) - 25]
  sig_marker_y <- names(ordered_markers)[which(names(ordered_markers) == sig_marker) + 25]
  # adjust label
  sig_marker_vjust <- -0.5
  sig_marker_hjust <- 0.5
}
if (sig_marker_pos > quantile(ordered_markers)[4]) {
  # adjust arrow
  sig_marker_x <- names(ordered_markers)[which(names(ordered_markers) == sig_marker) - 25]
  sig_marker_y <- sig_marker
  # adjust label
  sig_marker_vjust <- 0.5
  sig_marker_hjust <- 1
}
# get closest markers to quantile boundaries
# i'm doing this because x-axis is not numerical and i need names of markers to guide my plotting
quant_marker_2 <- names(ordered_markers)[which(ordered_markers >= quantile(ordered_markers)[2])[1]]
quant_marker_3 <- names(ordered_markers)[which(ordered_markers >= quantile(ordered_markers)[3])[1]]
quant_marker_4 <- names(ordered_markers)[which(ordered_markers >= quantile(ordered_markers)[4])[1]]

# plot ld triangle
plot_triangle <- ggplot(ld_matrix, aes(x = as.character(SNP_B), y = as.character(SNP_A), fill = R2)) +
  # add heatmap
  geom_raster() +
  scale_fill_distiller(palette = "Spectral") +
  # add labels for significant marker
  annotate("segment", x = sig_marker_x, xend = sig_marker_xend, y = sig_marker_y,
           yend = sig_marker_yend, arrow = arrow(length = unit(0.1, "in"))) +
  annotate("text", x = sig_marker_x, y = sig_marker_y, label = sig_marker,
           vjust = sig_marker_vjust, hjust = sig_marker_hjust, fontface = 2) +
  # add labels for physical position
  annotate("segment", x = quant_marker_2, xend = quant_marker_2, y = quant_marker_2, yend = -Inf) +
  annotate("segment", x = quant_marker_3, xend = quant_marker_3, y = quant_marker_3, yend = -Inf) +
  annotate("segment", x = quant_marker_4, xend = quant_marker_4, y = quant_marker_4, yend = -Inf) +
  annotate("text", x = quant_marker_2, y =  -Inf, label = quant_marker_2, vjust = 1.2, hjust = 0.5) +
  annotate("text", x = quant_marker_3, y =  -Inf, label = quant_marker_3, vjust = 1.2, hjust = 0.5) +
  annotate("text", x = quant_marker_4, y =  -Inf, label = quant_marker_4, vjust = 1.2, hjust = 0.5) +
  # adjust theme
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(t = 0.1, r = 0.1, b = 0.2, l = 0.1, unit = "in"))
# plot legend separately
plot_legend <- ggplot(ld_matrix, aes(x = as.character(SNP_B), y = as.character(SNP_A), fill = R2)) +
  geom_raster() +
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
plot_legend <- get_legend(plot_legend)
# plot title separately
plot_title <- ggplot() +
  labs(title = "Pairwise LD",
       subtitle = paste0("chr1:", format(ordered_markers[1], big.mark= ","),
                         "..",format(rev(ordered_markers)[1], big.mark= ","))) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 15))



# open graphic device
plot_name <- paste0(outfolder, "/ld_triangle_end_chr.pdf")
pdf(plot_name, width = 10, height = 7)
grid.newpage()
grid.draw(arrangeGrob(plot_triangle))
grid.draw(arrangeGrob(plot_title))
grid.draw(arrangeGrob(plot_legend, vp = viewport(x = 0.09, y = 0.8)))
# close graphic device for saving plot
dev.off()



#### markers in high ld to sig marker ----

# determine position of marker in ld to significant marker
ld_sig_marker$not_sig_marker_pos <- as.numeric(apply(ld_sig_marker, MARGIN = 1, function(marker_info, sig_marker) {
  
  if (marker_info["SNP_A"] == sig_marker) {
    return(marker_info["BP_B"])
  } else {
    return(marker_info["BP_A"])
  }
  
}, sig_marker = sig_marker))
ld_sig_marker$dist_to_sig_marker <- ld_sig_marker$BP_B - ld_sig_marker$BP_A

# get position of significant marker (for plotting)
sig_marker_pos <- ld_sig_marker[ld_sig_marker$SNP_A == sig_marker, "BP_A"][1]
if (is.na(sig_marker_pos)) sig_marker_pos <- ld_sig_marker[ld_sig_marker$SNP_B == sig_marker, "BP_B"][1]

# plot r2 values along chr
plot_sig_hit <- ggplot(ld_sig_marker, aes(x = not_sig_marker_pos, y = R2)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_vline(xintercept = sig_marker_pos, color = "firebrick") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(labels = function(x) x / 1e6) +
  labs(x = "Physical position (Mb)", y = "R2")
plot_sig_hit_name <- paste0(outfolder, "/ld_to_sig-hit.pdf")
ggsave(plot_sig_hit, filename = plot_sig_hit_name, device = "pdf", width = 10, height = 7)

# create bins for plotting
ld_sig_marker$bins <- cut_width(ld_sig_marker$not_sig_marker_pos, width = 1e5, boundary = 0)
levels(ld_sig_marker$bins) <- sapply(levels(ld_sig_marker$bins), function(level) {
  new_level <- unlist(strsplit(level, split = ","))[1]
  new_level <- as.numeric(substring(new_level, 2))
  new_level <- new_level / 1e6
  return(new_level)
})
# plot ld but in bins
plot_sig_hit_bins <- ggplot(ld_sig_marker, aes(x = bins, y = R2)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Physical position (Mb)", y = "R2")
plot_sig_hit_bins_name <- paste0(outfolder, "/ld_to_sig-hit_bins.pdf")
ggsave(plot_sig_hit_bins, filename = plot_sig_hit_bins_name, device = "pdf", width = 10, height = 7)


# get table with markers in highest ld to signgicant marker
list_markers_high_ld <- subset(ld_sig_marker, R2 >= 0.9)
list_markers_high_ld <- apply(list_markers_high_ld, MARGIN = 1, function(marker_info, sig_marker) {
  
  if (marker_info["SNP_A"] == sig_marker) {
    return(c(marker_info["SNP_B"], marker_info["CHR_B"], marker_info["BP_B"], marker_info["R2"]))
  } else {
    return(c(marker_info["SNP_A"], marker_info["CHR_A"], marker_info["BP_A"], marker_info["R2"]))
  }
  
}, sig_marker = sig_marker)
# reformat table and columns
list_markers_high_ld <- data.frame(t(list_markers_high_ld), stringsAsFactors = FALSE)
colnames(list_markers_high_ld) <- c("marker", "chr", "pos", "r2")
list_markers_high_ld$chr <- as.numeric(list_markers_high_ld$chr)
list_markers_high_ld$pos <- as.numeric(list_markers_high_ld$pos)
list_markers_high_ld$r2 <- as.numeric(list_markers_high_ld$r2)
# reorder table
list_markers_high_ld <- list_markers_high_ld[order(list_markers_high_ld$chr, list_markers_high_ld$pos), ]

# write table
table_name <- paste0(outfolder, "/list_markers_high_ld.txt")
fwrite(list_markers_high_ld, file = table_name, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
