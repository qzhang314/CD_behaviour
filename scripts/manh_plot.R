args <- commandArgs(T)
cov_type <- as.character(args[1])

library(ggplot2)
library(data.table)
library(tidyverse)

gwasResults <- fread(paste0("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_",cov_type,"_maf01.regenie"),fill=TRUE)

gwasResults <- as.data.frame(gwasResults)

gwasResults$logp <- gwasResults$LOG10P
gwasResults$CHR <- as.numeric(gwasResults$CHROM)
gwasResults$BP <- as.numeric(gwasResults$GENPOS)


don <- gwasResults %>%

  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%

  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%

  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%

  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


axisdf$CHR_lab <- ifelse(axisdf$CHR == 23, "X", as.character(axisdf$CHR))

axisdf$CHR_lab2 <- ifelse(
  axisdf$CHR %% 2 == 1,
  paste0(axisdf$CHR_lab, "\n"),   # odd chromosomes on first line
  paste0("\n", axisdf$CHR_lab)    # even chromosomes on second line
)

p <- ggplot(don, aes(x = BPcum, y = LOG10P)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.6) +
  scale_color_manual(values = rep(c("grey", "deepskyblue1"), 22)) +
  scale_x_continuous(
    labels = axisdf$CHR_lab2,
    breaks = axisdf$center,
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  geom_hline(yintercept = -log10(5e-8), col = "violetred") +
  xlab("Chromosome") +
  ylab("-log10(P)") +
  ylim(c(0, 25))



#p <- ggplot(don, aes(x=BPcum, y=LOG10P)) +
#
#    # Show all points
#    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
#    scale_color_manual(values = rep(c("grey", "deepskyblue1"), 22 )) +
#
#    # custom X axis:
#    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
#    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
#
#    # Custom the theme:
#    theme_bw() +
#    theme(
#      legend.position="none",
#      axis.line = element_line(),
#      axis.line.y.right = element_blank(),
#      axis.line.x.top = element_blank(),
#      panel.grid.major.x = element_blank(),
#      panel.grid.minor.x = element_blank(),
#      axis.text=element_text(size=12,face="bold"),
#      axis.title=element_text(size=12,face="bold")
#    )+geom_hline(yintercept=-log10(5e-8),col="violetred")+
#    xlab("Chromosome")+
#    ylab("-log10(P)")

ggsave(file=paste0("~/CD_behaviour/figures/manh_plot_",cov_type,".png"),dpi=300,width=15,height=6)
