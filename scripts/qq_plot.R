args <- commandArgs(T)
cov_type <- as.character(args[1])

library(ggplot2)
library(data.table)

gwas <- fread(paste0("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_",cov_type,"_maf01.regenie"),fill=TRUE)

gwas <- as.data.frame(gwas)

pvals <- 10^-gwas$LOG10P

# Clean the p-values
pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals < 1]

# Chi-square transformation and lambda calculation
chisq_obs <- qchisq(pvals, df = 1, lower.tail = FALSE)
lambda_gc <- median(chisq_obs) / qchisq(0.5, df = 1)

print(lambda_gc)

# Prepare expected and observed -log10(p)
n <- length(pvals)
expected <- -log10(ppoints(n))
observed <- -log10(sort(pvals))

# Data for plot
qq_df <- data.frame(expected = expected, observed = observed)

# Generate QQ plot
ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_bw()+
  theme(axis.text=element_text(size=12,face="bold"),
  	axis.title=element_text(size=12,face="bold"))+
  xlab("Expected -log10(P)")+
  ylab("Observed -log10(P)")

ggsave(file=paste0("~/CD_behaviour/figures/qq_plot_",cov_type,".png"),dpi=300,width=8,height=8)
