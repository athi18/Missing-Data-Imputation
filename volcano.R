if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
  install.packages('dplyr')
  
library(ggplot2)
library(ggrepel)
library(dplyr)

df <- read.csv("dfev_dhf.csv",header=TRUE)
df$Expression <- 'NO'
df$Expression[df$log2FC >= 0.15 & df$p_adjusted<=0.05] <- 'UP'
df$Expression[df$log2FC <= -0.15 & df$p_adjusted<=0.05] <- 'DOWN'
  
p2 <- ggplot(df, aes(log2FC, -log(p_adjusted,10))) +
    geom_point(aes(color = Expression), size = 3/5) +
    xlab("log2FC") +
    ylab("-log10P") +
    geom_vline(xintercept = c(0.15,-0.15),linetype="dotted")+
    ggtitle("Dengue Fever vs Dengue Haemorrhagic Fever")+
    scale_color_manual(values = c("gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5)))

top_genes <- bind_rows(
  df %>% 
    filter(Expression == 'UP') %>% 
    arrange(p_adjusted, desc(abs(log2FC))) %>% 
    head(10),
  df %>% 
    filter(Expression == 'DOWN') %>% 
    arrange(p_adjusted, desc(abs(log2FC))) %>% 
    head(10)
)
p3 <-  p2 +
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FC, -log(p_adjusted,10), label = symbol),
                   size = 2)
p3
ggsave("DFev_vs_DHF.png",p3,width=16,height=11.5,dpi=300,units="cm")

#df_dfev <- read.csv("hc_dfev.csv",header=TRUE)
#df_dfev$Expression <- 'NO'
#df_dfev$Expression[df_dfev$log2FC >= 0.15 & df_dfev$p_adjusted<=0.05] <- 'UP'
#df_dfev$Expression[df_dfev$log2FC <= -0.15 & df_dfev$p_adjusted<=0.05] <- 'DOWN'

#df_dhf <- read.csv("hc_dhf.csv",header=TRUE)
#df_dhf$Expression <- 'NO'
#df_dhf$Expression[df_dhf$log2FC >= 0.15 & df_dhf$p_adjusted<=0.05] <- 'UP'
#df_dhf$Expression[df_dhf$log2FC <= -0.15 & df_dhf$p_adjusted<=0.05] <- 'DOWN'
