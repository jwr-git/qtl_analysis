library(httr)
library(readxl)
library(ggplot2)
library(dplyr)
library(ieugwasr)
library(scales)
library(ggsignif)
library(patchwork)
library(grid)
library(gridExtra)

rm(list = ls(all.names = TRUE))

res <- read.table("../results/pqtls_eqtlgen.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
gtex_res <- read.table("../results/pqtls_gtex_concord.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
hpa_res <- read.table("../results/hpa_genes_expressed_tissues.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")

concord_colour <- "#FF7400"
concord_darker <- "#BF7130"
discord_colour <- "#009999"
discord_darker <- "#1D7373"

gtex_res$ppval <- as.double(gtex_res$ppval)
gtex_res$epval <- as.double(gtex_res$epval)
gtex_res <- gtex_res[!startsWith(gtex_res$tissue, "rs"), ] %>%
  tidyr::separate(tissue, c("tissue", "ensg2", "rsid2"), sep = "\\.")

# STable 1
tab <- data.frame(tissue = character(),
                  all = character(),
                  nominal = character())
                  #genome = character())
for (t in unique(gtex_res$tissue)) {
  tab[nrow(tab) + 1, ] <- c(
    t,
    paste0(nrow(gtex_res[gtex_res$tissue == t,]), " (", nrow(gtex_res[gtex_res$tissue == t & gtex_res$signal == 1,]), " / ", nrow(gtex_res[gtex_res$tissue == t & gtex_res$signal != 1,]), ")"),
    paste0(nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 2.95e-5 & gtex_res$ppval < 2.95e-5,]), " (", nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 2.95e-5 & gtex_res$ppval < 2.95e-5 & gtex_res$signal == 1,]), " / ", nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 2.95e-5 & gtex_res$ppval < 2.95e-5 & gtex_res$signal != 1,]), ")")
    #paste0(nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 5e-8 & gtex_res$ppval < 5e-8,]), " (", nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 5e-8 & gtex_res$ppval < 5e-8 & gtex_res$signal == 1,]), " / ", nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 5e-8 & gtex_res$ppval < 5e-8 & gtex_res$signal != 1,]), ")")
  )
}

write.table(tab, "../results/gtex_numbers_tab.txt",
            sep = "\t", quote = F, row.names = F)

# Figure 3
gtex_pval <- 2.95e-5
p1 <- gtex_res %>%
  dplyr::filter(ppval < gtex_pval & epval < gtex_pval) %>%
  dplyr::filter(signal == 1) %>%
  dplyr::group_by(tissue, type) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  ggplot(aes(x = tissue, y = n, fill = type, label = n)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(
    #"QTL pair category",
    "QTL pair category\n(primary signals only)",
    labels = c("Concordant", "Discordant"),
    values = c("concordant" = concord_colour, "discordant" = discord_colour)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 16)) +
  labs(
    #title = "Concordance of QTL pairs using pQTLs and GTEx tissue-specific eQTLs",
    subtitle = "i) Primary QTL pairs only",
    x = "GTEx Tissue",
    #y = "Percentage of concordance"
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_discrete(limits = rev)

p2 <- gtex_res %>%
  dplyr::filter(ppval < gtex_pval & epval < gtex_pval) %>%
  dplyr::group_by(tissue, type) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  ggplot(aes(x = tissue, y = n, fill = type, label = n)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(
    #"QTL pair category",
    "QTL pair category",
    labels = c("Concordant", "Discordant"),
    values = c("concordant" = concord_colour, "discordant" = discord_colour)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 16),
        plot.margin = unit(c(0, 0, 0, 0.7), "cm")) +
  labs(
    #title = "Concordance of QTL pairs using pQTLs and GTEx tissue-specific eQTLs",
    subtitle = "ii) All QTL pairs",
    #x = "GTEx Tissue",
    #y = "Percentage of concordance"
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_discrete(limits = rev)

p <- cowplot::plot_grid(
  p1 + 
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 12)),
  p2 + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 12)),
  align = "h",
  ncol = 2,
  rel_widths = c(1, 0.85)
)
title.grob <- textGrob("Concordance of QTL pairs using plasma pQTLs and GTEx eQTLs", gp=gpar(fontsize=22))
x.grob <- textGrob("Percentage of QTL pairs enriched", gp=gpar(fontsize=16))
y.grob <- textGrob("GTEx Tissue", gp=gpar(fontsize=16), rot=90)
grid.arrange(arrangeGrob(p, top = title.grob, bottom = x.grob, left = y.grob))

# Heuristics and data exploration
gtex_rates <- gtex_res %>%
  dplyr::filter(ppval < gtex_pval & epval < gtex_pval) %>%
  #dplyr::filter(signal == 1) %>%
  dplyr::group_by(tissue, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::mutate(perc = n / tot * 100) %>%
  dplyr::filter(type == "concordant")

# Box plot of concordance in tissues
histo <- gtex_res %>%
  dplyr::filter(ppval < gtex_pval & epval < gtex_pval) %>%
  dplyr::filter(signal == 1) %>%
  dplyr::group_by(ensg, rsid, type) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::mutate(lab = paste0(ensg, "_", rsid)) %>%
  dplyr::filter(tot > 1) %>%
  dplyr::mutate(perc = n / tot) %>%
  dplyr::arrange(lab) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::mutate(concord = ifelse(type == "concordant", perc, 1 - perc))

ggplot(data = histo, aes(x = perc)) +
  geom_histogram(binwidth = 0.05,
                 colour = "black",
                 fill = "grey") +
  scale_x_continuous(breaks = seq(0, 100, 10) / 100, labels = scales::percent, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Distribution of concordance across tissues for GTEx eQTLs paired with plasma pQTLs",
       x = "Percentage of concordant eQTLs across tissues",
       y = "Count of eQTLs") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 16))

###############################################################################
# GTEx linkage with HPA
library(hash)
secretome <- read.table("../results\\secretome_genes_of_interest.txt", sep = "\t", header = T)
gtex_secretome <- base::merge(gtex_res, secretome, by.x = "ensg", by.y = "Gene", all.x = T)

h <- hash::hash()
# HPA <-> GTEx (be sure to use fuzzy on GTEx)
h[["breast"]] <- "Breast_Mammary_Tissue"
h[["bronchus"]] <- "Lung"
h[["cervix, uterine"]] <- "Uterus"
h[["endometrium 1"]] <- "Uterus"
h[["endometrium 2"]] <- "Uterus"
h[["esophagus"]] <- "Esophagus"
h[["fallopian tube"]] <- "NA"
h[["nasopharynx"]] <- "NA"
h[["placenta"]] <- "NA"
h[["seminal vesicle"]] <- "NA"
h[["skin 1"]] <- "Skin"
h[["tonsil"]] <- "NA"
h[["urinary bladder"]] <- "NA"
h[["vagina"]] <- "Vagina"
h[["adrenal gland"]] <- "Adrenal_Gland"
h[["appendix"]] <- "NA"
h[["colon"]] <- "Colon"
h[["duodenum"]] <- "Colon"
h[["gallbladder"]] <- "NA"
h[["heart muscle"]] <- "Heart"
h[["kidney"]] <- "Kidney"
h[["pancreas"]] <- "Pancreas"
h[["rectum"]] <- "NA"
h[["skeletal muscle"]] <- "Muscle_Skeletal"
h[["small intestine"]] <- "Small_Intestine_Terminal_Ileum"
h[["smooth muscle"]] <- "NA"
h[["stomach 1"]] <- "Stomach"
h[["stomach 2"]] <- "Stomach"
h[["testis"]] <- "Testis"
h[["caudate"]] <- "Brain_Caudate_basal_ganglia"
h[["cerebellum"]] <- "Brain_Cerebellum"
h[["cerebral cortex"]] <- "Cortex"
h[["epididymis"]] <- "NA"
h[["hippocampus"]] <- "Brain_Hippocampus"
h[["parathyroid gland"]] <- "Thyroid"
h[["prostate"]] <- "Prostate"
h[["salivary gland"]] <- "Minor_Salivary_Gland"
h[["skin 2"]] <- "Skin"
h[["thyroid gland"]] <- "Thyroid"
h[["oral mucosa"]] <- "NA"
h[["ovary"]] <- "Ovary"
h[["soft tissue 2"]] <- "NA"
h[["lung"]] <- "Lung"
h[["adipose tissue"]] <- "Adipose"
h[["liver"]] <- "Liver"
h[["bone marrow"]] <- "NA"
h[["spleen"]] <- "Spleen"
h[["soft tissue 1"]] <- "NA"
h[["lymph node"]] <- "NA"
h[["hypothalamus"]] <- "Brain_Hypothalamus"
h[["hair"]] <- "NA"
h[["retina"]] <- "NA"
h[["skin"]] <- "Skin"
h[["thymus"]] <- "NA"
h[["cartilage"]] <- "NA"
h[["eye"]] <- "NA"
h[["pituitary gland"]] <- "Pituitary"
h[["dorsal raphe"]] <- "NA"
h[["choroid plexus"]] <- "NA"
h[["lactating breast"]] <- "Breast"
h[["intestine"]] <- "Intestine"
h[["lymphoid tissue"]] <- "NA"
h[["tongue"]] <- "NA"
h[["brain"]] <- "Brain"
h[["blood"]] <- "Whole_Blood"
h[["ductus deferens"]] <- "NA"

# Linkage to HPA tissue expression and secretome
gtex_hpa_pval <- 5e-2
gtex_linked_to_hpa_expression <- data.frame(ensg = character(),
                                            hgnc = character(),
                                            tissue = character(),
                                            signal = integer(),
                                            type = character(),
                                            tissue_hpa = character(),
                                            cell.type = integer())

for (i in 1:nrow(gtex_res)) {
  dat <- gtex_res[i, ]
  dat <- dat[dat$ppval < gtex_hpa_pval & dat$epval < gtex_hpa_pval, ]
  
  if (nrow(dat) < 1) {
    next
  }
  
  # Expression linkage
  if (dat$ensg %in% hpa_res$Gene) {
    hpa_this <- hpa_res[hpa_res$Gene == dat$ensg, ]
    
    for (k in 1:nrow(hpa_this)) {
      tissue <- hpa_this[k, ]$Tissue
      gtex_tissue_lookup <- h[[tissue]]
      if (gtex_tissue_lookup == "NA" || grepl(gtex_tissue_lookup, dat$tissue, fixed = T) == FALSE) {
        next
      }
      
      # Matched -- expressed tissue evidence and GTEx eQTL are the same!
      gtex_linked_to_hpa_expression[nrow(gtex_linked_to_hpa_expression) + 1, ] <- c(
        dat$ensg,
        dat$hgnc,
        dat$tissue,
        dat$signal,
        dat$type,
        hpa_this[k, ]$Tissue,
        hpa_this[k, ]$Cell.type
      )
    }
  }
}

gtex_express_counts <- gtex_linked_to_hpa_expression %>%
  dplyr::count(tissue, type)

gtex_linked_to_hpa_expression %>%
  dplyr::filter(signal == 1) %>%
  dplyr::group_by(tissue, type) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  ggplot(aes(x = tissue, y = n, fill = type, label = n)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(
    #"QTL pair category",
    "QTL pair category\n(primary signals only)",
    labels = c("Concordant", "Discordant"),
    values = c("concordant" = concord_colour, "discordant" = discord_colour)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(
    title = "Concordance between pQTLs and GTEx tissue-specific eQTLs with evidence of expression from HPA",
    subtitle = paste0("QTL pair selection threshold: P < ", gtex_hpa_pval),
    x = "GTEx Tissue",
    y = "Percentage of concordance"
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0))


# Secretome
hpa_pval <- 2.95e-5
gtex_linked_to_secretome <- data.frame(ensg = character(),
                                       hgnc = character(),
                                       snp = character(),
                                       tissue = character(),
                                       signal = integer(),
                                       type = character(),
                                       Gene.Name = character(),
                                       Uniprot = character(),
                                       sec_tissue = character(),
                                       sec_score = integer())
for (i in 1:nrow(gtex_secretome)) {
  dat <- gtex_secretome[i, ]
  
  dat <- dat[dat$ppval < hpa_pval & dat$epval < hpa_pval, ]
  
  if (nrow(dat) < 1) {
    next
  }
  
  if (is.na(dat$X1)) {
    next
  }
  
  # Secretome linkage
  for (j in 14:18) {
    if (is.na(dat[, j])) {
      break
    }
    
    tissue_score <- strsplit(dat[, j], ":")[[1]]
    gtex_tissue_lookup <- h[[tissue_score[1]]]
    if (gtex_tissue_lookup == "NA" || grepl(gtex_tissue_lookup, dat$tissue, fixed = T) == FALSE) {
      next
    }
    
    # Matched -- secretome and GTEx eQTL are the same!
    gtex_linked_to_secretome[nrow(gtex_linked_to_secretome) + 1, ] <- c(
      dat$ensg,
      dat$hgnc,
      dat$rsid,
      dat$tissue,
      dat$signal,
      dat$type,
      dat$Gene.name,
      dat$Uniprot,
      tissue_score[1],
      tissue_score[2]
    )
  }
}

write.table(gtex_linked_to_secretome,
            file = "../results\\gtex_linked_to_secretome.txt",
            sep = "\t", quote = F, row.names = F)

gtex_secret_counts <- gtex_linked_to_secretome %>%
  #dplyr::filter(signal == 1) %>%
  dplyr::group_by(sec_tissue, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::mutate(perc = n / tot * 100)

gtex_linked_to_secretome$sec_tissue <- forcats::fct_rev(factor(gtex_linked_to_secretome$sec_tissue))
gtex_linked_to_secretome %>%
  #dplyr::filter(signal == 1) %>%
  dplyr::group_by(sec_tissue, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  purrr::map_df(rev) %>%
  ggplot(aes(x = sec_tissue, y = n, fill = type, label = n)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(
    "QTL pair category",
    #"QTL pair category\n(primary signals only)",
    labels = c("Concordant", "Discordant"), 
    values = c("concordant" = concord_colour, "discordant" = discord_colour)) +
  labs(
    title = "Concordance of QTL pairs using pQTLs and GTEx eQTLs\nCombined with secretome evidence from HPA",
    subtitle = paste0("QTL pair selection threshold: P < ", hpa_pval),
    x = "HPA Tissue",
    y = "Percentage of concordance"
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  coord_flip() +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 16))
