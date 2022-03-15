library(httr)
library(readxl)
library(ggplot2)
library(dplyr)
library(ieugwasr)
library(scales)
library(ggsignif)
library(patchwork)

rm(list = ls(all.names = TRUE))

res <- read.table("../results/pqtls_eqtlgen.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
go_res <- read.table("../results/go_reactome_res.txt", sep = "\t", header = T, stringsAsFactors = F, fill = T, quote = "")
dgidb_res <- read.table("../results/dgidb_res.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
gtex_res <- read.table("../results/pqtls_gtex_concord.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
hpa_res <- read.table("../results/hpa_genes_expressed_tissues.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")

concord_colour <- "#FF7400"
concord_darker <- "#BF7130"
discord_colour <- "#009999"
discord_darker <- "#1D7373"

################################################################################
# Step 1
# Link to DGId terms
dgidb_pval <- 2.95e-5
res_dgidb <- data.frame("term" = character(),
                         "concordant" = numeric(),
                         "discordant" = numeric(),
                         "concordant_std" = numeric(),
                         "discordant_std" = numeric(),
                         "concordant_pri" = numeric(),
                         "discordant_pri" = numeric()
)

res_drugs <- data.frame("concordant_drugs" = 0,
                        "concordant_drugsn" = 0,
                        "discordant_drugs" = 0,
                        "discordant_drugsn" = 0,
                        "concordant_pri_drugs" = 0,
                        "concordant_pri_drugsn" = 0,
                        "discordant_pri_drugs" = 0,
                        "discordant_pri_drugsn" = 0
                        )

# Non-primary signals only
# Have to do some data re-jigging
#res <- res %>%
#  dplyr::filter(is.na(type_primary))

dgidb_linkage <- function(row, cutoff = 5, ppval = 0.05, epval = 0.05)
{
  ensg <- row[1]
  dgidb_sub <- dgidb_res[dgidb_res$ensg %in% ensg, ]
  if (!nrow(dgidb_sub)) {
    return(NA)
  }
  
  if (as.double(row[5]) > ppval) {
    return(NA)
  }
  
  if (as.double(row[6]) > epval) {
    return(NA)
  }
  
  terms <- strsplit(dgidb_sub$types, ",", fixed = T)[[1]]
  drugs <- strsplit(dgidb_sub$drugnames, ";", fixed = T)[[1]]
  scores <- strsplit(dgidb_sub$scores, ";", fixed = T)[[1]]
  for (j in 1:length(terms)) {
    if (length(terms) < 1) {
      break
    }
    
    this_term <- terms[j]
    if (is.na(this_term)) {
      next
    }
    
    if (startsWith(this_term, " ")) {
      this_term <- sub(".", "", this_term)
    }
    
    if (this_term %in% res_dgidb$term) {
      if (!is.na(row[7])) {
        res_dgidb[res_dgidb$term == this_term, ][[row[7]]] <<- as.numeric(res_dgidb[res_dgidb$term == this_term, ][[row[7]]]) + 1
      }
      if (!is.na(row[8])) {
        res_dgidb[res_dgidb$term == this_term, ][[paste0(row[8], "_std")]] <<- as.numeric(res_dgidb[res_dgidb$term == this_term, ][[paste0(row[8], "_std")]]) + 1
      }
      if (!is.na(row[9])) {
        res_dgidb[res_dgidb$term == this_term, ][[paste0(row[9], "_pri")]] <<- as.numeric(res_dgidb[res_dgidb$term == this_term, ][[paste0(row[9], "_pri")]]) + 1
      }
    } else {
      res_dgidb[nrow(res_dgidb) + 1, ] <<- c(this_term, 
                                             ifelse("concordant" %in% row[6], 1, 0), 
                                             ifelse("discordant" %in% row[6], 1, 0),
                                             ifelse("concordant" %in% row[7], 1, 0),
                                             ifelse("discordant" %in% row[7], 1, 0),
                                             ifelse("concordant" %in% row[8], 1, 0),
                                             ifelse("discordant" %in% row[8], 1, 0)
      )
    }
  }
  
  if (identical(drugs, character(0))) {
    return(NA)
  }
  
  # Drugs
  if (!is.na(row[7])) {
    res_drugs[[paste0(row[7], "_drugsn")]] <<- as.numeric(res_drugs[[paste0(row[7], "_drugsn")]]) + 1
  }
  if (!is.na(row[9])) {
    res_drugs[[paste0(row[9], "_pri_drugsn")]] <<- as.numeric(res_drugs[[paste0(row[9], "_pri_drugsn")]]) + 1
  }
  for (j in 1:length(drugs)) {
    if (is.na(drugs[j]) || is.na(scores[j]) || identical(drugs[j], character(0))) {
      next
    }
    
    if (as.numeric(scores[j]) < cutoff) {
      next
    }
    
    if (!is.na(row[7])) {
      res_drugs[[paste0(row[7], "_drugs")]] <<- as.numeric(res_drugs[[paste0(row[7], "_drugs")]]) + 1
    }
    if (!is.na(row[9])) {
      res_drugs[[paste0(row[9], "_pri_drugs")]] <<- as.numeric(res_drugs[[paste0(row[9], "_pri_drugs")]]) + 1
    }
  }
}

apply(res, 1, dgidb_linkage, ppval = dgidb_pval, epval = dgidb_pval)
res_dgidb[, 2:7] <- sapply(res_dgidb[, 2:7], as.numeric)

# Next, some plotting helpers
# Plot how P value changes concordance
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

scientific_10 <- function(x) {
  ifelse(x >= 0.05, 
         x,
         parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  )
}

# This plot shows how discordance varies with P value.
# This was not included in the final paper but given here for posterity.
plot_pval <- function(dat, epval = 0.05, type_ = "") {
  type_var <- paste0("type", type_)
  
  counts <- data.frame(epval = double(),
                       concordant = integer(),
                       discordant = integer(),
                       concord_ci_lo = double(),
                       concord_ci_up = double(),
                       discord_ci_lo = double(),
                       discord_ci_up = double())
  
  dat <- dat[!is.na(dat[type_var]), ]
  
  for (i in 1:length(epval)) {
    part <- dat[dat$epval <= epval[i] & dat$ppval <= epval[i], ]
    
    concord <- part %>%
      dplyr::group_by((!!rlang::sym(type_var))) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(tot = sum(n)) %>%
      dplyr::mutate(perc = n / tot) %>%
      dplyr::mutate(ci_lo = perc - 1.96 * sqrt((n / tot) * ((1 - (n / tot)) / tot))) %>%
      dplyr::mutate(ci_up = perc + 1.96 * sqrt((n / tot) * ((1 - (n / tot)) / tot)))
    
    success <- part %>%
      dplyr::mutate(success = ifelse((!!rlang::sym(type_var)) == "concordant", 1, 0)) %>%
      dplyr::count(success) %>%
      '['(2, 2)
      
    ans <- prop.test(x = success, n = nrow(part))
    disc <- prop.test(x = nrow(part) - success, n = nrow(part))

    counts <- rbind(counts,
                    c(epval[i], success, nrow(part) - success, ans$estimate, disc$estimate, ans$conf.int[1], ans$conf.int[2], disc$conf.int[1], disc$conf.int[2]))
  }
  
  names(counts) <- c("epval", "con_count", "dis_count", "concordant", "discordant", "concord_ci_lo", "concord_ci_up", "discord_ci_lo", "discord_ci_up")

  title <- paste0("How percentage of concordant pQTLs and eQTLs\nchange with an increasingly stringent P value",
                  ifelse(type_ == "_primary", " (primary signals only)", " (all signals)"))
  p1 <- ggplot(data = counts, aes(x = 1:length(epval), y = concordant)) +
    geom_point(color = concord_colour, size = 3) +
    #geom_point(data = counts, aes(x = epval, y = discordant), color = discord_colour, size = 2) +
    geom_errorbar(aes(ymin = concord_ci_lo, ymax = concord_ci_up, width = 0.1), 
                  color = concord_darker, 
                  alpha = 0.8,
                  size = 1) +
    #geom_errorbar(data = counts, aes(ymin = discord_ci_lo, ymax = discord_ci_up, width = 0.1), color = discord_darker, alpha = 0.8) +
    #scale_x_continuous(trans = reverselog_trans(10), breaks = epval) +
    scale_x_continuous(breaks = 1:length(epval), expand = c(0, 0), limits = c(0.5, 11.5),
                       labels = scientific_10(epval)) +
    scale_y_continuous(limits = c(0.6, 0.8), labels = scales::percent_format(accuracy = 1), 
      expand = c(0, 0)
      #, sec.axis = sec_axis(~.*1, name = "Percentage of discordance")
      ) +
    labs(
      title = title,
      #x = "P value threshold for QTL pair",
      x = NULL,
      y = "Percentage of concordance"
    ) +
    #annotation_custom(gridExtra::tableGrob(
    #    t(counts[c("con_count")])
    #  ), 
    #  xmin = 0.5, xmax = 13.5,
    #  ymin = 0.60, ymax = 0.64
    #) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          text = element_text(size = 20),
          axis.text.x = element_blank()
    )

  p2 <- counts[c("epval", "con_count", "dis_count")] %>%
    tidyr::pivot_longer(c(epval, con_count, dis_count), names_to = "layer", values_to = "label") %>%
    ggplot(aes(x = rep(1:length(epval), each = 3))) +
    geom_text(aes(y = factor(layer, c("dis_count", "con_count", "epval")), label = label), size = 5) +
    labs(y = NULL, x = NULL) +
    scale_x_continuous(breaks = 1:length(epval), expand = c(0, 0), limits = c(0.5, 11.5),
                       labels = scientific_10(epval)) +
    scale_y_discrete(labels = c("Disc pairs", "Conc pairs", "P value")) +
    theme_minimal() +
    theme(axis.line = element_blank(), 
          axis.ticks = element_blank(), axis.text.x = element_blank(),
          panel.grid = element_blank(), strip.text = element_blank(),
          text = element_text(size = 20))
  
  p <- p1 / p2 + plot_layout(heights = c(8, 2))
  #pp <- list(p1, p2)
  #p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v")

  return(p)
}

plot_pval(res, type_ = "_primary", epval = c(0.05, 5e-3, 5e-4, 5e-5, 5e-6, 5e-7, 5e-8, 5e-9, 5e-10, 5e-11, 5e-12))

# Concordance in primary vs non-primary signals
pri_vs_non <- res %>%
  dplyr::filter(ppval < 2.95e-5 & epval < 2.95e-5) %>%
  dplyr::mutate(signal_group = ifelse(is.na(type_primary), 2, 1)) %>%
  dplyr::group_by(signal_group, type) %>%
  dplyr::summarise(n = n())

pvn <- data.frame(
  "concord" = c(pri_vs_non$n[1], pri_vs_non$n[3]),
  "discord" = c(pri_vs_non$n[2], pri_vs_non$n[4]),
  row.names = c("pri", "non"),
  stringsAsFactors = F
)

p <- fisher.test(pvn)$p.value

###############################################################################
# Step 2
# Plotting GTEx results
gtex_res$ppval <- as.double(gtex_res$ppval)
gtex_res$epval <- as.double(gtex_res$epval)
gtex_res <- gtex_res[!startsWith(gtex_res$tissue, "rs"), ] %>%
  tidyr::separate(tissue, c("tissue", "ensg2", "rsid2"), sep = "\\.")

# STable 1
tab <- data.frame(tissue = character(),
                  all = character(),
                  nominal = character(),
                  genome = character())
for (t in unique(gtex_res$tissue)) {
  tab[nrow(tab) + 1, ] <- c(
    t,
    paste0(nrow(gtex_res[gtex_res$tissue == t,]), " (", nrow(gtex_res[gtex_res$tissue == t & gtex_res$signal == 1,]), " / ", nrow(gtex_res[gtex_res$tissue == t & gtex_res$signal != 1,]), ")"),
    paste0(nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 0.05 & gtex_res$ppval < 0.05,]), " (", nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 0.05 & gtex_res$ppval < 0.05 & gtex_res$signal == 1,]), " / ", nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 0.05 & gtex_res$ppval < 0.05 & gtex_res$signal != 1,]), ")"),
    paste0(nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 5e-8 & gtex_res$ppval < 5e-8,]), " (", nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 5e-8 & gtex_res$ppval < 5e-8 & gtex_res$signal == 1,]), " / ", nrow(gtex_res[gtex_res$tissue == t & gtex_res$epval < 5e-8 & gtex_res$ppval < 5e-8 & gtex_res$signal != 1,]), ")")
  )
}

write.table(tab, "../results/gtex_numbers_tab.txt",
            sep = "\t", quote = F, row.names = F)

# Same heuristics
gtex_pval <- 2.95e-5
gtex_rates <- gtex_res %>%
  dplyr::filter(ppval < gtex_pval & epval < gtex_pval) %>%
  dplyr::group_by(tissue, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::mutate(perc = n / tot * 100)

# Box plot of concordance in tissues
gtex_res %>%
  dplyr::filter(ppval < gtex_pval & epval < gtex_pval) %>%
  dplyr::filter(signal == 1) %>%
  dplyr::group_by(ensg, rsid, type) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::mutate(lab = paste0(ensg, "_", rsid)) %>%
  dplyr::filter(tot > 1) %>%
  ggplot(aes(x = n, fill = type)) +
  geom_histogram()
  
###############################################################################
# GTEx linkage with HPA secretome evidence
library(hash)
secretome <- read.table("../results/secretome_genes_of_interest.txt", sep = "\t", header = T)
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
gtex_hpa_pval <- 2.95e-5
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
            file = "../results/gtex_linked_to_secretome.txt",
            sep = "\t", quote = F, row.names = F)

# Heuristics
gtex_secret_counts <- gtex_linked_to_secretome %>%
  #dplyr::filter(signal == 1) %>%
  dplyr::group_by(sec_tissue, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(tot = sum(n)) %>%
  dplyr::mutate(perc = n / tot * 100)

# Figure 4
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
    title = "Concordance of QTL pairs using pQTLs and GTEx tissue-specific eQTLs\nCombined with secretome evidence from HPA",
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
