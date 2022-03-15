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
go_res <- read.table("../results/go_reactome_res.txt", sep = "\t", header = T, stringsAsFactors = F, fill = T, quote = "")
dgidb_res <- read.table("../results/dgidb_res.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")

concord_colour <- "#FF7400"
concord_darker <- "#BF7130"
discord_colour <- "#009999"
discord_darker <- "#1D7373"

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

dgidb_linkage <- function(row, cutoff = 5, repl = "discordant", ppval = 0.05, epval = 0.05)
{
  ensg <- row[1]
  dgidb_sub <- dgidb_res[dgidb_res$ensg %in% ensg, ]
  if (!nrow(dgidb_sub)) {
    return(NA)
  }
  
  if (repl != "drop") {
    row <- gsub("partially discordant", repl, row)
  } else {
    row <- gsub("partially discordant", NA, row)
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

# Figure 2
plot_results <- function(dat, res, term, type = "", pval = 0.05) {
  con_var <- paste0("concordant", type)
  dis_var <- paste0("discordant", type)
  
  res <- res[res$epval < pval & res$ppval < pval, ]
  
  if (type != "_pri") {
    con_count <- nrow(res[res[[paste0("type", type)]] == "concordant", ])
    dis_count <- nrow(res[res[[paste0("type", type)]] == "discordant", ])
  } else {
    res <- res[!is.na(res$type_primary), ]
    con_count <- nrow(res[res$type_primary == "concordant", ])
    dis_count <- nrow(res[res$type_primary == "discordant", ])
  }
  
  if (type == "") {
    num <- 2
  } else if (type == "_z") {
    num <- 4
  } else {
    num <- 6
  }
  
  # P value calculation
  dat$pval <- apply(dat, 1, function(x) {
    fdf <- data.frame(
      "concord" = c(as.numeric(x[num]), con_count - as.numeric(x[num])),
      "discord" = c(as.numeric(x[num+1]), dis_count - as.numeric(x[num+1])),
      row.names = c("linked", "not"),
      stringsAsFactors = F
    )
    
    p <- try(fisher.test(fdf)$p.value)
    if (inherits(p, "error")) {
      return(1)
    }
    return(p)
  })
  
  # Plotting
  dat$con_percent <- dat[[num]] / con_count * 100
  dat$dis_percent <- dat[[num+1]] / dis_count * 100
  dat$diff_percent <- abs(dat$con_percent - dat$dis_percent)
  dat$diff <- abs(dat[num] - dat[num+1])
  dat$pval <- format(dat$pval, digits = 3)
  
  #res.m <- reshape2::melt(dat[dat$diff_percent >= 1, c(term, paste0("concordant", type), paste0("discordant", type), "pval")], id.vars = c(term, "pval"))
  #res.m2 <- reshape2::melt(dat[(dat[[paste0("concordant", type)]] >= 5 | dat[[paste0("discordant", type)]] >= 5) & dat$pval < 0.05, c(term, "con_percent", "dis_percent", "pval")], id.vars = c(term, "pval"))
  res.m2 <- reshape2::melt(dat[(dat[[paste0("concordant", type)]] >= 1 | dat[[paste0("discordant", type)]] >= 1), c(term, "con_percent", "dis_percent", "pval")], id.vars = c(term, "pval"))
  res.m2 <- reshape2::melt(dat[, c(term, "con_percent", "dis_percent", "pval")], id.vars = c(term, "pval"))
  
  pvals <- unique(res.m2$pval)
  fdrs <- prettyNum(p.adjust(pvals, method = "fdr"), digits = 2)
  pvals_to_merge <- as.data.frame(cbind(pvals, fdrs))
  names(pvals_to_merge) <- c("pval", "fdr")
  res.m2 <- base::merge(res.m2, pvals_to_merge, by = "pval", all.x = T)
  
  # Prettyify
  #res.m2$pval <- prettyNum(as.double(res.m2$pval), digits = 2)
  #res.m2$fdr <- prettyNum(as.double(res.m2$fdr), digits = 2)
  res.m2$p_label <- ifelse(res.m2$variable == "con_percent" | (as.double(res.m2$pval) > 0.05 & as.double(res.m2$fdr) > 0.05),
                           "",
                           paste0("P = ", prettyNum(as.double(res.m2$pval), digits = 2), "; FDR = ", prettyNum(as.double(res.m2$fdr), digits = 2)))
  #res.m2[res.m2$variable == "con_percent", "p_label"] <- ""
  #res.m2[res.m2$pval > 0.05 & res.m2$fdr > 0.05, "p_label"] <- ""
  
  terms <- as.data.frame(unique(res.m2[order(res.m2$term), ]$term))
  res.m2$term <- forcats::fct_rev(factor(res.m2$term))
  
  # For signif lines
  # First find the x locations for the bars
  signif_x <- terms %>%
    dplyr::mutate(x_min = 1 * dplyr::row_number()) %>%
    dplyr::mutate(x_max = x_min + 0.3)
  names(signif_x) <- c("term", "x_min", "x_max")
  
  # Now find the y locations
  signif_y <- res.m2 %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(y_end = max(value) + 1) %>%
    dplyr::select(term, y_end) %>%
    dplyr::distinct(term, y_end)
  
  signif <- base::merge(signif_x, signif_y, by = "term")
  signif$x_min <- rev(signif$x_min)
  signif$x_max <- rev(signif$x_max)
  
  signif <- base::merge(signif, res.m2[c("term", "p_label", "fdr")], by = "term", all = T)
  
  title <- paste0("Enrichment analysis using DGIdb terms.")
  #subtitle <- paste0("QTL pair selection threshold: P < ", pval)
  #legend <- paste0("QTL pair category", ifelse(type == "_pri", "\n(Primary signals only)", ""))
  legend <- "QTL pair category"
  
  p <- ggplot(data = res.m2, aes_string(x = term, y = "value")) +
    geom_bar(aes(fill = variable), position = "dodge", stat = "identity") +
    geom_text(aes(label = p_label), position=position_dodge(width=1), hjust = -0.25, vjust = -0.25) +
    scale_fill_manual(legend, 
                      #labels = c(paste0("Concordant (n = ", con_count, ")"), paste0("Discordant (n = ", dis_count, ")")), 
                      labels = c("Concordant", "Discordant"),
                      values = c("con_percent" = concord_colour, "dis_percent" = discord_colour)) +
    #geom_signif(y = signif$y_end, y_end = signif$y_end,
    #            xmin = signif$x_min, xmax = signif$x_max,
    #            annotation = signif$pval) +
    #theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(title) +
    #scale_x_discrete(limits = rev(levels(res.m2[[term]]))) +
    scale_y_continuous(limits = c(0, 70), expand = expansion(0, c(0, 0))) +
    coord_flip() +
    labs(y = "Percentage of QTL pairs enriched",
         x = "DGIdb term",
         title = title
         #subtitle = subtitle
    ) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          #plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 16))
  
  return(p)
}

p1 <- plot_results(res_dgidb, res, term = "term", pval = dgidb_pval)
p2 <- plot_results(res_dgidb, res, term = "term", type = "_pri", pval = dgidb_pval)

p <- cowplot::plot_grid(
  p1 + 
    labs(title = "i) Primary QTL pairs only") +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 12)),
  p2 + 
    labs(title = "ii) All QTL pairs") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 12)),
  align = "h",
  ncol = 2,
  rel_widths = c(1, 0.85)
)
title.grob <- textGrob("Enrichment analysis using DGIdb terms", gp=gpar(fontsize=22))
x.grob <- textGrob("Percentage of QTL pairs enriched", gp=gpar(fontsize=16))
y.grob <- textGrob("DGIdb term", gp=gpar(fontsize=16), rot=90)
grid.arrange(arrangeGrob(p, top = title.grob, bottom = x.grob, left = y.grob))
