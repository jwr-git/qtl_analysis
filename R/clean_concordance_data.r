library(httr)
library(readxl)
library(ggplot2)
library(dplyr)
library(ieugwasr)
library(biomaRt)

rm(list = ls(all.names = TRUE))

grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                 host="grch37.ensembl.org", 
                 path="/biomart/martservice", 
                 dataset="hsapiens_gene_ensembl")

################################################################################
# Step 1
# Load, clean and combine the Fenland and Zhang datasets
zhang <- readxl::read_xlsx("../data/EUR-pQTL-Zhang.xlsx", 1)
fenland <- readxl::read_xlsx("../data/fenland_pqtl_instruments.xlsx", 1)

# Zhang dataset needs some more annotation
zhang$gene.name <- paste0("'", zhang$Full_sum_gene)
zhang <- zhang[!(zhang$gene.name %in% fenland$HGNC.symbol.protein), ]
zhang$cis.trans <- "cis"
zhang$study <- "Zhang"
table <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
               filters = "hgnc_symbol", 
               values = zhang$Full_sum_gene, 
               mart = grch37)
zhang <- merge(zhang, table, by.x = "Full_sum_gene", by.y = "hgnc_symbol", all.x = T)

zhang_less <- zhang[c("Conditional independent cis-SNP", "UniProt", "EA", "AA", "EAF", "BETA", "SE", "P", "N", "ensembl_gene_id", "gene.name", "cis.trans", "study")]
zhang_less <- zhang_less %>% dplyr::group_by(gene.name) %>% dplyr::arrange(P, .by_group = T) %>% dplyr::mutate(signal = row_number())

# Manual mapping for some targets
zhang_less[zhang_less$gene.name == "'ADGRE2", ]$ensembl_gene_id <- "ENSG00000127507"
zhang_less[zhang_less$gene.name == "'ADGRF5", ]$ensembl_gene_id <- "ENSG00000069122"
zhang_less[zhang_less$gene.name == "'DCI", ]$ensembl_gene_id <- "ENSG00000167969"
zhang_less[zhang_less$gene.name == "'G6B", ]$ensembl_gene_id <- "ENSG00000204420"
zhang_less[zhang_less$gene.name == "'HSPC159", ]$ensembl_gene_id <- "ENSG00000119862"
zhang_less[zhang_less$gene.name == "'PACAP", ]$ensembl_gene_id <- "ENSG00000141433"
zhang_less[zhang_less$gene.name == "'SELM", ]$ensembl_gene_id <- "ENSG00000198832"

fenland$n <- 10708
fenland$study <- "Fenland"

fenland_less <- fenland[c("rsID", "Target", "EA", "NEA", "EAF", "Effect", "SE", "P.value", "n", "ENSEMBL.id.protein", "HGNC.symbol.protein", "cis.trans", "study", "signal.locus")]
fenland_less <- fenland_less[fenland_less$cis.trans == "cis", ]

header <- c("rsID", "Target", "EA", "NEA", "EAF", "Effect", "SE", "P.value", "n", "ENSG", "HGNC.symbol.protein", "cis.trans", "study", "signal.locus")
colnames(zhang_less) <- header
colnames(fenland_less) <- header
dat <- rbind(fenland_less, zhang_less)
write.table(dat, "../data/fenland_zheng_pqtls.txt", row.names = F, quote = F, sep = "\t")

################################################################################
# Step 2
# Load, clean and combine the pQTLs with GTEx
dat <- read.table("../data/fenland_zheng_pqtls.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
gtex <- read.table("../data/pqtl_gtex_all_snps.txt", sep = "\t", header = F, stringsAsFactors = F)
names(gtex) <- c("ensg", "id", "tss_distance", "ma_samples", "ma_count", "maf", "pval", "slope", "slope_se", "file")

gtex_rsids_map <- read.table("../data/pqtls_linked_gtex_ids.txt", sep = "\t", header = F, stringsAsFactors = F)
names(gtex_rsids_map) <- c("id", "chr", "pos", "ref", "alt", "num_alt", "rsID", "b37")

samplesizes_file <- read.csv("../data/gtex_sample_sizes.csv", header = T, stringsAsFactors = F)
names(samplesizes_file) <- c("tissue", "id", "n", "ignore")
ssn <- samplesizes_file[c("tissue", "n")]

# GTEx data needs mapped to rsIDs
gtex_rsids_map <- gtex_rsids_map[c("id", "rsID", "ref", "alt")]
gtex <- base::merge(gtex, gtex_rsids_map, by = "id", all.x = T)
gtex$ensg <- sub("\\..*", "", gtex$ensg)
gtex$file <- sub("\\..*", "", gtex$file)

lookups <- paste0(dat$ENSG, ".", dat$rsID)
gtex$lookup <- paste0(gtex$ensg, ".", gtex$rsID)
gtex <- gtex[which(gtex$lookup %in% lookups), ]

gtex_res <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(gtex_res) <- c("ensg", "hgnc", "rsid", "tissue", "type", "ppval", "epval", "signal", "indel")
all_gtex_harmonised <- list()
g_list <- unique(gtex$ensg)
for (i in 1:length(g_list)) {
  g <- g_list[i]
  if (is.na(g)) {
    next
  }

  subset_dat <- dat[which(dat$ENSG %in% g), ]
  subset_dat <- subset_dat[subset_dat$cis.trans == "cis", ]
  rsids <- unique(subset_dat$rsID)
  
  if (nrow(subset_dat) < 1) {
    next
  }
  subset_gtex <- gtex[which(gtex$ensg %in% g), ]
  if (nrow(subset_gtex) < 1) {
    next
  }
  subset_gtex <- subset_gtex[!duplicated(subset_gtex), ]
  subset_gtex$uid <- paste0(subset_gtex$file, ".", subset_gtex$lookup)
  subset_gtex <- subset_gtex %>%
    tidyr::separate(id, c("chr", "pos", "a1", "a2", "b"), sep = "_")
  subset_gtex$chr <- gsub("^.{0,3}", "", subset_gtex$chr)
  subset_gtex <- base::merge(subset_gtex, ssn, by.x = "file", by.y = "tissue")

  hgnc <- unique(subset_dat$HGNC.symbol.protein)[[1]]
  subset_dat <- subset_dat %>% TwoSampleMR::format_data(., type = "exposure",
                                                        snp_col = "rsID", beta_col = "Effect",
                                                        se_col = "SE", eaf_col = "EAF",
                                                        effect_allele_col = "EA", other_allele_col = "NEA",
                                                        pval_col = "P.value", phenotype_col = "HGNC.symbol.protein",
                                                        samplesize_col = "n", info_col = "signal.locus")
  subset_gtex <- subset_gtex %>% TwoSampleMR::format_data(., type = "outcome",
                                                          snp_col = "rsID", beta_col = "slope",
                                                          se_col = "slope_se",
                                                          effect_allele_col = "alt", other_allele_col = "ref",
                                                          pval_col = "pval",
                                                          phenotype_col = "uid", id_col = "uid",
                                                          chr_col = "chr", pos_col = "pos",
                                                          samplesize_col ="n")
  
  
  
  harmonised <- TwoSampleMR::harmonise_data(subset_dat, subset_gtex)
  if (nrow(harmonised) < 1) {
    next
  }

  # Type
  harmonised$type <- ifelse(
    (harmonised$beta.exposure >= 0 & harmonised$beta.outcome >= 0) | (harmonised$beta.exposure < 0 & harmonised$beta.outcome < 0),
    "concordant",
    "discordant"
  )
  
  # Merge into harmonised df for saving
  if (!length(all_gtex_harmonised)) {
    all_gtex_harmonised <- harmonised
  } else {
    all_gtex_harmonised <- dplyr::bind_rows(all_gtex_harmonised, harmonised)
  }
  
  for (t in unique(harmonised$id.outcome)) {
    tis_dat <- harmonised[harmonised$id.outcome == t, ]
    
    if (all("discordant" == tis_dat$type)) {
      type <- "discordant"
    } else if (all("concordant" == tis_dat$type)) {
      type <- "concordant"
    } else {
      type <- "partially discordant"
    }
    
    rsid <- tis_dat$SNP
    tiss <- tis_dat$outcome
    ppval <- tis_dat$pval.exposure
    epval <- tis_dat$pval.outcome
    signal <- tis_dat$info.exposure
    indel <- ifelse(tis_dat$effect_allele.exposure == "D" | tis_dat$effect_allele.exposure == "I", 1, 0)
    
    gtex_res[nrow(gtex_res) + 1, ] <- c(g, hgnc, rsid, tiss, type, ppval, epval, signal, indel)
  }
}

write.table(gtex_res, file = "../results/pqtls_gtex_concord.txt",
            row.names = F, quote = F, sep = "\t")
write.table(all_gtex_harmonised, "../results/pqtls_gtex_harmonised.txt",
            row.names = F, quote = F, sep = "\t")

################################################################################
# Step 3
# Merge pQTLs with eQTLGen and run the concordance tests
ensg.ids <- unique(dat$ENSG) # dat from step 1
ao <- ieugwasr::gwasinfo()
res <- data.frame(ensg = character(),
                  hgnc = character(),
                  snp = character(),
                  indel = integer(),
                  ppval = double(),
                  epval = double(),
                  type = character(),
                  type_std = character(),
                  type_primary = character(),
                  type_std_primary = character()
)
all_harmonised <- list()

# Function to standardise the effect sizes.
# Originally included as some SNPs in some QTL datasets have different sample
# sizes, so we sought to standardise the effects across all SNPs.
# Results from this were unused in the final paper/analyses but are kept here
# for posterity.
std_effect <- function(snp_freq, b, se, n) {
  # make sure they are numeric numbers
  snpfreq = as.numeric(as.character(snp_freq))
  b = as.numeric(as.character(b))
  se = as.numeric(as.character(se))
  n = as.numeric(as.character(n))
  
  zscore = b / se
  b_p = zscore / sqrt(2*snp_freq*(1-snp_freq)*(n+zscore^2))
  se_p = 1 / sqrt(2*snp_freq*(1-snp_freq)*(n+zscore^2))
  
  return(list(b=b_p,se=se_p))
}

# Not the most efficient way of doing this -- but it works!
for (i in 1:length(ensg.ids)) {
  ensg <- ensg.ids[i]
  if (is.na(ensg) || ensg %in% res$ensg) {
    next
  }
  
  id <- paste0("eqtl-a-", ensg)
  
  if (!(id %in% ao$id)) {
    next
  }
  
  # Subset data
  pqtldat <- dat[grepl(ensg, dat$ENSG, fixed = T) & dat$cis.trans == "cis", ]
  if (nrow(pqtldat) < 1) {
    next
  }
  
  eqtldat <- TwoSampleMR::extract_outcome_data(pqtldat$rsID, id)
  
  if (is.null(eqtldat)) {
    next
  }
  if (nrow(pqtldat) < 1 || nrow(eqtldat) < 1 || all(is.na(pqtldat)) || all(is.na(eqtldat))) {
    next
  }
  
  hgnc <- unique(pqtldat$HGNC.symbol.protein)[[1]]
  pqtldat <- pqtldat %>% TwoSampleMR::format_data(., type = "exposure",
                                                  snp_col = "rsID", beta_col = "Effect",
                                                  se_col = "SE", eaf_col = "EAF",
                                                  effect_allele_col = "EA", other_allele_col = "NEA",
                                                  pval_col = "P.value", phenotype_col = "HGNC.symbol.protein",
                                                  samplesize_col = "n", info_col = "signal.locus")
  pqtldat$exposure <- hgnc
  harmonised <- TwoSampleMR::harmonise_data(pqtldat, eqtldat)
  
  # Standardised concordance check
  stds_exp <- std_effect(harmonised$eaf.exposure, harmonised$beta.exposure, harmonised$se.exposure, harmonised$samplesize.exposure)
  stds_out <- std_effect(harmonised$eaf.outcome, harmonised$beta.outcome, harmonised$se.outcome, harmonised$samplesize.outcome)
  harmonised$type_std <- ifelse(
    (stds_exp$b >= 0 & stds_out$b >= 0) | (stds_exp$b < 0 & stds_out$b < 0),
    "concordant",
    "discordant"
  )
  
  # Naive concordance check
  harmonised$type <- ifelse(
    #as.numeric(harmonised$pval.exposure) <= 0.05 & as.numeric(harmonised$pval.outcome) <= 0.05 & 
    (harmonised$beta.exposure >= 0 & harmonised$beta.outcome >= 0) | (harmonised$beta.exposure < 0 & harmonised$beta.outcome < 0),
    "concordant",
    "discordant"
  )
  
  # Merge into harmonised df for saving
  if (!length(all_harmonised)) {
    all_harmonised <- harmonised
  } else {
    all_harmonised <- dplyr::bind_rows(all_harmonised, harmonised)
  }
  
  # Treat each SNP as independent signal
  for (j in 1:nrow(harmonised)) {
    hdat <- harmonised[j, ]
    
    snp <- hdat$SNP
    indel <- ifelse(hdat$effect_allele.exposure == "D" | hdat$effect_allele.exposure == "I", 1, 0)
    ppval <- hdat$pval.exposure
    epval <- hdat$pval.outcome
    type <- hdat$type
    type_std <- hdat$type_std
    
    # Check if this is the primary signal
    if (identical(hdat[hdat$info.exposure == 1, ]$type, character(0))) {
      type_primary <- NA
    } else {
      type_primary <- hdat[hdat$info.exposure == 1, ]$type
    }
    
    if (identical(hdat[hdat$info.exposure == 1, ]$type_std, character(0))) {
      type_std_primary <- NA
    } else {
      type_std_primary <- hdat[hdat$info.exposure == 1, ]$type_std
    }
    
    res[nrow(res) + 1, ] <- c(ensg, unique(pqtldat$exposure)[[1]], snp, indel, ppval, epval, type, type_std, type_primary, type_std_primary)
  }
}

write.table(res, file = "../results/pqtls_eqtlgen.txt",
            row.names = F, quote = F, sep = "\t")
write.table(all_harmonised, "../results/pqtls_eqtlgen_harmonised.txt",
            row.names = F, quote = F, sep = "\t")

################################################################################
# Step 4
# Winner's curse analysis
dat <- read.table("../results/pqtls_eqtlgen_harmonised.txt",
                  sep = "\t", header = T, stringsAsFactors = F, quote = "")
dat <- dat[dat$pval.exposure < 2.95e-5 & dat$pval.outcome < 2.95e-5, ]
dat <- dat[dat$info.exposure == 1,]

# This is the winner's curse function
# A file is also bundled in this repo with only this code (and with an
# explanation to boot)
prop_overlap = function(b_disc, b_rep, se_disc, se_rep, alpha)
{
  p_sign <- pnorm(-abs(b_disc) / se_disc) * pnorm(-abs(b_disc) / se_rep) + ((1 - pnorm(-abs(b_disc) / se_disc)) * (1 - pnorm(-abs(b_disc) / se_rep)))
  p_sig <- pnorm(-abs(b_disc) / se_rep + qnorm(alpha / 2)) + (1 - pnorm(-abs(b_disc) / se_rep - qnorm(alpha / 2)))
  p_rep <- pnorm(abs(b_rep)/se_rep, lower.tail=FALSE)
  res <- dplyr::tibble(
    nsnp=length(b_disc),
    metric=c("Sign", "Sign", "P-value", "P-value"),
    datum=c("Expected", "Observed", "Expected", "Observed"),
    value=c(sum(p_sign, na.rm=TRUE), sum(sign(b_disc) == sign(b_rep)), sum(p_sig, na.rm=TRUE), sum(p_rep < alpha, na.rm=TRUE))
  )
  return(list(res=res, variants=dplyr::tibble(
    sig=p_sig, 
    sign=p_sign, 
    fdisc=b_disc^2/se_disc^2, 
    frep=b_rep^2/se_rep^2, 
    p_rep=p_rep, 
    b_disc,
    b_rep,
    same_sign=sign(b_disc) == sign(b_rep))))
}

prop_res <- prop_overlap(dat$beta.exposure, dat$beta.outcome,
                         dat$se.exposure, dat$se.outcome,
                         0.05)
prop_res

norep <- prop_res$variants$sig > 0.99 & prop_res$variants$p_rep > 0.1

wrongdir <- (! prop_res$variants$same_sign) & 
  prop_res$variants$p_rep < 2.95e-5 & 
  prop_res$variants$sig > 0.99 &
  prop_res$variants$sign > 0.99
table(wrongdir)

discdir <- dat[wrongdir,] %>% dplyr::select(SNP, exposure, outcome, beta.exposure, beta.outcome, pval.outcome, type)
