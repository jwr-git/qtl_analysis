library(metafor)
gtex <- read.table("../results/pqtls_gtex_harmonised.txt",
                   sep = "\t", header = T, stringsAsFactors = F, quote = "")
gtex <- gtex[gtex$pval.exposure < 2.95e-5 & gtex$pval.outcome < 2.95e-5, ]

gtex <- gtex %>%
  tidyr::separate(outcome, c("tissue", "gene", "rs"), sep = "\\.")

#######
# Supp Table 5
hetero_res <- data.frame(
  term = character(),
  QE = double(),
  QEp = double(),
  k = numeric(),
  I2 = double(),
  H2 = double()
)
#for (t in unique(hetero$tissue)) {
#  hetdat <- hetero[hetero$tissue == t,]
#  res <- rma(hetdat$beta.outcome, hetdat$samplesize.outcome, method="FE")
#  temp <- as.data.frame(c(t, res[c("QE", "QEp", "k", "I2", "H2")]))
#  names(temp) <- c("term", "QE", "QEp", "k", "I2", "H2")
#  hetero_res <- dplyr::bind_rows(hetero_res, temp)
#}
for (snp in unique(gtex$rs)) {
  hetdat <- gtex[gtex$rs == snp,]
  res <- rma(hetdat$beta.outcome, hetdat$samplesize.outcome, method="FE")
  temp <- as.data.frame(c(snp, res[c("QE", "QEp", "k", "I2", "H2")])) %>%
    dplyr::mutate(k = k - 1) %>%
    dplyr::filter(k > 1)
  names(temp) <- c("term", "QE", "QEp", "k", "I2", "H2")
  hetero_res <- dplyr::bind_rows(hetero_res, temp)
}

write.table(hetero_res, 
            "../results/hetero_gtex_res.txt",
            sep = "\t", row.names = F, quote = F)
