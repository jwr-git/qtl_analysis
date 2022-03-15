library(metafor)
gtex <- read.table("../results/pqtls_gtex_harmonised.txt",
                   sep = "\t", header = T, stringsAsFactors = F, quote = "")
gtex <- gtex[gtex$pval.exposure < 2.95e-5 & gtex$pval.outcome < 2.95e-5, ]

gtex <- gtex %>%
  tidyr::separate(outcome, c("tissue", "gene", "rs"), sep = "\\.")

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

#######
# Supp Table 2
prop_res <- data.frame(
  t <- character(),
  pairs <- double(),
  expected.p <- double(),
  expected.sign <- double(),
  observed.p <- double(),
  observed.sign <- double()
)

for (t in unique(gtex$tissue)) {
  dat <- gtex[gtex$tissue == t, ]
  prop_res_gtex <- prop_overlap(dat$beta.exposure, dat$beta.outcome,
                                dat$se.exposure, dat$se.outcome,
                                0.05)
  
  r <- c(t, 
         prop_res_gtex$res$nsnp[1],
         prop_res_gtex$res$value[3],
         prop_res_gtex$res$value[1],
         prop_res_gtex$res$value[4],
         prop_res_gtex$res$value[2])
  
  prop_res <- rbind(prop_res, r)
}
names(prop_res) <- c("tissue", "npairs", "expected.p", "expected.sign", "observed.p", "observed.sign")
write.table(prop_res, "../results/prop_res_gtex.txt",
            sep = "\t", row.names = F, quote = F)
