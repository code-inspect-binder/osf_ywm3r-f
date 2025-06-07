#@@@@@@@@@@@@@#
## Readme ----
#@@@@@@@@@@@@@#

# Here we provide data (with codebook) and code for reproducing the analyses of
# Costantini, G., Saraulli, D., Perugini, M. (in press), Uncovering the 
# Motivational Core of Traits: The Case of Conscientiousness.
# European Journal of Personality. https://doi.org/10.1002/per.2237 


#@@@@@@@@@@@@@@@@@@@@#
## load packages ----
#@@@@@@@@@@@@@@@@@@@@#

# load R packages using the pacman tool
if(!require(pacman)) install.packages("pacman")
pacman::p_load("dplyr", "psych", "GPArotation",
               "qgraph", "bootnet", "corpcor",
               "lavaan", "semTools",
               "Partiallyoverlapping",
               "lm.beta", "AutoModel",
               "censReg", "huge")

#@@@@@@@@@@@@@@#
## Study 1 ----
#@@@@@@@@@@@@@@#

# load data
dt1 <- read.csv("Study1.csv")

#### Hypothesis 1.1 ####
# Hypothesis 1 was tested via paired-sample t.tests
# (alphas need to be Bonferroni-corrected)

# Examples for conscientious goal classes (first two rows of
# Table 1 - Conscientious Goals). Remaining tests are simialr
t.test(x = dt1$CP_Class_01, y = dt1$CN_Class_01, paired = TRUE,
       alternative = "greater")
t.test(x = dt1$CP_Class_02, y = dt1$CN_Class_02, paired = TRUE,
       alternative = "greater")

# Examples for unconscientious goal classes (first two rows of
# Table 1 - Unconscientious Goals)
t.test(x = dt1$CN_Class_04, y = dt1$CP_Class_04, paired = TRUE,
       alternative = "greater")
t.test(x = dt1$CN_Class_06, y = dt1$CP_Class_06, paired = TRUE,
       alternative = "greater")

#### Hypothesis 1.2 ####
# Hypothesis 1.2 was tested via partially-overlapping samples
# t-tests. We report the analyses for the first panel of Table
# S2 (and for p-values reported in Table 2). The remaining
# analyses are identical, but were performed on other goal classes
Partover.test(dt1$CP_Class_01, dt1$HP_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$EP_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$XP_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$AP_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$OP_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$HN_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$EN_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$XN_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$AN_Class_01,
              stacked = TRUE, alternative = "greater")
Partover.test(dt1$CP_Class_01, dt1$ON_Class_01,
              stacked = TRUE, alternative = "greater")

#@@@@@@@@@@@@@@#
## Study 2 ----
#@@@@@@@@@@@@@@#

dt2 <- read.csv("Study2.csv")

#### Hypothesis 2.1 ####
# correlations between conscientiousness and goal importance
select(dt2, HEXACO_C, Importance_GC, Importance_GU) %>%
  corr.test()

# overall goal importance regressed on HEXACO traits
lm(Importance_GC ~ ., data = select(dt2, Importance_GC,
                                    HEXACO_H:HEXACO_O)) %>%
     lm.beta %>% summary

lm(Importance_GU ~ ., data = select(dt2, Importance_GU,
                                    HEXACO_H:HEXACO_O)) %>%
  lm.beta %>% summary

# Examples of Tobit regressions of the subjective importance of
# specific goals on HEXACO traits (first two lines of Table 3).
# P-values need to be corrected using the Bonferroni method,
# considering that we performed 21 multiple regressions
# (multiply them by 21 and round values larger than 1 to 1)

predictors <- select(dt2, HEXACO_H:HEXACO_O) %>% scale()
fit1 <- censReg(dt2$G_C_01_DimostrareQlcQln_Importance ~ predictors,
        left = 1, right = 9)
summary(fit1)

fit2 <- censReg(dt2$G_C_02_EssereDegnoFiducia_Importance ~ predictors,
        left = 1, right = 9)
summary(fit2)


# A function to compute the McKelvey & Zavoina R2
# following McBee (2010) and Veall & Zimmermann (1994)
McKelveyZavoinaR2 <- function(fit, IVs, DV, left, right)
{
  # fit: a fitted censReg object
  # IVs: predictor matrix
  # DV: dependent variable
  # left, right: limits for the censored dependent variable (see censReg)
  
  # compute continuous predicted variables (y*)
  prd <- cbind(1, IVs, 0) %*%  fit$estimate
  # compute truncated predicted variables
  prdcut <- prd
  prdcut[prdcut > right] <- right
  prdcut[prdcut < left] <- left
  # compute errors
  err <- (DV - prdcut)
  # compute R2
  R2 <- var(prd, na.rm = T)/(var(prd, na.rm = T) + var(err, na.rm = T))
  R2
}  


McKelveyZavoinaR2(fit1, IVs = predictors, DV = dt2$G_C_01_DimostrareQlcQln_Importance,
                   left = 1, right = 9)
McKelveyZavoinaR2(fit2, IVs = predictors, DV = dt2$G_C_02_EssereDegnoFiducia_Importance,
                   left = 1, right = 9)

#### Hypothesis 2.2 ####
# Most individuals were willing to become more conscientious
# according to BFTGI
table(dt2$BFTGI_C)
t.test(dt2$BFTGI_C, mu = 2)
# and according to CBFI
t.test(dt2$CBFI_C, mu = 3)

# Correlations of the willingness to change conscientiousness
# (C-BFI)
corr.test(select(dt2, CBFI_C),
          select(dt2, HEXACO_C, BFI2_C, POS, Importance_GC,
       Importance_GU))

# Separate logistic regressions of the willingness to increase
# conscientiousness according to BF-TGI
glm(I(BFTGI_C == 3) ~ scale(HEXACO_C), family = "binomial",
    data = dt2) %>% summary

glm(I(BFTGI_C == 3) ~ scale(BFI2_C), family = "binomial",
    data = dt2) %>% summary

glm(I(BFTGI_C == 3) ~ scale(POS), family = "binomial",
    data = dt2) %>% summary

glm(I(BFTGI_C == 3) ~ scale(Importance_GC), family = "binomial",
    data = dt2) %>% summary

glm(I(BFTGI_C == 3) ~ scale(Importance_GU), family = "binomial",
    data = dt2) %>% summary

# Hierarchical regression predicting the willingness to change
# conscientiousness according to CBFI (Table 4). For ease
# of formatting, I used the R package AutoModel
Data4HierarchicalReg <- select(dt2,
                               CBFI_C, POS,
                               HEXACO_H:HEXACO_O,
                               BFI2_O:BFI2_N,
                               Importance_GC, Importance_GU) %>%
  scale() %>%
  data.frame()
Data4HierarchicalReg$BFTGI_C_bin <- as.numeric(dt2$BFTGI_C == 3)

# fit the model
invisible(capture.output(
  mdls <- run_model(outcome = "CBFI_C",
                    c("HEXACO_H", "HEXACO_E", "HEXACO_X",
                      "HEXACO_A", "HEXACO_O"),
                    "HEXACO_C",
                    "POS",
                    c("Importance_GC", "Importance_GU"),
                    dataset = Data4HierarchicalReg,
                    assumptions.check = FALSE,
                    type = "gaussian")
))

# format the results
mdls$SummaryDF[,c("R2", "DeltaR2", "SigFch", "pval") ]
Tbl4_CBFI <- mdls$CoefficientsDF[c(1:6,13, 21, 30:31),
                            c("Model", "term", "estimate", "sig")]
Tbl4_CBFI[,"estimate"] <- round(Tbl4_CBFI[,"estimate"], 2)
Tbl4_CBFI

# Hierarchical regression predicting the willingness to change
# conscientiousness according to BFTGI (Table 4).

# fit model
invisible(capture.output(
mdls <- run_model(outcome = "BFTGI_C_bin",
                  c("HEXACO_H", "HEXACO_E", "HEXACO_X",
                    "HEXACO_A", "HEXACO_O"),
                  "HEXACO_C",
                  "POS",
                  c("Importance_GC", "Importance_GU"),
                  dataset = Data4HierarchicalReg,
                  assumptions.check = FALSE,
                  type = "binomial")
))

frm <- create_formula_objects(outcome = "BFTGI_C_bin",
                              c("HEXACO_H", "HEXACO_E",
                                "HEXACO_X",
                                "HEXACO_A", "HEXACO_O"),
                              "HEXACO_C",
                              "POS",
                              c("Importance_GC",
                                "Importance_GU"))
mdl <- create_model_objects(frm, type = "binomial",
                            dataset = Data4HierarchicalReg)

# Coefficients
bind_rows(mdls$Coefficients[1:6,c(1, 2, 5)],
      mdls$Coefficients[7,c(1, 6, 9)],
      mdls$Coefficients[8,c(1, 10, 13)],
      mdls$Coefficients[9:10,c(1, 14, 17)])

# Model comparisons
anova(mdl[[1]], mdl[[2]], mdl[[3]], mdl[[4]], 
      test = "Chisq")
# Nagelkerke R2
mdls$Summary[,"Nagelkerke"]

# LRT p-value for first model
1-pchisq(mdl[[1]]$deviance, mdl[[1]]$df.residual)
# LRT p-value for full model
1-pchisq(mdl[[4]]$deviance, mdl[[4]]$df.residual)


#@@@@@@@@@@@@@@#
## Study 3 ----
#@@@@@@@@@@@@@@#

dt3 <- read.csv("Study3.csv")

#### Hypothesis 3.1 ####
# factor models for lavaan
fmdl <- "
f8 =~ Goal_8_Agire_secondo_le_regole +
 Goal_8_Rispettare_le_regole +
 Goal_8_Evitare_di_infrangere_le_regole +
 Goal_8_Seguire_le_leggi +
 Goal_8_Rispettare_il_contesto_in_cui_mi_trovo +
 Goal_8_Rispettare_un.autorita +
 Goal_8_Non_uscire_dagli_schemi +
 Goal_8_Agire_secondo_i_miei_valori
f10 =~ Goal_10_Avere_il_pieno_controllo +
 Goal_10_Avere_tutto_sotto_controllo +
 Goal_10_Avere_la_situazione_sotto_controllo +
 Goal_10_Evitare_il_caos +
 Goal_10_Non_fare_confusione +
 Goal_10_Tenere_tutto_sotto_controllo +
 Goal_10_Ritrovare_le_cose_quando_ne_ho_bisogno +
 Goal_10_Trovare_le_mie_cose
f11 =~ Goal_11_Ottenere_buoni_risultati +
 Goal_11_Superare_me_stesso +
 Goal_11_Avere_successo +
 Goal_11_Raggiungere_risultati_scolastici_accademici_o_lavorativi +
 Goal_11_Realizzarmi + 
 Goal_11_Avere_un_futuro +
 Goal_11_Avere_sicurezza_economica +
 Goal_11_Ottenere_il_massimo
f12  =~ Goal_12_Fare_le_cose_per_bene +
 Goal_12_Fare_al_meglio_cio_che_faccio +
 Goal_12_Fare_un_ottimo_lavoro +
 Goal_12_Fare_le_cose_nel_migliore_dei_modi +
 Goal_12_Evitare_di_fare_errori +
 Goal_12_Non_sbagliare +
 Goal_12_Non_tralasciare_alcun_dettaglio +
 Goal_12_Non_commettere_sviste
f13 =~ Goal_13_Portare_a_termine_qualcosa +
 Goal_13_Portare_a_termine_un_progetto +
 Goal_13_Terminare_un_lavoro_in_tempo +
 Goal_13_Fare_tutto_nei_tempi_prestabiliti +
 Goal_13_Mantenere_gli_impegni_presi +
 Goal_13_Non_lasciare_le_cose_a_meta +
 Goal_13_Non_lasciare_un_lavoro_a_meta +
 Goal_13_Rispettare_un_impegno_preso
f17 =~ Goal_17_Non_correre_rischi + 
 Goal_17_Non_mettermi_nei_guai +
 Goal_17_Non_mettermi_in_pericolo +
 Goal_17_Non_farmi_del_male +
 Goal_17_Proteggermi +
 Goal_17_Stare_al_sicuro +
 Goal_17_Stare_in_salute +
 Goal_17_Salvaguardarmi 
f25 =~ Goal_25_Analizzare_bene_le_situazioni +
 Goal_25_Riordinare_le_idee +
 Goal_25_Prevedere_le_conseguenze_delle_mie_azioni +
 Goal_25_Valutare_tutte_le_opzioni_prima_di_prendere_una_decisione +
 Goal_25_Prendere_buone_decisioni +
 Goal_25_Prendere_decisioni_giuste +
 Goal_25_Trovare_la_soluzione_migliore_quando_ho_un_problema +
 Goal_25_Non_lasciare_nulla_al_caso
f16 =~ Goal_16_Evitare_di_fare_qualcosa +
 Goal_16_Evitare_un_compito_faticoso +
 Goal_16_Evitare_un_impegno +
 Goal_16_Far_vedere_che_qualcosa_non_mi_interessa +
 Goal_16_Far_capire_a_qualcuno_che_non_puo_contare_su_di_me +
 Goal_16_Non_assumermi_responsabilita +
 Goal_16_Non_fare_cose_che_non_ho_voglia_di_fare +
 Goal_16_Dedicare_il_mio_tempo_solo_a_cio_che_mi_piace_fare
f26 =~ Goal_26_Allontanare_pensieri_problematici +
 Goal_26_Non_pensare +
 Goal_26_Non_ragionare_troppo +
 Goal_26_Distrarmi_dai_problemi +
 Goal_26_Essere_meno_razionale +
 Goal_26_Evitare_di_pensare_troppo +
 Goal_26_Non_pensare_alle_conseguenze_di_cio_che_faccio +
 Goal_26_Non_preoccuparmi_di_cio_che_faccio"
OneSecondOrder <- "
GOAL =~ f8 +  f10 + f11 + f12 + f13 + f17 + f25 + f16 + f26"
TwoSecondOrder <- "
CGOAL =~ f8 +  f10 + f11 + f12 + f13 + f17 + f25 
UGOAL =~ f16 + f26"

fits <- list()
# first order oblique
fits[[1]] <- lavaan::cfa(model = fmdl, data = dt3, orthogonal = FALSE)
# Two second-order, correlated
fits[[2]] <- lavaan::cfa(model = paste(fmdl, ";", TwoSecondOrder), data = dt3)
# One second-order 
fits[[3]] <- lavaan::cfa(model = paste(fmdl, ";", OneSecondOrder), data = dt3)
# first order orthogonal
fits[[4]] <- lavaan::cfa(model = fmdl, data = dt3, orthogonal = TRUE)

# baseline RMSEA is too small for interpreting CFI
semTools::moreFitIndices(fits[[1]], "baseline.rmsea")

# Fit indices and LRT (presented in Table 6 in the manuscript)
sapply(fits, fitmeasures, c("chisq", "df", "pvalue",
                            "rmsea",
                            "rmsea.ci.lower",
                            "rmsea.ci.upper",
                            "rmsea.pvalue",
                            "srmr",
                            "AIC",
                            "BIC")) %>% t()

anova(fits[[1]], fits[[2]], fits[[3]], fits[[4]])

# CFA Coefficients for the final model (presented in Table 5 and 7
# in the manuscript)
summary(fits[[1]], std = TRUE)

#### Hypothesis 3.2 ####
# multiple regressions predicting each goal class from HEXACO
# traits (we report examples reproducing the first two rows
# of Table 8, the code for the others is similar, save for the
# goal class). P-values need to be corrected using the
# Bonferroni method, considering 9 multiple regressions
# (multiply them by 9 and transform values > 1 to 1)
lm(G08Rules ~ ., data = select(dt3, G08Rules,
                               HEXACO_H:HEXACO_O)) %>%
     lm.beta %>% summary
  
lm(G10Control ~ ., data = select(dt3, G10Control,
                               HEXACO_H:HEXACO_O)) %>%
  lm.beta %>% summary



#### Hypothesis 3.3 ####
# Compute facet scores (see Materials section -
# Conscientiousness facets)
dt3$IMC <- select(dt3, IPIP_C_ca, ACC_ci, HEXACOC_PRU) %>% 
  scale() %>% rowMeans(na.rm = TRUE)
dt3$IND <- select(dt3, IPIP_C_ac, ACC_tp, HEXACOC_DIL) %>%
  scale() %>% rowMeans(na.rm = TRUE)
dt3$ORD <- select(dt3, IPIP_C_or, ACC_od, HEXACOC_ORD) %>%
  scale() %>% rowMeans(na.rm = TRUE)

# select data for the network
dtNet <- select(dt3, IMC, IND, ORD, SCS, CFC, FU,
                G08 = G08Rules,
                G10 = G10Control,
                G11 = G11Realiz,
                G12 = G12DoWell,
                G13 = G13Accomplish,
                G17 = G17Safety,
                G25 = G25Think,
                G16 = G16Avoid,
                G26 = G26DoNotThink)

# apply a nonparanormal transformation
dtNet <- huge.npn(dtNet)

# estimate the network using EBICglasso
net <- EBICglasso(cor(dtNet), n = nrow(dtNet))
# deifne colors
clr <- c(rep("lightblue", 3), rep("pink", 3),
         rep("palegreen", 7), rep("firebrick1", 2))
# estimate predictability using a function available in the 
# "R" folder
source("predictability.R")
prd <- predictability(net, dt = dtNet)
# visualize the network with qgraph
qg <- qgraph(net, layout = "spring", color = clr, pie = prd,
             pieColor = "black", pieBorder = .20,
             negDashed = TRUE)

# density
sum(net[lower.tri(net)] != 0) / (nrow(net)*(ncol(net)-1)/2)


# CS stability coefficient for predictability
# (be advised: bootnet takes a long time)
# first preform case-dropping bootstrap using custom
# parameters (drop from 5% to 95% of the samples, in steps
# of 5%, 200 replications for each bootstrapped sample size)
fun <- function(x) EBICglasso(cor(x), n = nrow(x))
suppressWarnings(
  suppressMessages(
    {set.seed(1);
      btnt <- bootnet(dtNet,
                      type = "case",
                      statistics = c("strength"),
                      prepFun = identity,
                      prepArgs = NULL,
                      estFun = fun,
                      estArgs = NULL,
                      nBoots = 18200,
                      caseMin = .05,
                      caseMax = .95,
                      caseN = 91,
                      verbose = FALSE,
                      nCores = 1)}
  ))

# compute predictability in the bootstrapped networks
bootNets <- lapply(btnt$boots, function(x) x$graph)
bootPred <- lapply(bootNets, predictability, dtNet)

# compute correlations between original and boostraped
#  predictabilities
cors <- t(sapply(bootPred, cor, prd, use = "p"))
nBoot <- sapply(btnt$boots, function(x) x$nPerson)
n <- btnt$sampleSize

# compute CS index as the proportion of sample that can be
# dropped and still have a correlation > .70 between the
# predictability on the original and on the bootstrapped
# sample in 95% of cases
myboot <- data.frame("prop" = 1-nBoot/n, "cors" = cors[1,])
prop <- group_by(myboot, prop) %>%
  summarize(P = mean(cors > .7))
max(prop[prop[,"P"] > .95, "prop"])


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
## Output of SessionInfo() for reproducibility of results ----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 17763)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=Italian_Italy.1252  LC_CTYPE=Italian_Italy.1252   
# [3] LC_MONETARY=Italian_Italy.1252 LC_NUMERIC=C                  
# [5] LC_TIME=Italian_Italy.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] huge_1.3.3               censReg_0.5-26           maxLik_1.3-6            
# [4] miscTools_0.6-22         AutoModel_0.4.9          lm.beta_1.5-1           
# [7] Partiallyoverlapping_2.0 semTools_0.5-2           lavaan_0.6-5            
# [10] corpcor_1.6.9            bootnet_1.2.4            ggplot2_3.2.1           
# [13] qgraph_1.6.3             GPArotation_2014.11-1    psych_1.8.12            
# [16] dplyr_0.8.3              pacman_0.5.1            
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1         backports_1.1.4      Hmisc_4.2-0          BDgraph_2.61        
# [5] plyr_1.8.4           igraph_1.2.4.1       lazyeval_0.2.2       splines_3.6.1       
# [9] candisc_0.8-0        digest_0.6.21        foreach_1.4.7        htmltools_0.3.6     
# [13] matrixcalc_1.0-3     gdata_2.18.0         magrittr_1.5         checkmate_1.9.4     
# [17] cluster_2.1.0        doParallel_1.0.15    ROCR_1.0-7           aod_1.3.1           
# [21] etm_1.0.5            openxlsx_4.1.0.1     longitudinal_1.1.12  wordcloud_2.6       
# [25] R.utils_2.9.0        sandwich_2.5-1       bdsmatrix_1.3-3      jpeg_0.1-8          
# [29] colorspace_1.4-1     mitools_2.4          haven_2.1.1          pan_1.6             
# [33] xfun_0.9             crayon_1.3.4         networktools_1.2.1   lme4_1.1-21         
# [37] rowr_1.1.3           zeallot_0.1.0        survival_2.44-1.1    zoo_1.8-6           
# [41] iterators_1.0.12     glue_1.3.1           relaimpo_2.2-3       gtable_0.3.0        
# [45] nnls_1.4             NetworkToolbox_1.3.0 car_3.0-3            weights_1.0         
# [49] ggm_2.3              jomo_2.6-9           abind_1.4-5          scales_1.0.0        
# [53] mvtnorm_1.0-11       DBI_1.0.0            bibtex_0.4.2         Rcpp_1.0.2          
# [57] plotrix_3.7-6        cmprsk_2.2-8         mgm_1.2-7            htmlTable_1.13.2    
# [61] foreign_0.8-71       Formula_1.2-3        stats4_3.6.1         heplots_1.3-5       
# [65] survey_3.36          glmnet_2.0-18        htmlwidgets_1.3      gplots_3.0.1.1      
# [69] RColorBrewer_1.1-2   acepack_1.4.1        IsingFit_0.3.1       mice_3.6.0          
# [73] pkgconfig_2.0.3      R.methodsS3_1.7.1    plm_2.2-0            nnet_7.3-12         
# [77] tidyselect_0.2.5     rlang_0.4.0          reshape2_1.4.3       polynom_1.4-0       
# [81] munsell_0.5.0        cellranger_1.1.0     tools_3.6.1          generics_0.0.2      
# [85] IsingSampler_0.2     broom_0.5.2          glmmML_1.1.0         fdrtool_1.2.15      
# [89] stringr_1.4.0        knitr_1.25           zip_2.0.4            caTools_1.17.1.2    
# [93] purrr_0.3.2          mitml_0.3-7          GeneNet_1.2.13       glasso_1.10         
# [97] pbapply_1.4-2        nlme_3.1-140         whisker_0.4          R.oo_1.22.0         
# [101] smacof_2.0-0         BaylorEdPsych_0.5    compiler_3.6.1       rstudioapi_0.10     
# [105] curl_4.1             png_0.1-7            tibble_2.1.3         pbivnorm_0.6.0      
# [109] stringi_1.4.3        Epi_2.38             forcats_0.4.0        eigenmodel_1.11     
# [113] lattice_0.20-38      Matrix_1.2-17        nloptr_1.2.1         vctrs_0.2.0         
# [117] parcor_0.2-6         pillar_1.4.2         lifecycle_0.1.0      lmtest_0.9-37       
# [121] Rdpack_0.11-0        bitops_1.0-6         data.table_1.12.2    gbRd_0.4-11         
# [125] R6_2.4.0             latticeExtra_0.6-28  KernSmooth_2.23-15   gridExtra_2.3       
# [129] rio_0.5.16           codetools_0.2-16     boot_1.3-22          MASS_7.3-51.4       
# [133] gtools_3.8.1         assertthat_0.2.1     rjson_0.2.20         withr_2.1.2         
# [137] mnormt_1.5-5         ppls_1.6-1.1         mgcv_1.8-28          parallel_3.6.1      
# [141] hms_0.5.1            grid_3.6.1           rpart_4.1-15         tidyr_1.0.0         
# [145] minqa_1.2.4          carData_3.0-2        d3Network_0.5.2.1    numDeriv_2016.8-1.1 
# [149] base64enc_0.1-3      ellipse_0.4.1       