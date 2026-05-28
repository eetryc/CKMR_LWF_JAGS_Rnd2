#=======================#
#==== 2.CKMRsim =====####
#=======================#

#=====Emily Tryc=====#
#====May 28, 2026====#



## >Load libraries ========================================================####

library(vcfR)
library(tidyverse)
library(DBI)
library(odbc)
library(dbplyr)
library(adegenet)
library(CKMRsim)
 # with so many packages make sure to [package]:: if needed, queen

## >Simulated pairs =======================================================####
# This is if you are testing the POP/HSP probability NOT true pairs! 
# Unnecesary when using actual genotyping data ! Skip to [] if using gt data.

# Only run one of the following two assignment functions!

  ### Random assignment  ==================================================####
  set.seed(777)
  kinships <- kinships %>%
    mutate(TruePairs = case_when(RObase == 2 ~ rbinom(nrow(.), size = 1, prob = 0.001),
      RObase == 4 ~ rbinom(nrow(.), size = 1, prob = 0.0001),
      TRUE ~ 0))
  
  
  ### * Alt: AgeDif Based  ================================================####
  kinships <- kinships %>%
    mutate(prob = case_when(RObase == 2 ~ 0.01 * exp(-0.06 * AgeDif),
      RObase == 4 ~ 0.0001 * exp(-0.06 * AgeDif),
      TRUE ~ 0),
    TruePairs = rbinom(nrow(.), size = 1, prob = prob))
  
  
    #### Check rates=======================================================####
    true_per_year <- kinships %>% 
      group_by(Cohort_2,AgeDif) %>% 
      dplyr::summarise(Pair_true_count = as.integer(sum(TruePairs > 0,
                                                        na.rm = TRUE))) %>% 
      filter(Pair_true_count > 0) %>%     # keep only cohorts with >0
      arrange(Cohort_2)
  
  
  
  
  
  
  
  
## >CKMRsim Pairing =======================================================####

  ### Load genetic data  ==================================================####
  vcf <- read.vcfR("[vcf]", verbose = FALSE)
  popmap <- read.table("[].pop.map.txt", header = FALSE, 
                       col.names = c("IND", "POP"))
  
  genind_obj <- vcfR2genind(vcf)
  pop(genind_obj) <- popmap$POP[match(indNames(genind_obj), popmap$IND)]
  
  # Keep as genind object for individual-based processing
  dat <- genind_obj
  
  
# EDIT THIS TO ADD BIODATA FROM ACTUAL SAMPLES
  # # Import Biodata from Access
  # db_path <- "C:/Users/eetry/Michigan State University/MEER Lab - Documents/Meer_Lab_Database/backup/meerLab_database_be_02022026.accdb"
  # 
  # con <- dbConnect(odbc(),
  #                  Driver = "Microsoft Access Driver (*.mdb, *.accdb)",
  #                  Dbq = db_path)
  # 
  # samples <- tbl(con, "samples") %>% as_tibble()
  # collections <- tbl(con, "collections") %>% as_tibble()
  # locations <- tbl(con, "locations") %>% as_tibble()
  # 
  # dbDisconnect(con)
  # 
  # # Set up biodata & map safely to individuals
  # db_metadata <- samples %>%  
  #   filter(sample_id %in% indNames(dat)) %>%  
  #   left_join(collections, join_by("collection_id")) %>%  
  #   left_join(locations, join_by("location_id")) %>%  
  #   select(sample_id, location_name) %>% 
  #   as_tibble()
  
  pop_ordering <- tibble(sample_id = indNames(dat)) %>% 
    left_join(db_metadata, by = "sample_id")
  
  pop(dat) <- factor(replace_na(pop_ordering$location_name, "Unknown"))
  
  ### Extract Allele Frequencies  =========================================####
  tab_mat <- tab(dat, freq = TRUE)
  allele_freqs <- colMeans(tab_mat, na.rm = TRUE)
  
  df <- data.frame(LocusAllele = names(allele_freqs),
    Freq = allele_freqs)
  
  df <- df %>% 
    separate(LocusAllele, into = c("Locus", "Allele"), sep = "\\.") %>% 
    mutate(Chrom = 1,
      Locus = as.character(Locus),
      Allele = as.character(Allele)) %>% 
    group_by(Locus) %>% 
    dplyr::mutate(Pos = dplyr::cur_group_id()) %>% 
    ungroup()
  
  long_markers <- df %>%
    select(Chrom, Locus, Pos, Allele, Freq) %>%
    as_tibble() %>%
    dplyr::group_by(Locus) %>%
    dplyr::mutate(Freq = Freq / sum(Freq)) %>%
    dplyr::ungroup() %>%
    arrange(Locus, Allele) %>%
    dplyr::group_by(Locus) %>%
    dplyr::mutate(LocIdx = dplyr::cur_group_id(),
      AlleIdx = row_number()) %>%
    dplyr::ungroup()
  
  ### Genind -> long_markers   ============================================####
  genind_to_long_markers <- function(genind_obj) {
    
    tab     <- genind_obj@tab
    loc_fac <- genind_obj@loc.fac
    loci    <- levels(loc_fac)
    
    out <- list()
    loc_idx <- 1L
    
    for (L in loci) {
      cols <- which(loc_fac == L)
      
      allele_names <- colnames(tab)[cols]
      allele_values <- sub(".*\\.", "", allele_names)
      
      ac <- colSums(tab[, cols, drop = FALSE], na.rm = TRUE)
      if (sum(ac) == 0) next
      
      freqs <- ac / sum(ac)
      
      # Critical filter: drop monomorphic loci
      if (sum(freqs > 0) <= 1) next
      
      # Convert character alleles (A,C,T,G) to safe positional integers (1, 2, 3...)
      allele_integers <- as.integer(factor(allele_values, levels = unique(allele_values)))
      
      out[[loc_idx]] <- data.frame(
        Chrom  = 1,
        Locus  = L,
        Pos    = loc_idx,
        Allele = allele_integers, 
        LocIdx = loc_idx,
        AlleIdx = seq_along(freqs),
        Freq   = as.numeric(freqs)
      )
      
      loc_idx <- loc_idx + 1L
    }
    
    dplyr::as_tibble(do.call(rbind, out))
  }
  
      ### Convert  ========================================================####
      dat.ckmrSim <- genind_to_long_markers(dat)
      
      # Extract genotype Matrix 
      geno <- tab(dat, NA.method = "asis")
      
      geno_df <- as.data.frame(geno) %>% 
        rownames_to_column(var = "Sample") %>% 
        as_tibble()
      
      
      
  
 
      
      
## >CKMRsim ===============================================================####
  ### Compute allele freq =================================================####
      geno_df <- as.tibble(geno_df)
      nc <- ncol(geno_df)
      
  ### Extract loci  =======================================================####
      long_genos <- geno_df %>% 
        tidyr::pivot_longer(-Sample, names_to = "loc", values_to = "count") %>% 
        tidyr::separate(loc, into = c("Locus", "Allele"), sep = "\\.") %>% 
        dplyr::mutate(
          count = as.integer(count),
          Allele = as.integer(Allele)
        ) %>% 
        dplyr::filter(!is.na(count) & count > 0) %>% 
        tidyr::uncount(count) %>% 
        dplyr::group_by(Sample, Locus) %>% 
        dplyr::mutate(gene_copy = row_number()) %>% 
        dplyr::ungroup() %>% 
        dplyr::rename(Indiv = Sample)


      # Notes from my lord and saviour, Eric Anderson:
      # It has four columns: Indiv, Locus, gene_copy, and Allele.
      # Alleles are named with characters (even if they are numbers, they must 
      # be coerced to characters).
      
      # The gene_copy column contains either a 1 or a 2 in every row, telling 
      # us which copy of the gene (in a diploid) is which allele.
      
      # Missing data in the Allele column is given by NA.
            
      
      long_genos_clean <- long_genos %>%
        filter(!is.na(Allele)) %>%
        mutate(gene_copy = as.integer(gene_copy),
          Allele = as.integer(Allele))
      
      long_genos_clean <- long_genos_clean %>%
        group_by(Indiv, Locus) %>%
        filter(n() == 2) %>%   # keep only proper diploid loci
        ungroup()
      
      long_genos_clean <- long_genos_clean %>%
        mutate(Allele = as.integer(gsub("^a", "", Allele)))
      
      loci <- locNames(dat)
      
      alle_freqs <- long_genos_clean %>%
        count(Locus, Allele) %>%
        group_by(Locus) %>%
        mutate(Freq = n / sum(n),
               Chrom = "Unk",
               Pos = as.integer(factor(Locus, levels = loci))) %>%
        ungroup() %>%
        select(Chrom, Pos, Locus, Allele, Freq) %>%
        arrange(Pos, desc(Freq)) %>%
        mutate(AlleIdx = NA,
               LocIdx = NA) %>%
        filter(!is.na(Allele)) # very important step!
      
      # see what that looks like:
      alle_freqs
      
      
      # Reindex 
      afreqs_ready <- reindex_markers(alle_freqs)
      
      
      
      # Check to make sure these allign
      setdiff(unique(long_genos_clean$Allele),
        unique(afreqs_ready$Allele))
      
      
      
      
  ### Create CKMR object  =================================================####
      
      kappas # this is the probability of relationship having 0,1,2 identical by descent genes
      
      ex1_ckmr <- create_ckmr(D = afreqs_ready, # tib of freq
        kappa_matrix = kappas[c("PO", "FS", "HS", "U"), ], # made, select the relationships you want
        ge_mod_assumed = ge_model_TGIE, # gt error model
        ge_mod_true = ge_model_TGIE, # gt error model true, lets you look at potential differences
        ge_mod_assumed_pars_list = list(epsilon = 0.01), 
        ge_mod_true_pars_list = list(epsilon = 0.01))
      
  
      # Examine    
      ex1_ckmr
      
      
      ex1_Qs <- simulate_Qij(ex1_ckmr, 
                             calc_relats = c("PO", "FS", "U","HS"),
                             sim_relats = c("PO", "FS", "HS", "U") )
      
      
      ex1_Qs
      
      
      #### Extract log lik ================================================####
      PO_U_logls <- extract_logls(ex1_Qs,
                                  numer = c(PO = 1),
                                  denom = c(U = 1)) # select only POP and unrelated
      ggplot(PO_U_logls,
             aes(x = logl_ratio, fill = true_relat)) +
        geom_density(alpha = 0.25)
      
      
      #### False positive and neg rates ===================================####
      
        ##### Calculate target FPR ========================================####
        # Based on number of samples (n *(n-1))
        target_fpr <- 0.1 * (190 * 189) ^ (-1)
         
        
        
        ##### Importance Sampling =========================================####
        
        ex1_PO_is <- mc_sample_simple(ex1_Qs, 
                                      nu = "PO",
                                      de = "U",
                                      lambda_stars = 5)
        
        
        ex1_HS_is <- mc_sample_simple(ex1_Qs, 
                                      nu = "HS",
                                      de = "U",
                                      lambda_stars = 5)
        
        
        
      
  ### Comparisons =========================================================####
      
    #### Duplicates  ======================================================####

      
    matchers <- find_close_matching_genotypes(LG = long_genos_clean,
                                                CK = ex1_ckmr,
                                                max_mismatch = 6)
    # If some matchers... take them out!
    long_genos.edited <- long_genos %>% 
      filter(long_genos$Indiv %in% matchers$indiv_2)
      
      
      
    #### Pair Likelihoods  ================================================####
     ##### POPs   =======================================================####
    
      
      po_pairwise_logls <- pairwise_kin_logl_ratios(
        D1 = long_genos.edited,
        D2 = long_genos.edited,
        CK = ex1_ckmr,
        numer = "PO",
        denom = "U",
        num_cores = 1
      )
      
    # The one true god, Eric Anderson
      # note, num_cores is set to 1 because more than that breaks
      # it on CRAN checks.  But, leave it blank (i.e., don't use the
      # num_core option) to, by default, use all
      # your cores in parallel on a non-windows machine.
      
    ###### Pair Inference  =================================================####
    
      inferred_POPs <- po_pairwise_logls %>%
        filter(logl_ratio > 55 ) %>% # choose this value based on the ratios from
                                      # the importance sampling step
        arrange(desc(logl_ratio))
      
      
      
  
      
    ##### HSPs   =======================================================####
    hs_pairwise_logls <- pairwise_kin_logl_ratios(D1 = long_genos_clean,
                                                    D2 = long_genos_clean,
                                                    CK = ex1_ckmr,
                                                    numer = "HS",
                                                    denom = "U", 
                                                    num_cores = 1)
    ###### Pair Inference  =================================================####
    
      inferred_HSPs <- hs_pairwise_logls %>%
        filter(logl_ratio > 4) %>%  # choose this value based on the ratios from
                                    # the importance sampling step
        arrange(desc(logl_ratio))
      
      

    ##### FSPs   =======================================================####
    
    # You want full sibs to make sure no FSPs were mistaken for POPs 
      fs_pairwise_logls <- pairwise_kin_logl_ratios(D1 = long_genos_clean,
                                                    D2 = long_genos_clean,
                                                    CK = ex1_ckmr,
                                                    numer = "FS",
                                                    denom = "U", 
                                                    num_cores = 1)
    ###### Pair Inference  =================================================####
    
      inferred_FSPs <- fs_pairwise_logls %>%
        filter(logl_ratio > 4) %>% # choose this value based on the ratios from
                                    # the importance sampling step
        arrange(desc(logl_ratio))
      
      
      
      
      
      