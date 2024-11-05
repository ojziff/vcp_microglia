# run this script in Rstudio server OnDemand Terminal Tab (not with conda env r4.0.3) with:
# cd /nemo/lab/patanir/home/shared/patani-collab/microglia-bulk-rnaseq
# sbatch -N 1 -c 12 --mem=0 -t 3:00:00 --wrap="Rscript /nemo/lab/patanir/home/users/ziffo/projects/microglia-bulk-rnaseq/scripts/microglia_deseq_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects
# sbatch -N 1 -c 12 --mem=72G -t 3:00:00 --wrap="Rscript /nemo/lab/patanir/home/users/ziffo/projects/microglia-bulk-rnaseq/scripts/microglia_deseq_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects
# sbatch -N 1 -c 6 --mem=50G -t 3:00:00 --wrap="Rscript /nemo/lab/patanir/home/users/ziffo/projects/microglia-bulk-rnaseq/scripts/microglia_deseq_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects 
# sbatch --cpus-per-task=10 --mem=1500G -N 1 --partition=hmem -t 3:00:00 --wrap="Rscript /nemo/lab/patanir/home/users/ziffo/projects/microglia-bulk-rnaseq/scripts/microglia_deseq_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects
library(here)
# nemo_path = here("/Volumes/lab-patanir/home/users/ziffo")
nemo_path = here("/nemo/lab/patanir/home/users/ziffo")
# shared_path = here("/Volumes/lab-patanir/home/users/ziffo/proj-luscombn-patani/working")
shared_path = here("/nemo/project/proj-luscombn-patani/working")
# collab_path = here("/Volumes/lab-patanir/home/shared/patani-collab")
collab_path = here("/nemo/lab/patanir/home/users/ziffo/patani-collab")
# proj_path = here("/Volumes/lab-patanir/home/users/ziffo/projects/microglia-bulk-rnaseq")
proj_path = here("/nemo/lab/patanir/home/users/ziffo/projects/microglia-bulk-rnaseq")
# load(here(proj_path,"scripts/microglia-bulk-rnaseq.RData")) # overwrites with old functions
# load(here(proj_path,"sample-details/metadata.RData"))
source(here(nemo_path,"scripts/functions/R_workspace.R"))

# # Metadata ----------------------------------------------------------------
metadata <- read_csv(here(collab_path, "microglia-bulk-rnaseq/sample-details/samplesheet.csv")) %>% clean_names() %>%
  mutate(database_dir = "/nemo/lab/patanir/home/shared/patani-collab/microglia-bulk-rnaseq", condition = paste(mutation, treatment, sep = "_"),
         mutation = factor(mutation, levels = c("ctrl", "vcp")), treatment = factor(treatment, levels = c("untx", "lps")),
         parent_cellline = case_when(cellline %in% c("ncrm1", "ncrme6") ~ "ncrm", cellline %in% c("vcpf10", "cb1d") ~ "cb1", TRUE ~ cellline), parent_cellline = as.factor(parent_cellline), cellline = as.factor(cellline),
         isogenic = case_when(cellline == "ncrme6" ~ "mutation inserted", cellline == "vcpf10" ~ "mutation corrected", TRUE ~ "no"), name = sample, sample = paste0(condition, "_R", replicate),
         celltype = "microglia", age = "ipsc", author = "Clarke", celltype_age=paste(celltype, age, sep="_"), celltype_author=paste(celltype, author, sep="_"), age_author=paste(age, author, sep="_"), celltype_age_author=paste(celltype,age,author,sep="_")) %>%
  select(sample, condition, mutation, treatment, replicate, cellline, name, genotype, parent_cellline, isogenic, database_dir, everything())
  
metadata.noiso <- metadata %>% filter(!cellline %in% c("vcpf10", "ncrme6"))
metadata.iso <- metadata %>% filter(parent_cellline %in% c("ncrm", "cb1"))
untx.metadata <- metadata %>% filter(treatment == "untx") %>% group_by(mutation) %>% mutate(denom=-1/(n()-1), corrected_pair = ifelse(parent_cellline=="cb1", 1, denom), inserted_pair =  ifelse(parent_cellline=="ncrm", 1, denom)) %>% ungroup() # account for isogenic pairing
untx.metadata.noiso <- untx.metadata %>% filter(!cellline %in% c("vcpf10", "ncrme6"))
untx.metadata.iso <- untx.metadata %>% filter(parent_cellline %in% c("ncrm", "cb1"))
lps.metadata <- metadata %>% filter(treatment == "lps") %>% group_by(mutation) %>% mutate(denom=-1/(n()-1), corrected_pair = ifelse(parent_cellline=="cb1", 1, denom), inserted_pair =  ifelse(parent_cellline=="ncrm", 1, denom)) %>% ungroup() # account for isogenic pairing
lps.metadata.noiso <- lps.metadata %>% filter(!cellline %in% c("vcpf10", "ncrme6")) 
lps.metadata.iso <- lps.metadata %>% filter(parent_cellline %in% c("ncrm", "cb1"))
ctrl.metadata <- metadata %>% filter(mutation == "ctrl")
ctrl.metadata.noiso <- ctrl.metadata %>% filter(!cellline %in% c("vcpf10", "ncrme6"))
ctrl.metadata.iso <- ctrl.metadata %>% filter(parent_cellline %in% c("ncrm", "cb1"))
vcp.metadata <- metadata %>% filter(mutation == "vcp")
vcp.metadata.noiso <- vcp.metadata %>% filter(!cellline %in% c("vcpf10", "ncrme6"))
vcp.metadata.iso <- vcp.metadata %>% filter(parent_cellline %in% c("ncrm", "cb1"))

# # Public datasets ----------------------------------------------------------------
# Comparative datasets: Zhang adult & fetal cell types. Abud Adult & fetal microglia. Reich ipsc microglia.
zhang <- read_csv(here(shared_path, "public-data/purified-cortex-human-zhang-2016/sample-details/samplesheet.csv")) %>%
  filter(!group %in% c("astrocyte_fetal_cortex", "astrocyte_mouse_cortex", "astrocyte_adult_tumour", "astrocyte_adult_hipponemous", "whole_adult_temporal_lobe")) %>%
  mutate(database_dir = "/nemo/project/proj-luscombn-patani/working/public-data/purified-cortex-human-zhang-2016", age = "adult", author = "Zhang", celltype_age=paste(celltype, age, sep="_"),
         celltype_author=paste(celltype, author, sep="_"), age_author=paste(age, author, sep="_"), celltype_age_author=paste(celltype,age,author,sep="_"), age = "biopsy", sample = paste0(group,"_R",replicate)) %>%
  select(group, replicate, age, celltype, author, celltype_age, celltype_author, age_author, celltype_age_author, sample, database_dir)
reich <- read_csv(here(shared_path, "public-data/ipsc-microglia-reich-2020/sample-details/samplesheet.csv")) %>%
  mutate(database_dir = "/nemo/project/proj-luscombn-patani/working/public-data/ipsc-microglia-reich-2020", name = differentiation, sample = paste0(group,"_R",replicate), 
         celltype = "microglia", age = "ipsc", author = "Reich", celltype_age=paste(celltype, age, sep="_"), celltype_author=paste(celltype, author, sep="_"), age_author=paste(age, author, sep="_"),
         celltype_age_author=paste(celltype,age,author,sep="_")) %>%
  select(group, replicate, age, celltype, author, celltype_age, celltype_author, age_author, celltype_age_author, sample, database_dir)
abud <- read_csv(here(shared_path, "public-data/ipsc-microglia-abud-2019/sample-details/samplesheet.csv")) %>%
  mutate(database_dir = "/nemo/project/proj-luscombn-patani/working/public-data/ipsc-microglia-abud-2019", sample = paste0(group,"_R",replicate), group = gsub("ipsc-microglia-abud-2019", "microglia_ipsc", group),
         age = case_when(group == "microglia_ipsc" ~ "ipsc", group == "adult-microglia" ~ "adult",  group == "fetal-microglia" ~ "fetal",), celltype = "microglia",  author = "Abud",
         celltype_age=paste(celltype, age, sep="_"), celltype_author=paste(celltype, author, sep="_"), age_author=paste(age, author, sep="_"), celltype_age_author=paste(celltype,age,author,sep="_")) %>%
  select(group, replicate, age, celltype, author, celltype_age, celltype_author, age_author, celltype_age_author, sample, database_dir)

clarke_zhang_reich_abud.metadata <- bind_rows(select(untx.metadata, group, replicate, age, celltype, author, celltype_age, celltype_author, age_author, celltype_age_author, sample, database_dir), 
                                              zhang, reich, abud)

# iPSC motor neurons
ipsc_motorneuron.metadata <- read_csv(here(collab_path,"motor-neuron-vcp-luisier-2018/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%  filter(day == 35) %>%
  mutate(database_dir = here(collab_path,"motor-neuron-vcp-luisier-2018"), name = gsub("d35_","",sample), author = "Luisier", celltype = "motorneuron", age = "ipsc", age_author=paste(age, author, sep="_"), celltype_age_author=paste(celltype,age,author,sep="_"))

# # iPSC astrocytes
# ipsc_astrocyte.metadata <- read_csv(here(collab_path,"astrocyte-fractionation-bulk-rnaseq/sample-details/sample_table.csv")) %>% distinct(fraction_cellline, .keep_all = TRUE) %>% filter(fraction == "whole", cellline != "ctrl2") %>% 
#   distinct(fraction_cellline, .keep_all = TRUE) %>%
#   mutate(database_dir = here(collab_path,"astrocyte-fractionation-bulk-rnaseq"), condition = vcp, genotype = tolower(genotype), group = vcp,
#          replicate = rep(1:2, times = 2),  sample = paste0("ac_", fraction,"_",condition, "_R",replicate), sample = gsub("whole_","who_",sample), celltype = "astrocyte", author = "Ziff", mutation = condition,
#          age = "ipsc", age_author=paste(age, author, sep="_"), celltype_age_author=paste(celltype,age,author,sep="_")) %>%
#   select(sample, group, database_dir, replicate, cellline, condition, genotype, celltype, age, author, age_author, celltype_age_author) 

# ipsc_mn_ac_mg.metadata = bind_rows(ipsc_motorneuron.metadata, ipsc_astrocyte.metadata, untx.metadata, reich, abud) %>% select(sample, celltype, condition, mutation, replicate, genotype, author, database_dir)

# clarke_zhang_reich_abud.ac.mn.metadata = bind_rows(clarke_zhang_reich_abud.metadata, ipsc_motorneuron.metadata, ipsc_astrocyte.metadata)


### Chiu et al mouse: 3 factors - LPS, disease stage, mutation
chiu.metadata <- read_csv(here(shared_path, "public-data/microglia-sod1-mouse-chiu-2013/sample-details/samplesheet.csv")) %>%
  mutate(database_dir = "/nemo/project/proj-luscombn-patani/working/public-data/microglia-sod1-mouse-chiu-2013", sample = paste0(group, "_R", replicate), treatment = factor(treatment, levels = c("untx", "lps")), 
         mutation = factor(mutation, levels = c("ctrl", "sod1")), background = tolower(gsub("/","_", `genetic background`))) %>%
  select(group, replicate, day, stage, mutation, treatment, sample, background, database_dir, run_accession)
end_stage.untx.chiu.metadata <- chiu.metadata %>% filter(stage == "end-stage", treatment == "untx")
# symptomatic.untx.chiu.metadata <- chiu.metadata %>% filter(stage == "symptomatic", treatment == "untx")
# pre_symptomatic.untx.chiu.metadata <- chiu.metadata %>% filter(treatment == "untx", stage == "pre-symptomatic")
end_stage.ctrl.chiu.metadata <- chiu.metadata %>% filter(stage == "end-stage", mutation == "ctrl")


# # DESeq2 ------------------------------------------------------------------
# clarke = DESeq.analysis(metadata = metadata, design = ~ 1)
# clarke_zhang_reich_abud = DESeq.analysis(metadata = clarke_zhang_reich_abud.metadata, design = ~ 1)
# ipsc_mn_ac_mg = DESeq.analysis(metadata = ipsc_mn_ac_mg.metadata, design = ~ 1)
# clarke_zhang_reich_abud_ac_mn = DESeq.analysis(metadata = clarke_zhang_reich_abud.ac.mn.metadata, design = ~ 1)

untx_vcp_vs_ctrl = DESeq.analysis(metadata = mutate(untx.metadata, "condition" = mutation), design = ~condition + corrected_pair + inserted_pair, contrast = "condition_vcp_vs_ctrl", plot.figures = TRUE, plot.path.prefix = "untx_vcp_vs_ctrl")#,
                                  # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ corrected_pair + inserted_pair + mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation") # Untreated samples only, VCP vs CTRL
untx_vcp_vs_ctrl.noiso = DESeq.analysis(metadata = mutate(untx.metadata.noiso, "condition" = mutation), design = ~ condition, contrast = "condition_vcp_vs_ctrl", plot.figures = TRUE, plot.path.prefix = "untx_vcp_vs_ctrl.noiso")#,
                                        # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation") # Untreated samples only, VCP vs CTRL - isogenics removed
untx_vcp_vs_ctrl.iso = DESeq.analysis(metadata = mutate(untx.metadata.iso, "condition" = mutation), design = ~ parent_cellline + condition, contrast = "condition_vcp_vs_ctrl", plot.figures = TRUE, plot.path.prefix = "untx_vcp_vs_ctrl.iso")#,
                                      # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ parent_cellline + mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation") # Untreated samples only, VCP vs CTRL - isogenics removed
lps_vcp_vs_ctrl = DESeq.analysis(metadata = mutate(lps.metadata, "condition" = mutation), design = ~condition + corrected_pair + inserted_pair, contrast = "condition_vcp_vs_ctrl", plot.figures = TRUE, plot.path.prefix = "lps_vcp_vs_ctrl")#,
                                 # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ corrected_pair + inserted_pair + mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation") # LPS treated samples only, VCP vs CTRL
lps_vcp_vs_ctrl.noiso = DESeq.analysis(metadata = mutate(lps.metadata.noiso, "condition" = mutation), design = ~ condition, contrast = "condition_vcp_vs_ctrl", plot.figures = TRUE, plot.path.prefix = "lps_vcp_vs_ctrl.noiso")#,
                                       # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation") # LPS treated samples only, VCP vs CTRL - isogenics removed
lps_vcp_vs_ctrl.iso = DESeq.analysis(metadata = mutate(lps.metadata.iso, "condition" = mutation), design = ~ parent_cellline + condition, contrast = "condition_vcp_vs_ctrl", plot.figures = TRUE, plot.path.prefix = "lps_vcp_vs_ctrl.iso")#,
                                     # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ parent_cellline + mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation") # LPS treated samples only, VCP vs CTRL - isogenics removed
ctrl_lps_vs_untx = DESeq.analysis(metadata = mutate(ctrl.metadata, "condition" = treatment), design = ~ cellline + condition, contrast = "condition_lps_vs_untx", plot.figures = TRUE, plot.path.prefix = "ctrl_lps_vs_untx")#,
                                  # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ cellline + treatment + treatment:IRFinder, irfinder_contrast_variable = "treatment") # CTRL samples only, LPS vs untx
ctrl_lps_vs_untx.noiso = DESeq.analysis(metadata = mutate(ctrl.metadata.noiso, "condition" = treatment), design = ~ cellline + condition, contrast = "condition_lps_vs_untx", plot.figures = TRUE, plot.path.prefix = "ctrl_lps_vs_untx.noiso")#,
                                        # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ cellline + treatment + treatment:IRFinder, irfinder_contrast_variable = "treatment") # CTRL samples only, LPS vs untxL - isogenics removed
ctrl_lps_vs_untx.iso = DESeq.analysis(metadata = mutate(ctrl.metadata.iso, "condition" = treatment), design = ~ cellline + condition, contrast = "condition_lps_vs_untx", plot.figures = TRUE, plot.path.prefix = "ctrl_lps_vs_untx.iso")#,
                                      # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ cellline + treatment + treatment:IRFinder, irfinder_contrast_variable = "treatment") # CTRL samples only, LPS vs untxL - isogenics removed
vcp_lps_vs_untx = DESeq.analysis(metadata = mutate(vcp.metadata, "condition" = treatment), design = ~ cellline + condition, contrast = "condition_lps_vs_untx", plot.figures = TRUE, plot.path.prefix = "vcp_lps_vs_untx")#,
                                 # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ cellline + treatment + treatment:IRFinder, irfinder_contrast_variable = "treatment") # VCP samples only, LPS vs untx
vcp_lps_vs_untx.noiso = DESeq.analysis(metadata = mutate(vcp.metadata.noiso, "condition" = treatment), design = ~ cellline + condition, contrast = "condition_lps_vs_untx", plot.figures = TRUE, plot.path.prefix = "vcp_lps_vs_untx.noiso")#,,
                                       # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ cellline + treatment + treatment:IRFinder, irfinder_contrast_variable = "treatment") # VCP samples only, LPS vs untx# lps_vs_untx.vcp_vs_ctrl = DESeq.analysis(metadata = metadata, design = ~ treatment * mutation, contrast = "treatmentlps.mutationvcp")
vcp_lps_vs_untx.iso = DESeq.analysis(metadata = mutate(vcp.metadata.iso, "condition" = treatment), design = ~ cellline + condition, contrast = "condition_lps_vs_untx", plot.figures = TRUE, plot.path.prefix = "vcp_lps_vs_untx.iso")#,
                                     # run_irfinder = TRUE, irfinder_names = "name", irfinder_design = ~ cellline + treatment + treatment:IRFinder, irfinder_contrast_variable = "treatment") # VCP samples only, LPS vs untx# lps_vs_untx.vcp_vs_ctrl = DESeq.analysis(metadata = metadata, design = ~ treatment * mutation, contrast = "treatmentlps.mutationvcp")
# lps_vs_untx.vcp_vs_ctrl = DESeq.analysis(metadata = metadata, design = ~ corrected_pair + inserted_pair + treatment * mutation, contrast = "treatmentlps.mutationvcp") #  + treatment:parent_cellline + treatment:mutation
# ir.lps_vs_untx.vcp_vs_ctrl = IRFinder.analysis(metadata = metadata, sample.names = "sample", variable.name = c("mutation","treatment"), design = ~ corrected_pair + inserted_pair + treatment * mutation * IRFinder, file.var = "sample", contrast = "treatmentlps.mutationvcp.IRFinderIR", ge.res = lps_vs_untx.vcp_vs_ctrl$res, irfinder.dir = irfinder_dir)
# lps_vs_untx.vcp_vs_ctrl.noiso = DESeq.analysis(metadata = metadata.noiso, design = ~ treatment * mutation, contrast = "treatmentlps.mutationvcp")
# lps_vs_untx.vcp_vs_ctrl.iso = DESeq.analysis(metadata = metadata.iso, design = ~ parent_cellline + treatment * mutation, contrast = "treatmentlps.mutationvcp")
# 
# ipsc_astrocyte.vcp_vs_ctrl = DESeq.analysis(metadata = ipsc_astrocyte.metadata, design = ~condition, contrast = "condition_vcp_vs_ctrl") 
#                                   
# chiu_end_stage_untx_sod1_vs_ctrl = DESeq.analysis(metadata = end_stage.untx.chiu.metadata, design = ~ background + mutation, contrast = "mutation_sod1_vs_ctrl", species = "mouse",
#                                                   run_irfinder = TRUE, irfinder_names = "run_accession", irfinder_design = ~ background + mutation + mutation:IRFinder, irfinder_contrast_variable = "mutation")
# # chiu_symptomatic_untx_sod1_vs_ctrl = DESeq.analysis(metadata = symptomatic.untx.chiu.metadata, design = ~ mutation, contrast = "mutation_sod1_vs_ctrl", species = "mouse") 
# # chiu_presymptomatic_untx_sod1_vs_ctrl = DESeq.analysis(metadata = pre_symptomatic.untx.chiu.metadata, design = ~ mutation, contrast = "mutation_sod1_vs_ctrl", species = "mouse") 
# chiu_ctrl_lps_vs_untx =  DESeq.analysis(metadata = end_stage.ctrl.chiu.metadata, design = ~ treatment, contrast = "treatment_lps_vs_untx", species = "mouse",
#                                         run_irfinder = TRUE, irfinder_names = "run_accession", irfinder_design = ~ treatment + treatment:IRFinder, irfinder_contrast_variable = "treatment")

# # # Mass Spectrometry DEP ----------------------------------------------------------------
# print("running DEP")
# # # https://www.bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html
# # # NCRCM_2 == NCRMC2 line. The NCRMC2 replicates are an isogenic inserted R191Q clone and can be treated in the same way as NCRME6. So the groups we can compare here are: 1) CNTRL vs VCP (cntrl1, cntrl3, cntrl4, NCRM1 vs CB1D, GliA, GliB). 2) Ctrl vs LPS (ctrl3, ctrl4 vs ctrl3LPS, ctrl4LPS). 3) Parent line vs isogenic insertion (ncrm1 vs ncrmc2, ncrme6). 4) isogenic corrected vs mutant (vcpf10 vs cb1d). 
# # all untreated & ml240 samples for PCA
# proteinGroups <- read_tsv(here(proj_path, "mass-spectrometry/proteinGroups.txt")) %>% clean_names() %>% rename(protein_ids = protein_i_ds)
# proteinGroups$gene_names %>% duplicated() %>% any()
# proteinGroups %>% group_by(gene_names) %>% dplyr::summarize(frequency = dplyr::n()) %>%  arrange(desc(frequency)) %>% filter(frequency > 1) # for genes without an annotated gene name use the Uniprot ID
# proteinGroups_all <- make_unique(proteinGroups, "gene_names", "protein_ids")
# proteinGroups_all$name %>% duplicated() %>% any()
# 
# # Generate a SummarizedExperiment object using an experimental design
# LFQ_columns <- grep("lfq_", colnames(proteinGroups_all)) # get LFQ column numbers
# experimental_design = proteinGroups_all %>% select(starts_with("lfq_intensity_")) %>% colnames() %>% gsub("lfq_intensity_", "", .) %>% as_tibble_col(column_name = "label") %>%
#   mutate(condition = case_when(grepl("lps", label) ~ "lps", grepl("ctrl|vpf10|ncrcm1", label) ~ "ctrl", grepl("cb1d|gli|ncrcme6|ncrcm2", label) ~ "vcp"), 
#          cellline = str_split_fixed(label,"_",2)[,1], technical_replicate = str_split_fixed(label,"_",2)[,2]) %>%
#   group_split(condition) %>% purrr::map_df(~.x %>% group_by(label) %>% mutate(replicate = cur_group_id())) %>% ungroup
# # data_se <- make_se(proteinGroups_all, LFQ_columns, experimental_design) # create summarised experiment - this log2 transforms assay data
# 
# 
# experimental_design_vcp_vs_ctrl = proteinGroups_vcp_vs_ctrl %>% select(starts_with("lfq_intensity_")) %>% colnames() %>% gsub("lfq_intensity_", "", .) %>% as_tibble_col(column_name = "label") %>%
#   mutate(condition = case_when(grepl("lps", label) ~ "lps", grepl("ctrl|vpf10|ncrcm1", label) ~ "ctrl", grepl("cb1d|gli|ncrcme6|ncrcm2", label) ~ "vcp"), cellline = str_split_fixed(label,"_",2)[,1], technical_replicate = str_split_fixed(label,"_",2)[,2], 
#          parent_cellline = case_when(cellline %in% c("ncrcm1", "ncrcm2", "ncrcme6") ~ "ncrm", cellline %in% c("vpf10", "cb1d") ~ "cb1", cellline %in% c("glia", "glib") ~ "gli", TRUE ~ cellline),
#          denom=-1/(n()-1), corrected_pair = ifelse(parent_cellline=="cb1", 1, denom), inserted_pair =  ifelse(parent_cellline=="ncrm", 1, denom)) %>% ungroup() %>%
#   group_split(condition) %>% purrr::map_df(~.x %>% group_by(label) %>% mutate(replicate = cur_group_id())) %>% ungroup %>% group_by(cellline, technical_replicate) %>% mutate(replicate = cur_group_id()) %>% ungroup()
# 
# # plot_frequency(data_se) # 2841 proteins quantified in each replicate is  different - missing values need imputing
# # plot_coverage(data_se) # barplot of the protein identification overlap between samples
# # plot_numbers(data_se) #  barplot of the number of identified proteins per samples
# # Filter out proteins with missing values across samples
# # filter_missval(data_se, thr = 0) # Filter proteins identified in all replicates of at least one condition - too stringent - left with 1675 / 2841 proteins
# # filter_missval(data_se, thr = 1) # # Less stringent filtering - Filter proteins identified in 1/6 replicates of at least one condition 1892 / 2841
# # filter_missval(data_se, thr = 2) # 2059 / 2841
# # filter_missval(data_se, thr = 3) # 2233 / 2841
# # filter_missval(data_se, thr = 4) # 2390 / 2841
# # filter_missval(data_se, thr = 5) # 2547 / 2841
# # filter_missval(data_se, thr = 6) # 2841 / 2841 # no filtering
# 
# # VCP vs CTRL (cntrl1, cntrl3, cntrl4, NCRM1 vs CB1D, GliA, GliB)
# proteinGroups_vcp_vs_ctrl = proteinGroups_all %>% select(!contains(c("lps"))) # "ctrl1_","ctrl2_","ctrl3_","ncrcm1","cb1d","glia","glib" everything(),
# experimental_design_vcp_vs_ctrl = proteinGroups_vcp_vs_ctrl %>% select(starts_with("lfq_intensity_")) %>% colnames() %>% gsub("lfq_intensity_", "", .) %>% as_tibble_col(column_name = "label") %>%
#   mutate(condition = case_when(grepl("lps", label) ~ "lps", grepl("ctrl|vpf10|ncrcm1", label) ~ "ctrl", grepl("cb1d|gli|ncrcme6|ncrcm2", label) ~ "vcp"), cellline = str_split_fixed(label,"_",2)[,1], technical_replicate = str_split_fixed(label,"_",2)[,2], 
#          parent_cellline = case_when(cellline %in% c("ncrcm1", "ncrcm2", "ncrcme6") ~ "ncrm", cellline %in% c("vpf10", "cb1d") ~ "cb1", cellline %in% c("glia", "glib") ~ "gli", TRUE ~ cellline),
#          denom=-1/(n()-1), corrected_pair = ifelse(parent_cellline=="cb1", 1, denom), inserted_pair =  ifelse(parent_cellline=="ncrm", 1, denom)) %>% ungroup() %>%
#   group_split(condition) %>% purrr::map_df(~.x %>% group_by(label) %>% mutate(replicate = cur_group_id())) %>% ungroup
# dep.vcp_vs_ctrl = DEP.analysis(protein_groups = proteinGroups_vcp_vs_ctrl, experimental_design = experimental_design_vcp_vs_ctrl, threshold = 15, denomenator = "ctrl", numerator = "vcp", 
#                                design = ~0 + condition + corrected_pair + inserted_pair + technical_replicate, file_path = here(proj_path, "mass-spectrometry"), file_prefix = "untreated.vcp_vs_ctrl")
# 
# ## VCP vs CTRL sensitivity analysis without isogenics NCRME2, NCRME6, VPF10
# proteinGroups_vcp_vs_ctrl.noiso = proteinGroups_all %>% select(!contains(c("lps","ncrcm2","ncrcme6","vpf10"))) # "ctrl1_","ctrl2_","ctrl3_","ncrcm1","cb1d","glia","glib" everything(),
# experimental_design_vcp_vs_ctrl.noiso = proteinGroups_vcp_vs_ctrl.noiso %>% select(starts_with("lfq_intensity_")) %>% colnames() %>% gsub("lfq_intensity_", "", .) %>% as_tibble_col(column_name = "label") %>%
#   mutate(condition = case_when(grepl("lps", label) ~ "lps", grepl("ctrl|vpf10|ncrcm1", label) ~ "ctrl", grepl("cb1d|gli|ncrcme6|ncrcm2", label) ~ "vcp"), cellline = str_split_fixed(label,"_",2)[,1], technical_replicate = str_split_fixed(label,"_",2)[,2]) %>%
#   group_split(condition) %>% purrr::map_df(~.x %>% group_by(label) %>% mutate(replicate = cur_group_id())) %>% ungroup
# dep.vcp_vs_ctrl.noiso = DEP.analysis(protein_groups = proteinGroups_vcp_vs_ctrl.noiso, experimental_design = experimental_design_vcp_vs_ctrl.noiso, threshold = 9, denomenator = "ctrl", numerator = "vcp", 
#                                      design = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spectrometry"), file_prefix = "untreated.vcp_vs_ctrl.noiso") 
# 
# ## LPS vs CTRL (ctrl3, ctrl4 vs ctrl3LPS, ctrl4LPS)
# proteinGroups_lps_vs_untx = proteinGroups_all %>% select(!contains(c("cb1d","gli","ncrcm","vpf10","ctrl1"))) # "ctrl1_","ctrl2_","ctrl3_","ncrcm1","cb1d","glia","glib" everything(),
# experimental_design_lps_vs_untx = proteinGroups_lps_vs_untx %>% select(starts_with("lfq_intensity_")) %>% colnames() %>% gsub("lfq_intensity_", "", .) %>% as_tibble_col(column_name = "label") %>%
#   mutate(condition = case_when(grepl("lps", label) ~ "lps", grepl("ctrl[3-4]_[1-3]", label) ~ "ctrl"), cellline = case_when(grepl("ctrl3", label) ~ "ctrl3",grepl("ctrl4", label) ~ "ctrl4"), cellline = str_split_fixed(label,"_",2)[,1], technical_replicate = str_split_fixed(label,"_",2)[,2]) %>%
#   group_split(condition) %>% purrr::map_df(~.x %>% group_by(label) %>% mutate(replicate = cur_group_id())) %>% ungroup
# dep.lps_vs_untx = DEP.analysis(protein_groups = proteinGroups_lps_vs_untx, experimental_design = experimental_design_lps_vs_untx, threshold = 6, denomenator = "ctrl", numerator = "lps", 
#                                design = ~0 + condition + cellline + technical_replicate, file_path = here(proj_path, "mass-spectrometry"), file_prefix = "ctrl.lps_vs_untx") 
# 
# # ## Parent line vs isogenic insertion (ncrm1 vs ncrmc2, ncrme6)
# # proteinGroups_ncrmc2e6_vs_ncrm1 = proteinGroups_all %>% select(!contains(c("cb1d","gli","vpf10","ctrl"))) # keep ncrm1, ncrmc2, ncrme6
# # LFQ_columns_ncrmc2e6_vs_ncrm1 <- grep("lfq_", colnames(proteinGroups_ncrmc2e6_vs_ncrm1)) # get LFQ column numbers
# # experimental_design_ncrmc2e6_vs_ncrm1 = proteinGroups_ncrmc2e6_vs_ncrm1 %>% select(starts_with("lfq_intensity_")) %>% colnames() %>% gsub("lfq_intensity_", "", .) %>% as_tibble_col(column_name = "label") %>%
# #   mutate(condition = case_when(grepl("ncrcm1", label) ~ "ncrcm1", TRUE ~ "ncrmc2e6"), cellline = str_split_fixed(label,"_",2)[,1], technical_replicate = str_split_fixed(label,"_",2)[,2]) %>%
# #   group_split(condition) %>% purrr::map_df(~.x %>% group_by(label) %>% mutate(replicate = cur_group_id())) %>% ungroup
# # dep.ncrmc2e6_vs_ncrm1 = DEP.analysis(protein_groups = proteinGroups_ncrmc2e6_vs_ncrm1, experimental_design = experimental_design_ncrmc2e6_vs_ncrm1, threshold = 6, denomenator = "ncrcm1", numerator = "ncrmc2e6", 
# #                                      design = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spectrometry"), file_prefix = "ncrmc2e6_vs_ncrm1") 
# # 
# # ## isogenic corrected vs mutant (vpf10 vs cb1d)
# # proteinGroups_cb1d_vs_vpf10 = proteinGroups_all %>% select(!contains(c("lps","ncrcm","ctrl","gli"))) # keep vpf10, cb1d
# # LFQ_columns_cb1d_vs_vpf10 <- grep("lfq_", colnames(proteinGroups_cb1d_vs_vpf10)) # get LFQ column numbers
# # experimental_design_cb1d_vs_vpf10 = proteinGroups_cb1d_vs_vpf10 %>% select(starts_with("lfq_intensity_")) %>% colnames() %>% gsub("lfq_intensity_", "", .) %>% as_tibble_col(column_name = "label") %>%
# #   mutate(condition = case_when(grepl("cb1d", label) ~ "cb1d", grepl("vpf10", label) ~ "vpf10"), cellline = str_split_fixed(label,"_",2)[,1], technical_replicate = str_split_fixed(label,"_",2)[,2]) %>%
# #   group_split(condition) %>% purrr::map_df(~.x %>% group_by(label) %>% mutate(replicate = cur_group_id())) %>% ungroup
# # dep.cb1d_vs_vpf10 = DEP.analysis(protein_groups = proteinGroups_cb1d_vs_vpf10, experimental_design = experimental_design_cb1d_vs_vpf10, threshold = 3, denomenator = "vpf10", numerator = "vpf10",
# #                                  design = ~0 + condition + technical_replicate, file_path = here(proj_path, "mass-spectrometry"), file_prefix = "cb1d_vs_vpf10") 
# 
# print("saving DEP")
# save(dep.vcp_vs_ctrl, dep.vcp_vs_ctrl.noiso,  dep.lps_vs_untx, 
#   # dep.ncrmc2e6_vs_ncrm1, dep.cb1d_vs_vpf10,
#   file = here(proj_path, "scripts/microglia_mass_spec.RData"))
# 
# print("saving")
# save.image("/nemo/home/ziffo/home/projects/microglia-bulk-rnaseq/scripts/microglia-bulk-rnaseq.RData")
# print("complete")

