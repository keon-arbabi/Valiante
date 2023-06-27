# Bunchin.R
###########################################################################################################################- 

# load libraries ----
suppressPackageStartupMessages({
  library(scrattch.hicat)
  library(tidyverse)
  library(data.table)
  library(Seurat)
  library(glmGamPoi)
  library(matrixStats)
  library(ggsci)
  library(ggpubr)
  library(gridExtra)
})
options(stringsAsFactors=FALSE)

# functions ----
# convert gene names from MGI to HGNC
convertMouseGeneList = function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  mousex = unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(genesV2)
}

###########################################################################################################################- 
# prepare data ----
# load data 
matrix =  fread(file = "./data/Buchin/DG_RNAseq_data_processed.csv", data.table = FALSE) %>% column_to_rownames(var = "V1") 
# normalize 
matrix = logCPM(as.matrix(matrix))
anno = fread(file = "./data/Buchin/DG_RNAseq_data_processed_annotations.csv", data.table = FALSE) 
rownames(anno) = anno$title
# check 
all(rownames(anno) == colnames(matrix)) 

# read in predetermined genes from Allen Brain Map
high_dg   = intersect(scan("./data/Buchin/topDG_HBA.txt", what = "character", sep = "\n"), rownames(matrix))[1:100]
high_glut = intersect(scan("./data/Buchin/topGlut_CT.txt", what = "character", sep = "\n"), rownames(matrix))
high_gaba = intersect(scan("./data/Buchin/topGABA_CT.txt", what = "character", sep = "\n"), rownames(matrix))
high_glia = intersect(scan("./data/Buchin/topGlia_CT.txt", what = "character", sep = "\n"), rownames(matrix))
dex_genes = unique(c(high_dg, high_glut, high_gaba, high_glia)) 

# run K-means clustering on cells with enough genes detected
set.seed(1)  # For reproducibility
kpCells = anno$gene.counts.0 > 1000 
dat = t(matrix[dex_genes,kpCells])

rcl = kmeans(as.matrix(as.data.frame(dat)), 4)
rcl = kmeans(as.matrix(as.data.frame(dat)), 4)
rcl = kmeans(as.matrix(as.data.frame(dat)), 4)
# NOTE: this above line of code needs to be run at least 3 times to reproduce the results
# there is no reason why this should be the case, but it is.
cluster  = matrix[1,]*0
cluster[kpCells] = rcl$cluster

# assign clusters to appropriate classes based on marker gene expression
out_val   = NULL
for(i in 1:4) {
  out_val = cbind(out_val,Matrix::rowMeans(matrix[c("PROX1","GAD1","SLC17A7","SLC1A3"),cluster==i]>0))
}
val = 1
for (i in 1:4) {
  v   = which.max(out_val[i,])
  val = c(val, v+1)
  out_val[,v] = 0
}
classes  = c("Low quality", "Dentate Gyrus", "Inhibitory", "Excitatory", "Non-neuronal")[order(val)]
classes  = classes[cluster+1]    
anno$class_label = classes
print(table(classes[kpCells]))
# Dentate Gyrus    Excitatory    Inhibitory  Non-neuronal 
#        12011          6386          5302          3543 # these are the values you should get

# add Wyler grade 
anno = anno %>% 
  rename(donor_id = `Donor ID`, cell_id = title) %>%
  mutate(ep_state = case_when(donor_id %in% c("H16.06.008","H17.06.015") ~ "WG1",
                              donor_id %in% c("H16.06.009","H16.06.010") ~ "WG4")) %>%
  mutate(ep_state = factor(ep_state, levels = c("WG1","WG4"))) %>%
  mutate(class_label = factor(class_label, levels = c("Excitatory","Inhibitory","Non-neuronal","Dentate Gyrus","Low quality"), labels = c("Ex","Inh","Non-Neu","DG","Low quality"))) %>%
  relocate(c(donor_id, ep_state, class_label), .after = cell_id)
# save
#write_rds(anno, file = "./output/Buchin/anno.rds")

###########################################################################################################################- 
# seurat workflow ----

# load metadata 
anno = read_rds(file = "./output/Buchin/anno.rds") %>%
  # remove low quality cells
  filter(!class_label == "Low quality") 
# load expression data 
matrix = fread(file = "./data/Buchin/DG_RNAseq_data_processed.csv", data.table = FALSE) %>% column_to_rownames(var = "V1")
matrix[1:10,1:5]
# match and check 
matrix = matrix[, match(rownames(anno), colnames(matrix))]
all(rownames(anno) == colnames(matrix)) 

# create Seurat object
seu_obj = CreateSeuratObject(counts = matrix, meta.data = anno)

# # only consider cells with at least 200 detected genes and genes expressed in at least 5 cells
# selected_c = WhichCells(seu_obj, expression = nFeature_RNA > 200)
# selected_f = rownames(seu_obj)[rowSums(seu_obj) > 5]
# seu_obj = subset(seu_obj, features = selected_f, cells = selected_c)
# dim(seu_obj)

# basic qc check
seu_obj = PercentageFeatureSet(seu_obj, "^MT", col.name = "percent_mito")
seu_obj = PercentageFeatureSet(seu_obj, "^RP[SL]", col.name = "percent_ribo")
# violin plots 
feats = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "XIST")
VlnPlot(seu_obj, features = feats, group.by = "ep_state", pt.size = 0, ncol = length(feats)) + NoLegend()
#ggsave(filename = "./figures/Buchin/qc_violin_plots.png", width = 10, height = 5, dpi = 400)

# prepare for visualization 
seu_obj_vis = NormalizeData(seu_obj)
seu_obj_vis = FindVariableFeatures(seu_obj_vis, selection.method = "vst", nfeatures = 2000)
seu_obj_vis = ScaleData(seu_obj_vis)
seu_obj_vis = RunPCA(seu_obj_vis, features = VariableFeatures(object = seu_obj_vis))
#seu_obj_vis = subset(seu_obj_vis, subset = class_label == "DG")

Idents(seu_obj_vis) = "ep_state"
p1 = DimPlot(seu_obj_vis, reduction = "pca")
Idents(seu_obj_vis) = "class_label"
p2 = DimPlot(seu_obj_vis, reduction = "pca") 
# PCA plots 
p1 + p2 
#ggsave(filename = "./figures/Buchin/pca_plots.png", width = 12, height = 5, dpi = 400)

# normalize
# https://doi.org/10.1186/s13059-019-1874-1
seu_obj_norm = SCTransform(seu_obj, verbose = TRUE, return.only.var.genes = FALSE)
# save Seurat object
#write_rds(seu_obj_norm, file = "./output/Buchin/seu_obj_norm.rds")

###########################################################################################################################- 
# intraindividual gene-wise variance ----

# load metadata 
anno = read_rds(file = "./output/Buchin/anno.rds") %>% filter(!class_label == "Low quality") 
anno$class_label = droplevels(anno$class_label)
# load normalized Seurat object 
seu_obj_norm = readRDS(file = "./output/Buchin/seu_obj_norm.rds")
seu_obj_norm$class_label = anno$class_label
# get model attributes
sct_model_df = seu_obj_norm@assays$SCT@SCTModel.list$model1@feature.attributes

# visualize gene mean-variance relationship
p1 = sct_model_df %>%
  ggplot(., aes(x = log10(gmean), y = log10(variance))) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  geom_density_2d(size = 0.3) +
  labs(x = "log10 uncorrected mean UMI counts)",
       y = "log10 uncorrected variance") +
  ggtitle("uncorrected mean-variance relationship") +
  geom_smooth(method = "lm", se = F) +
  theme_classic()
p2 = sct_model_df %>%
  ggplot(., aes(x = log10(gmean), y = residual_variance)) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  geom_density_2d(size = 0.3) +
  labs(x = "corrected mean expression)",
       y = "corrected variance") +
  ggtitle("corrected mean-variance relationship") +
  geom_smooth(method = "lm", se = F) +
  theme_classic()
p1 + p2
#ggsave(filename = "./figures/Buchin/meanvar_plot.png", width = 12, height = 5, dpi = 400)

# initialize
res = NULL
# split cells by subclass and subject, then calculate gene-wise stats 
# (not very efficient, will fix later...)
for(CT in as.character(unique(anno$class_label))){
  for(INDV in unique(anno$donor_id)){
    # subset
    seu_obj_tmp = subset(seu_obj_norm, subset = class_label == CT & donor_id == INDV)
    # get normalized expression matrix (Pearson residuals)
    sct_matrix = seu_obj_tmp@assays$SCT@scale.data %>% as.matrix() 
    dim(sct_matrix)
    # get normalized counts matrix
    norm_counts_matrix = seu_obj_tmp@assays$SCT@counts %>% as.matrix()
    dim(norm_counts_matrix)
    # get original counts matrix
    counts_matrix = seu_obj_tmp@assays$RNA@counts %>% as.matrix()
    counts_matrix = counts_matrix[match(rownames(norm_counts_matrix),rownames(counts_matrix)), 
                                  match(colnames(norm_counts_matrix),colnames(counts_matrix))]
    dim(counts_matrix)
    # get gene-wise coefficient of variation (Stdev/Mean*100)
    tmp = data.frame(donor_id = INDV,
                     class_label = CT,
                     ep_state = unique(seu_obj_tmp$ep_state),
                     gene_symbol = rownames(sct_matrix),
                     n_cells = ncol(sct_matrix),
                     n_nonzero = rowSums(norm_counts_matrix != 0),
                     pct_nonzero = rowSums(norm_counts_matrix != 0)/ncol(norm_counts_matrix),
                     uncorrected_mean = rowMeans(counts_matrix),
                     uncorrected_var = rowVars(counts_matrix),
                     uncorrected_stdev = rowSds(counts_matrix),
                     count_mean = rowMeans(norm_counts_matrix),
                     count_var = rowVars(norm_counts_matrix),
                     count_stdev = rowSds(norm_counts_matrix),
                     resid_mean = rowMeans(sct_matrix),
                     resid_var = rowVars(sct_matrix),
                     resid_stdev = rowSds(sct_matrix)) %>%
      mutate(cv = count_stdev/(count_mean*100))
    # save
    res = rbind(res, tmp)
  }
}
res$ep_state = factor(res$ep_state, levels = c("WG1","WG4"))
# save results 
#write_rds(res, file = "./output/Buchin/res.rds")

###########################################################################################################################- 
# figures ----
# load 
res = read_rds(file = "./output/Buchin/res.rds") 

# distributions
p1 = ggplot(res, aes(x = resid_mean)) + # mean
  geom_histogram(bins = 100) + geom_vline(xintercept = median(res$resid_mean), color = "red") + scale_x_continuous(limits = c(-1.2, 1.2)) + theme_classic()
p2 = ggplot(res, aes(x = resid_var)) + # variance
  geom_histogram(bins = 100) + geom_vline(xintercept = median(res$resid_var), color = "red") + scale_x_continuous(limits = c(-1, 5)) + theme_classic()
p3 = ggplot(res, aes(x = resid_stdev)) + # standard deviation
  geom_histogram(bins = 100) + geom_vline(xintercept = median(res$resid_stdev), color = "red") + scale_x_continuous(limits = c(-1, 5)) + theme_classic()
p4 = ggplot(res, aes(x = cv)) + # coefficient of variation
  geom_histogram(bins = 150) + geom_vline(xintercept = median(res$count_cv), color = "red") + scale_x_continuous(limits = c(-0.2, 0.5)) + theme_classic()
p1 + p2 + p3 + p4
#ggsave(filename = "./figures/Buchin/stat_dists.png", width = 9, height = 7)

# get genes that are correlates of electrophysiology
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005814
ephys_genes = read.csv(file = "./data/Tripathy_2017_ephys_genes.csv")[,-1] %>% dplyr::rename(MGI.symbol = GeneSymbol)
# convert from mouse genes 
hgnc = convertMouseGeneList(ephys_genes$MGI.symbol)
ephys_genes = left_join(ephys_genes, hgnc, by = "MGI.symbol") %>% 
  filter(EphysProp %in% c("Rheo","Vrest")) %>% 
  pull(HGNC.symbol) %>% unique()

# get genes coding proteins with detected expression in human brain
# https://www.proteinatlas.org/search/NOT+tissue_category_rna%3Abrain%3BNot+detected
brain_genes = read_tsv(file = "./data/tissue_category_rna_brain_Tissue_detected.tsv", show_col_types = F) %>%
  filter(Gene %in% intersect(Gene, res$gene_symbol)) %>% pull(Gene) %>% unique()

# compare expression sd by cell type
p3 = res %>%
  filter(class_label == "DG") %>%
  filter(pct_nonzero > 0.1, 
         gene_symbol %in% ephys_genes) %>%
ggplot(., aes(x = ep_state, y = resid_stdev)) + 
  geom_violin(aes(fill = ep_state)) + 
  geom_boxplot(fill = "white", outlier.shape = NA, width = 0.4) +
  #ggbeeswarm::geom_quasirandom(na.rm = TRUE, shape = 16, alpha = 0.1, size = 0.6) + 
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(breaks = c("WG1","WG4"), values = c("darkgreen","darkorange")) +
  labs(x = "Wyler grade", y = "gene-wise standard deviation of expression") +
  #coord_flip() +
  facet_grid(~class_label, scales = "free") +
  theme_classic() +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.signif") 
ggsave(filename = paste0("./figures/Buchin/ephy_genes_sd_violin_ct.jpeg"), width = 5, height = 6.5)

# simulation ----

# initialize 
res_sim = NULL
plot_lst = NULL
# loop through each condition and cell type
for(ct in unique(res$class_label)){
  df = NULL
  for(cx in unique(res$ep_state)){
    sim = NULL
    tmp = res %>% filter(ep_state == cx, class_label == ct) 
    for(i in 1:10000){
      # simulate distribution of average standard deviations for brain expressed genes
      null_mean = tmp %>% filter(gene_symbol %in% sample(brain_genes, length(ephys_genes))) # random sampling 
      null_mean = mean(null_mean$resid_stdev)
      sim = c(sim, null_mean)
      print(paste(cx,ct,i))
    }
    # compare average stdev of ephys genes to simulated distribution 
    d_fun = ecdf(sim)
    h_mean = tmp %>% filter(gene_symbol %in% ephys_genes) %>% pull(resid_stdev) %>% mean
    tmp2 = data.frame(
      h_mean = h_mean,
      null_mean = mean(sim),
      null_sd = sd(sim),
      p_vals = round(d_fun(h_mean),4)
    )
    df = rbind(df, tmp2)
  }
  # save plot
  plot_lst[[paste(ct)]] = ggplot() +
    stat_function(fun = dnorm, args = list(mean = df$null_mean[1], sd = df$null_sd[1]), color = "darkorange") + 
    stat_function(fun = dnorm, args = list(mean = df$null_mean[2], sd = df$null_sd[2]), color = "darkgreen") + 
    geom_vline(xintercept = df$h_mean[1], color = "darkorange") +
    geom_vline(xintercept = df$h_mean[2], color = "darkgreen") +
    scale_x_continuous(limits = c(min(sim)-0.05, max(sim)+0.05)) +
    labs(x = "St. Dev", y = "Count", title = paste(ct,"\nWG4 p-value:",df$p_vals[1],"\nWG1 p-value:",df$p_vals[2])) +
    theme_classic()
}
png(filename = "./figures/Buchin/simulations_comb.png", width = 10, height = 4,  units = "in", res = 600)
grid.arrange(grobs = plot_lst, nrow = 1)
dev.off()











