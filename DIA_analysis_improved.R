
# From https://www.bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html

library("BiocStyle")
library("DEP")
library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("SummarizedExperiment")
library("tibble")
library("grid")
library("ComplexHeatmap")
library("openxlsx")
library("viridis")
library("readr")
library("ggthemes")
library("RColorBrewer")
library("openxlsx")
library("pals")
library("ggpubr")

#set a random seed -> todo check if DEP imputation uses this random seed
#set.seed(42)

data <- read.table(file.choose(), fill = TRUE, sep = "\t", header = TRUE, quote = "") # Choose file "BeforeImputation_with_MeCP2.txt"
colnames(data) <- gsub("T..", "", colnames(data), fixed = TRUE)
colnames(data)


# To visualize all duplicate values in order of frequency in the table
data %>% group_by(Genes) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

# To rename duplicate values
data$Genes %>% duplicated() %>% any()
data_unique <- data_unique <- make_unique(data, "Genes", "Protein.Ids", delim = ";")
data$name %>% duplicated() %>% any()

# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.intensity.", colnames(data_unique)) # get data column numbers
data_se <- make_se_parse(data_unique, LFQ_columns)
data_se
str(data_se)

plot_frequency(data_se)

# Less stringent filtering:
# Filter for proteins that are identified in all replicates -1 of at least one condition
data_filt <- filter_missval(data_se, thr = 1)
str(data_filt)
plot_numbers(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)
rmmeanSdPlot(data_norm)
meanSdPlot(data_filt)



# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# Perform a mixed imputation

colnames(data_norm)
str(data_norm)
data_norm_df <- get_df_wide(data_norm)
colnames(data_norm_df)

impute_group <- function(data_norm_df, #get_df_wide() of summerizedExperiment object
                         start, end, #start and end column indexes of group to impute
                         valid = 3, #number of missing values per group to be MNAR
                         shift = 1.8, #shift and scale of MAR imputation
                         scale = 0.3,
                         seed = 42 #random seed for MAR imputation
                         ){
  set.seed(seed)
  data_norm_df_x <- data_norm_df[, c(1, start:end)]
  
  #MNAR imputation
  protein_MNAR_x <- data_norm_df_x[rowSums(is.na(data_norm_df_x)) >= valid] %>%
    pull(name) %>%
    unique()
  
  MNAR_x <- names(data_norm) %in% proteins_MNAR_x
  
  data_norm_x <- data_norm[, start:end]
  
  data_imp_x <- DEP::impute(
    data_norm_x, 
    fun = "mixed",
    randna = !MNAR_x, # we have to define MAR which is the opposite of MNAR
    mar = "none", # imputation function for MAR
    mnar = "zero") # imputation function for MNAR
  
  data_imp_x[data_imp_x == 0] <- 1
  
  rm(MNAR_x,data_norm_x, data_norm_df_x)
  return(data_imp_x)
}

data_imp_1 <-impute_group(data_norm_df,2,5)
data_imp_2 <-impute_group(data_norm_df,6,9)
data_imp_3 <-impute_group(data_norm_df,10,13)
data_imp_4 <-impute_group(data_norm_df,14,17)
data_imp_5 <-impute_group(data_norm_df,18,21)
data_imp_6 <-impute_group(data_norm_df,22,25)
data_imp_7 <-impute_group(data_norm_df,26,29)
data_imp_8 <-impute_group(data_norm_df,30,33)
data_imp_9 <-impute_group(data_norm_df,34,37)
data_imp_10 <-impute_group(data_norm_df,38,41)
data_imp_11 <-impute_group(data_norm_df,42,45)
data_imp_12 <-impute_group(data_norm_df,46,49)
data_imp_13 <-impute_group(data_norm_df,50,53)
data_imp_14 <-impute_group(data_norm_df,54,57)
data_imp_15 <-impute_group(data_norm_df,58,61)

data_imp <- cbind(data_norm)
assay(data_imp) <- cbind(assay(data_imp_1),
                         assay(data_imp_2), 
                         assay(data_imp_3),
                         assay(data_imp_4),
                         assay(data_imp_5),
                         assay(data_imp_6),
                         assay(data_imp_7),
                         assay(data_imp_8), 
                         assay(data_imp_9), 
                         assay(data_imp_10),  
                         assay(data_imp_11), 
                         assay(data_imp_12), 
                         assay(data_imp_13), 
                         assay(data_imp_14), 
                         assay(data_imp_15))


data_imp_0 <- assay(data_imp)
data_imp_0[data_imp_0 == 0] <- 1 
assay(data_imp) <- data_imp_0

# We can impute the whole dataset and then impute all NAs remaining values from earlier with values from this new whole-imputation data frame

# -- SKIP this part if a prior imputation was already performed (and go directly to the next) --

data_imp_norm <- DEP::impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
data_imp_norm <- assay(data_imp_norm)

idx <- is.na(data_imp_0)
data_imp_0[idx] <- data_imp_norm[idx]

assay(data_imp) <- data_imp_0

plot_imputation(data_norm, data_imp)
plot_missval(data_imp) # This should result in "Error: No missing values in 'data_imp'"

assay(data_norm["MAGI1",]) # check one protein before...
assay(data_imp["MAGI1",]) # ... and after imputation

rm(data, data_filt, data_imp_0, data_norm, data_norm_df, data_imp_norm, data_se, data_unique, idx, LFQ_columns)


# -- SKIP the coming part if the script is run for the first time --

# With every round, imputation always changes! To avoid this, re-use previously imputed and exported data
# to avoid this define a random seed...


Original <- read.xlsx(file.choose()) # Choose file "Proteomics_iNGNs_DIA_R_DataOutput_impute_1+Gaussian_with-MeCP2.xlsx"
row.names(Original) <- Original[,2]
Original1 <- Original[, -c(1:2, 63:237)] # retain only gene names (row names) and replicate values
Original1 <- data.matrix(Original1)
assay(data_imp) <- Original1

is.matrix(Original1)
row.names(Original)
colnames(Original)


# Define the comparisons

assay(data_imp)

data_diff_all <- test_diff(data_imp, type = "manual", test=c("wt_d4._vs_wt_d0.", "wt_d4._vs_wt_d1.", "wt_d4._vs_wt_d2.", "wt_d4._vs_wt_d3.", 
                                                             "wt_d1._vs_wt_d0.", "wt_d2._vs_wt_d0.", "wt_d3._vs_wt_d0.", "wt_d2._vs_wt_d1.", "wt_d3._vs_wt_d2.",
                                                             "TET3KO_d4._vs_TET3KO_d0.", "TET3KO_d4._vs_TET3KO_d1.", "TET3KO_d4._vs_TET3KO_d2.", "TET3KO_d4._vs_TET3KO_d3.",
                                                             "TET3KO_d1._vs_TET3KO_d0.", "TET3KO_d2._vs_TET3KO_d0.", "TET3KO_d3._vs_TET3KO_d0.", "TET3KO_d2._vs_TET3KO_d1.", "TET3KO_d3._vs_TET3KO_d2.",
                                                             "MeCP2KO_d0._vs_TET3KO_d0.", "MeCP2KO_d1._vs_TET3KO_d1.", "MeCP2KO_d2._vs_TET3KO_d2.", "MeCP2KO_d3._vs_TET3KO_d3.", "MeCP2KO_d4._vs_TET3KO_d4."))

data_diff_wvT <- test_diff(data_imp, type = "manual", test=c("wt_d0._vs_TET3KO_d0.", "wt_d1._vs_TET3KO_d1.", "wt_d2._vs_TET3KO_d3.", "wt_d3._vs_TET3KO_d3.", "wt_d4._vs_TET3KO_d4."))

# Denote significant proteins based on user defined cutoffs

dep <- add_rejections(data_diff_wvT, alpha = 0.05, lfc = log2(1.5))

# Plot the first and second principal components: PDF 10 x 10 inches

library("RColorBrewer")
colnames(dep)

# -- Depending on what to show, the following can be used for the upcoming plots --

colnames(assay(dep))
dep <- dep[,-c(41:60)] # With this, MECP2 data is left out (for visualization)
colnames(assay(dep))


plot_pca_Eli <- function (dep, x = 1, y = 2, indicate = c("condition"), label = FALSE, n = 500, point_size = 4, 
                          label_size = 3, plot = TRUE) 
{
  if (is.integer(x)) 
    x <- as.numeric(x)
  if (is.integer(y)) 
    y <- as.numeric(y)
  if (is.integer(n)) 
    n <- as.numeric(n)
  if (is.integer(point_size)) 
    point_size <- as.numeric(point_size)
  if (is.integer(label_size)) 
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.numeric(x), length(x) == 1, is.numeric(y), length(y) == 
                            1, is.numeric(n), length(n) == 1, is.character(indicate), 
                          is.logical(label), is.numeric(point_size), length(point_size) == 
                            1, is.numeric(label_size), length(label_size) == 
                            1, is.logical(plot), length(plot) == 1)
  
  if (x > ncol(dep) | y > ncol(dep)) {
    stop(paste0("'x' and/or 'y' arguments are not valid\n", 
                "Run plot_pca() with 'x' and 'y' <= ", ncol(dep), 
                "."), call. = FALSE)
  }
  if (n > nrow(dep)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_pca() with 'n' <= ", 
                nrow(dep), "."), call. = FALSE)
  }
  
  columns <- colnames(colData(dep))
  
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_pca() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(paste0("'", paste0(indicate, collapse = "' and/or '"), 
                  "' column(s) is/are not present in ", deparse(substitute(dep)), 
                  ".\nValid columns are: '", paste(columns, 
                                                   collapse = "', '"), "'."), call. = FALSE)
    }
  }
  
  var <- apply(assay(dep), 1, sd)
  df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)], 
  ]
  pca <- prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame() %>% rownames_to_column() %>% 
    left_join(., data.frame(colData(dep), Sample = dep$condition), by = c(rowname = "ID"))
  percent <- round(100 * pca$sdev^2/sum(pca$sdev^2), 1)
  for (feat in indicate) {
    pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }
  p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC",y)))) +
    ggtitle(paste0("PCA plot - Top ", n, " variable proteins"))+
    labs(x = paste0("PC", x, ": ", percent[x], "%"), y = paste0("PC", y, ": ", percent[y], "%")) 
  
  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]]), 
                        size = point_size) + labs(col = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]], 
                            shape = pca_df[[indicate[2]]]), size = point_size) + 
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]], 
                            shape = pca_df[[indicate[2]]]), size = point_size) + 
      facet_wrap(~pca_df[[indicate[3]]])
    labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (plot) {
    return(p)
  }
  else {
    df <- pca_df %>% select(rowname, paste0("PC", c(x, 
                                                    y)), match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}

plot_pca_Eli(dep, x = 1, y = 2, n = 1000, point_size = 4) +
  theme_grey() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 13, face="bold"), legend.position = "bottom", legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  stat_ellipse(geom = "polygon",
               aes(fill = Sample), 
               alpha = 0.25, show.legend = FALSE) 




# Plot the Pearson correlation matrix: PDF 10 x 10 inches
plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")




#How many clusters should I set up for the k-means?
# Get an idea with Elbow method

library("factoextra")
library("NbClust")

fviz_nbclust(assay(dep), stats::kmeans, method = "wss") +
  labs(subtitle = "Elbow method")




# Plot a heatmap of all significant proteins with data centered per protein

plot_heatmap_Eli <- function (dep, type = c("contrast", "centered"), 
                              #centered -> data centered on mean
                              #contrast -> 
                              kmeans = FALSE, 
                              k = 6, 
                              col_limit = 6, 
                              indicate = NULL, 
                              #bool? for 
                              clustering_distance = c("euclidean","maximum", "manhattan", "canberra", 
                                                      "binary", "minkowski", "pearson", "spearman", 
                                                      "kendall", "gower"), 
                              row_font_size = 6, 
                              col_font_size = 10, 
                              plot = TRUE, ...) 
{
  #put these variables in correct datatype 
  k <- as.numeric(k)
  col_limit <- as.numeric(col_limit)
  row_font_size <- as.numeric(row_font_size)
  col_font_size <- as.numeric(col_font_size)
  
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.character(type), is.logical(kmeans), 
                          is.numeric(k), length(k) == 1, 
                          is.numeric(col_limit), length(col_limit) ==  1, 
                          is.numeric(row_font_size), length(row_font_size) == 1,
                          is.numeric(col_font_size), length(col_font_size) == 1, 
                          is.logical(plot), length(plot) == 1)
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)
  
  
  #create dfs for row and column data, check validity
  row_data <- rowData(dep, use.names = FALSE)
  col_data <- colData(dep) %>% as.data.frame()
  
  if (any(!c("label", "condition", "replicate") %in% 
          colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '", 
                deparse(substitute(dep)), "'"), call. = FALSE)
  }
  if (length(grep("_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."), 
         call. = FALSE)
  }
  if (!"significant" %in% colnames(row_data)) {
    stop(paste0("'significant' column is not present in '", 
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required column."), 
         call. = FALSE)
  }
  
  
  #set ha1
  if (!is.null(indicate) & type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'", 
            call. = FALSE)
  }
  
  if (!is.null(indicate) & type == "centered") {
    ha1 <- get_annotation(dep, indicate)
  }
  else {
    ha1 <- NULL
  }
  
  #get only significant values
  filtered <- dep[row_data$significant, ]
  
  if (any(is.na(assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dep)), 
            "'. ", "Using clustering_distance = 'gower'", 
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE
  }
  else {
    obs_NA <- FALSE
  }
  
  if (type == "centered") {
    rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
    rowData(filtered)$sd <- sd(assay(filtered), na.rm = TRUE)
    df <- (assay(filtered) - rowData(filtered, use.names = FALSE)$mean)/rowData(filtered, use.names = FALSE)$sd
  }
  
  if (type == "contrast") {
    df <- rowData(filtered, use.names = FALSE) %>% data.frame() %>% 
      column_to_rownames(var = "name") %>% select(ends_with("_diff"))
    colnames(df) <- gsub("_diff", "", colnames(df)) %>% 
      gsub("_vs_", " vs ", .)
    df <- as.matrix(df)
  }
  
  if (kmeans & obs_NA) {
    warning("Cannot perform kmeans clustering with missing values", 
            call. = FALSE)
    kmeans <- FALSE
  }
  
  if (kmeans & !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
    if (type == "centered") {
      order <- data.frame(df) %>% cbind(., cluster = df_kmeans$cluster) %>% 
        mutate(row = apply(.[, seq_len(ncol(.) - 1)], 
                           1, function(x) max(x))) %>% group_by(cluster) %>% 
        summarize(index = sum(row)/n()) %>% arrange(desc(index)) %>% 
        pull(cluster) %>% match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
    if (type == "contrast") {
      order <- data.frame(df) %>% cbind(df, cluster = df_kmeans$cluster) %>% 
        gather(condition, diff, -cluster) %>% group_by(cluster) %>% 
        summarize(row = mean(diff)) %>% arrange(desc(row)) %>% 
        pull(cluster) %>% match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
  }
  
  if (ncol(df) == 1) {col_clust = FALSE}
  else {col_clust = TRUE}
  
  if (nrow(df) == 1) {row_clust = FALSE}
  else {row_clust = TRUE}
  
  if (clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      return(dist)
    }
  }
  
  legend <- ifelse(type == "contrast", "log2 Fold change", 
                   "Z-score")
  f1 = circlize::colorRamp2(seq(-col_limit, 
                                col_limit, (col_limit/5)), rev(RColorBrewer::brewer.pal(11, 
                                                                                        "RdBu")))
  interesting = which(row.names(df) %in% c("NEUROG1", "NEUROG2", "TUBB3", "SOX2", "MKI67", "PCNA", "DCX", "NEFL", "NEFM", "CALB2", "SYT1", "SYT2", "CADM1", "JAK1", "SMAD1", "NANOG", "SYT4",
                                           "PAX6", "TET3", "MECP2", "DNMT3A", "DNMT3B", "DNMT1", "UHRF1"))
  ann = rowAnnotation(foo = anno_mark(at = interesting, labels = row.names(df[interesting,])))
  
  ht1 = Heatmap(df, col = f1, split = if (kmeans) {       df_kmeans$cluster
  }
  else {
    NULL
  }, cluster_rows = col_clust,  
  column_names_side = "top", 
  clustering_distance_rows = clustering_distance,  
  heatmap_legend_param = list(color_bar = "continuous", 
                              legend_direction = "horizontal", legend_width = unit(4, 
                                                                                   "cm"), title_position = "topcenter"), 
  name = legend, row_names_gp = gpar(fontsize = row_font_size), 
  column_names_gp = gpar(fontsize = col_font_size), top_annotation = ha1,
  right_annotation = ann,
  
  
  ...)
  
  
  if (plot) {
    draw(ht1, heatmap_legend_side = "top")
  }
  else {
    colnames(df) <- gsub(" ", "_", colnames(df))
    df <- df[, unlist(column_order(ht1))]
    if (kmeans) {
      df <- cbind(df, k = df_kmeans$cluster)
    }
    return <- df[unlist(row_order(ht1)), ]
    data.frame(protein = row.names(return), return) %>% mutate(order = row_number())
  }
}

plot_heatmap_Lukas <- function (dep,type = "centered", 
                                kmeans = TRUE, 
                                k = 8, 
                                col_limit = 2, 
                                indicate = NULL,
                                show_row_names = FALSE, 
                                clustering_distance = "euclidean", 
                                clustering_method_rows = "ward.D2", 
                                row_dend_side = "left", 
                                row_dend_gp = gpar(lwd = 1), 
                                row_title_side = "left",
                                row_font_size = 15, 
                                col_font_size = 10,
                                plot = TRUE,
                                cluster_columns=FALSE) 
{
  #put these variables in correct datatype 
  k <- as.numeric(k)
  col_limit <- as.numeric(col_limit)
  row_font_size <- as.numeric(row_font_size)
  col_font_size <- as.numeric(col_font_size)
  
  clustering_distance <- match.arg(clustering_distance)
  
  ha1 <- NULL

  #get only significant values
  filtered <- dep[row_data$significant,]
  

  rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
  rowData(filtered)$sd <- sd(assay(filtered), na.rm = TRUE)
  df <- (assay(filtered) - rowData(filtered, use.names = FALSE)$mean)/rowData(filtered, use.names = FALSE)$sd
    
  set.seed(1)
  k = 8
  df_kmeans <- kmeans(df, k)
  
  order <- data.frame(df) %>% cbind(., cluster = df_kmeans$cluster) %>% 
            mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>% group_by(cluster) %>% 
            summarize(index = sum(row)/n()) %>% 
            arrange(desc(index)) %>% 
            pull(cluster) %>% 
            match(seq_len(k), .)
  
  df_kmeans$cluster <- order[df_kmeans$cluster]
  
  legend <- ifelse(type == "contrast", "log2 Fold change", "Z-score")
  
  f1 = circlize::colorRamp2(seq(-col_limit, 
                                col_limit, (col_limit/5)), rev(RColorBrewer::brewer.pal(11, 
                                                                                        "RdBu")))
  interesting = which(row.names(df) %in% c("NEUROG1", "NEUROG2", "TUBB3", 
                                           "SOX2", "MKI67", "PCNA", 
                                           "DCX", "NEFL", "NEFM", 
                                           "CALB2", "SYT1", "SYT2", 
                                           "CADM1", "JAK1", "SMAD1", 
                                           "NANOG", "SYT4", "PAX6", 
                                           "TET3", "MECP2", "OTX2", 
                                           "ID4", "SYN3", "PROM1", 
                                           "PROX1", "BCL2", "BAX", 
                                           "CUX1", "BCL11B", "SATB2", 
                                           "CAMK4", "DNMT3A", "DNMT3B"
                                           "DNMT1", "UHRF1"))
  
  ann = rowAnnotation(foo = anno_mark(at = interesting, labels = row.names(df[interesting,])))
  
  ht1 = Heatmap(df, 
                col = f1, 
                split = if (kmeans) {df_kmeans$cluster} else {NULL}, 
                cluster_rows = col_clust,  
                column_names_side = "top", 
                clustering_distance_rows = clustering_distance,  
                heatmap_legend_param = list(color_bar = "continuous", 
                                            legend_direction = "horizontal", 
                                            legend_width = unit(4, "cm"), 
                                            title_position = "topcenter"), 
                name = legend, 
                row_names_gp = gpar(fontsize = row_font_size), 
                column_names_gp = gpar(fontsize = col_font_size), top_annotation = ha1,
                right_annotation = ann, 
                ...)
  
  

  draw(ht1, heatmap_legend_side = "top")

}



# Internal function to get ComplexHeatmap::HeatmapAnnotation object
get_annotation <- function(dep, indicate) {
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(indicate))
  
  # Check indicate columns
  col_data <- colData(dep) %>%
    as.data.frame()
  columns <- colnames(col_data)
  if(all(!indicate %in% columns)) {
    stop("'",
         paste0(indicate, collapse = "' and/or '"),
         "' column(s) is/are not present in ",
         deparse(substitute(dep)),
         ".\nValid columns are: '",
         paste(columns, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  if(any(!indicate %in% columns)) {
    indicate <- indicate[indicate %in% columns]
    warning("Only used the following indicate column(s): '",
            paste0(indicate, collapse = "', '"),
            "'")
  }
  
  # Get annotation
  anno <- dplyr::select(col_data, indicate)
  
  # Annotation color
  names <- colnames(anno)
  anno_col <- vector(mode="list", length=length(names))
  names(anno_col) <- names
  for(i in names) {
    var = anno[[i]] %>% unique() %>% sort()
    if(length(var) == 1)
      cols <- c("black")
    if(length(var) == 2)
      cols <- c("orangered", "cornflowerblue")
    if(length(var) < 7 & length(var) > 2)
      cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
    if(length(var) > 7)
      cols <- RColorBrewer::brewer.pal(length(var), "Set3")
    names(cols) <- var
    anno_col[[i]] <-  cols
  }
  
  # HeatmapAnnotation object
  ComplexHeatmap::HeatmapAnnotation(df = anno,
                                    col = anno_col,
                                    show_annotation_name = TRUE)
}

test <- plot_heatmap_Eli(dep, type = "centered", 
                 kmeans = TRUE, 
                 k = 8, 
                 col_limit = 2, 
                 show_row_names = FALSE, 
                 clustering_distance = "euclidean", 
                 clustering_method_rows = "ward.D2", 
                 row_dend_side = "left", 
                 row_dend_gp = gpar(lwd = 1), 
                 row_title_side = "left",
                 row_font_size = 15, 
                 col_font_size = 10, 
                 cluster_columns=FALSE,
                 plot = TRUE) 
# !!! WATCH OUT: the size of the plot panel in RStudio affects the correct indication of proteins. Reshape it BEFORE running the last command and DO NOT resize when exporting!




# Plot a volcano plot for a contrast: PDF 10 x 10 inches 


get_df_wide(dep)




plot_volcano_Eli <- function (dep, contrast, label_size = 3, y_limit = c(0, 50), x_limit = c(-40, 40), abs_lab_x = 15, lab_y = 35, adjusted = FALSE, 
                              plot = TRUE) 
{
  if (is.integer(label_size)) 
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.character(contrast), length(contrast) == 1, is.numeric(label_size), 
                          length(label_size) == 1, is.logical(adjusted), length(adjusted) == 1, is.logical(plot), 
                          length(plot) == 1)
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."), 
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 
      1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."), 
         call. = FALSE)
  }
  if (length(grep("_significant", colnames(row_data))) < 
      1) {
    stop(paste0("'[contrast]_significant' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."), 
         call. = FALSE)
  }
  if (length(grep(paste(contrast, "_diff", sep = ""), 
                  colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>% data.frame() %>% select(ends_with("_diff")) %>% 
      colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '", 
                                paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n", 
         valid_cntrsts_msg, call. = FALSE)
  }
  diff <- grep(paste(contrast, "_diff", sep = ""), 
               colnames(row_data))
  if (adjusted) {
    p_values <- grep(paste(contrast, "_p.adj", sep = ""), 
                     colnames(row_data))
  }
  else {
    p_values <- grep(paste(contrast, "_p.val", sep = ""), 
                     colnames(row_data))
  }
  signif <- grep(paste(contrast, "_significant", sep = ""), 
                 colnames(row_data))
  df <- data.frame(x = row_data[, diff], y = -log10(row_data[, 
                                                             p_values]), significant = row_data[, signif], name = row_data$name) %>% 
    filter(!is.na(significant)) %>% arrange(significant)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  
  
  Interesting <- ifelse(df$significant & (df$name == "NEFM" | df$name == "NEFL" | df$name == "TUBB3" | df$name ==  "DCX" | df$name ==  "NEUROG1" | df$name == "TET3" | df$name == "MECP2"  
                                          | df$name ==  "NEUROG2" | df$name ==  "NEFH" | df$name ==  "POU5F1" | df$name ==  "SOX2" | df$name ==  "NES" | df$name ==  "VIM" | df$name == "DNMT3A"
                                          | df$name == "DNMT3B" | df$name == "MAP2" | df$name == "MECP2" | df$name == "CDK2" | df$name == "SYT1" | df$name == "SYT4" | df$name == "NPTX1"
                                          | df$name == "OTX2" | df$name == "PAX6" | df$name == "ID4" | df$name == "PROM1" | df$name == "MKI67"), TRUE, FALSE)
  
  addline_format <- function(x,...){
    gsub("_vs_","vs ", x, fixed = TRUE)
    gsub("_", " ", x, fixed = TRUE)
  }
  
  p <- ggplot(df, aes(x, y, label = ifelse(Interesting == TRUE | (significant & (abs(x) > abs_lab_x & y > lab_y)), name,""))) +  
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_point(size=3.5, alpha=0.4, aes(col = ifelse(significant & x > 1.5, "Up-regulated", ifelse(significant & x < -1.5, "Down-regulated", ifelse(significant & x < 1.5 & x > -1.5, "Significant", "Non-significant"))))) + 
    ggrepel::geom_text_repel(box.padding = 0.5, nudge_x = 0.5, min.segment.length = 0, point.padding = unit(0.3, "cm"), max.overlaps = Inf,
                             size = label_size) +
    
    labs(title = addline_format(contrast), x = expression(log[2] ~ "Fold change")) + 
    theme_classic() +
    scale_colour_manual(name="Significance", breaks=c("Up-regulated", "Down-regulated", "Significant", "Non-significant"),values = c("firebrick3", "royalblue3","grey40","lightgrey")) +
    guides(colour = guide_legend(override.aes = list(size=4)))+
    theme(
      plot.title = element_text(size = 20, hjust=0.5),
      legend.title = element_text(size = 18, face="bold"), legend.position = "bottom",
      legend.text = element_text(size = 15),
      axis.text=element_text(size=18),
      axis.title=element_text(size=18),
      axis.line = element_line(colour = "black"),
      plot.margin = unit(c(1,1,1,1), "cm"))+
    scale_x_continuous(expand = c(0, 0), limits=x_limit) + scale_y_continuous(expand = c(0, 0), limits=y_limit)+
    geom_text(data = data.frame(), aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1.2, -0.2), vjust = c(-1, -1), label = c(name1, name2), 
                                       fontface = "bold"), size = 5)
  
  
  if (adjusted) {
    p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
  }
  else {
    p <- p + labs(y = expression(-log[10] ~ "P-value"))
  }
  if (plot) {
    return(p)
  }
  else {
    df <- df %>% select(name, x, y, significant) %>% arrange(desc(x))
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", 
                                  "p_value_-log10")
    if (adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}


v1 <- plot_volcano_Eli(dep, contrast = "wt_d4._vs_wt_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,80), lab_y = 50)
v2 <- plot_volcano_Eli(dep, contrast = "wt_d4._vs_wt_d1.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,80), lab_y = 50)
v3 <- plot_volcano_Eli(dep, contrast = "wt_d4._vs_wt_d2.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,60), lab_y = 40)
v4 <- plot_volcano_Eli(dep, contrast = "wt_d4._vs_wt_d3.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,60), lab_y = 40)

v5 <- plot_volcano_Eli(dep, contrast = "wt_d1._vs_wt_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,80), lab_y = 50)
v6 <- plot_volcano_Eli(dep, contrast = "wt_d2._vs_wt_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,80), lab_y = 50)
v7 <- plot_volcano_Eli(dep, contrast = "wt_d3._vs_wt_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,80), lab_y = 50)
v8 <- plot_volcano_Eli(dep, contrast = "wt_d2._vs_wt_d1.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,80), lab_y = 40)
v9 <- plot_volcano_Eli(dep, contrast = "wt_d3._vs_wt_d2.", label_size = 4.5, adjusted = FALSE, y_limit = c(0,60), lab_y = 40)

panel1 <- ggarrange(v1, v2, v3, v4, v5, v6, v7, v8, v9, nrow = 3, ncol = 3, common.legend = TRUE, legend = "bottom")
annotate_figure(plot1, top = text_grob("Wt comparisons", face = "bold", size = 18)) # Exported 20 x 20 inches (PDF)

v10 <- plot_volcano_Eli(dep, contrast = "TET3KO_d4._vs_TET3KO_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 50)
v11 <- plot_volcano_Eli(dep, contrast = "TET3KO_d4._vs_TET3KO_d1.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 50)
v12 <- plot_volcano_Eli(dep, contrast = "TET3KO_d4._vs_TET3KO_d2.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 40)
v13 <- plot_volcano_Eli(dep, contrast = "TET3KO_d4._vs_TET3KO_d3.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 40)

v14 <- plot_volcano_Eli(dep, contrast = "TET3KO_d1._vs_TET3KO_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 50)
v15 <- plot_volcano_Eli(dep, contrast = "TET3KO_d2._vs_TET3KO_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 50)
v16 <- plot_volcano_Eli(dep, contrast = "TET3KO_d3._vs_TET3KO_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 50)
v17 <- plot_volcano_Eli(dep, contrast = "TET3KO_d2._vs_TET3KO_d1.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 40)
v18 <- plot_volcano_Eli(dep, contrast = "TET3KO_d3._vs_TET3KO_d2.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 40)

panel2 <- ggarrange(v10, v11, v12, v13, v14, v15, v16, v17, v18, nrow = 3, ncol = 3, common.legend = TRUE, legend = "bottom")
annotate_figure(plot2, top = text_grob("TET3KO comparisons", face = "bold", size = 18))

v19 <- plot_volcano_Eli(dep, contrast = "TET3KO_d0._vs_wt_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 50)
v20 <- plot_volcano_Eli(dep, contrast = "TET3KO_d1._vs_wt_d1.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 40)
v21 <- plot_volcano_Eli(dep, contrast = "TET3KO_d2._vs_wt_d2.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 40)
v22 <- plot_volcano_Eli(dep, contrast = "TET3KO_d3._vs_wt_d3.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 45)
v23 <- plot_volcano_Eli(dep, contrast = "TET3KO_d4._vs_wt_d4.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 40)

panel3 <- ggarrange(v19, v20, v21, v22, v23, nrow = 3, ncol = 3, common.legend = TRUE, legend = "bottom")
annotate_figure(panel3, top = text_grob("TET3KO vs. wt comparisons", face = "bold", size = 18))

v24 <- plot_volcano_Eli(dep, contrast = "MeCP2KO_d0._vs_TET3KO_d0.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 100), lab_y = 50)
v25 <- plot_volcano_Eli(dep, contrast = "MeCP2KO_d1._vs_TET3KO_d1.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 50)
v26 <- plot_volcano_Eli(dep, contrast = "MeCP2KO_d2._vs_TET3KO_d2.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 50)
v27 <- plot_volcano_Eli(dep, contrast = "MeCP2KO_d3._vs_TET3KO_d3.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 60), lab_y = 50)
v28 <- plot_volcano_Eli(dep, contrast = "MeCP2KO_d4._vs_TET3KO_d4.", label_size = 4.5, adjusted = FALSE, y_limit = c(0, 80), lab_y = 50)

panel3 <- ggarrange(v24, v25, v26, v27, v28, nrow = 3, ncol = 3, common.legend = TRUE, legend = "bottom")
annotate_figure(panel3, top = text_grob("MeCP2KO vs. TET3KO comparisons", face = "bold", size = 18))


# Plot a barplot for selected proteins with the data centered

dep <- add_rejections(data_diff_all, alpha = 0.05, lfc = log2(1.5)) # If whole dataset is needed, otherwise go with the following
dep <- dep[,-c(41:60)] # MECP2 is now left out (for visualization)


plot_single_Eli <- function (dep, proteins, type = c("contrast", "centered"), 
                             plot = TRUE) 
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.character(proteins), is.character(type), is.logical(plot), 
                          length(plot) == 1)
  type <- match.arg(type)
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("label", "condition", "replicate") %in% 
          colnames(colData(dep)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '", 
         deparse(substitute(dep)), "'\nRun make_se() or make_se_parse() to obtain the required columns", 
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 
      1) {
    stop("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", 
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns", 
         call. = FALSE)
  }
  if (!"name" %in% colnames(row_data)) {
    stop("'name' column not present in '", deparse(substitute(dep)), 
         "'\nRun make_se() or make_se_parse() to obtain the required columns", 
         call. = FALSE)
  }
  if (all(!proteins %in% row_data$name)) {
    if (length(proteins) == 1) {
      rows <- grep(substr(proteins, 1, nchar(proteins) - 
                            1), row_data$name)
      possibilities <- row_data$name[rows]
    }
    else {
      rows <- lapply(proteins, function(x) grep(substr(x, 
                                                       1, nchar(x) - 1), row_data$name))
      possibilities <- row_data$name[unlist(rows)]
    }
    if (length(possibilities) > 0) {
      possibilities_msg <- paste0("Do you mean: '", 
                                  paste0(possibilities, collapse = "', '"), 
                                  "'")
    }
    else {
      possibilities_msg <- NULL
    }
    stop("please run `plot_single()` with a valid protein names in the 'proteins' argument\n", 
         possibilities_msg, call. = FALSE)
  }
  if (any(!proteins %in% row_data$name)) {
    proteins <- proteins[proteins %in% row_data$name]
    warning("Only used the following protein(s): '", 
            paste0(proteins, collapse = "', '"), "'")
  }
  subset <- dep[proteins]
  if (type == "centered") {
    means <- rowMeans(assay(subset), na.rm = TRUE)
    df_reps <- data.frame(assay(subset) - means) %>% rownames_to_column() %>% 
      gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset)), 
                                              by = "ID")
    df_reps$replicate <- as.factor(df_reps$replicate)
    df <- df_reps %>% group_by(condition, rowname) %>% summarize(mean = mean(val, 
                                                                             na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>% 
      mutate(error = sd/sqrt(n), CI.L = mean - error, CI.R = mean + error) %>% as.data.frame()
    df$rowname <- parse_factor(df$rowname, levels = proteins)
    p <- ggplot(df, aes(condition, mean)) + geom_hline(yintercept = 0) + 
      theme_DEP2() +
      geom_col(fill = "lightgrey", width = 0.8) + 
      geom_point(data = df_reps, aes(condition, val, fill = replicate), shape = 21, colour = "Black",
                 size = 5, position = position_dodge(width = 0.3)) + 
      scale_fill_manual(name= "Replicate", values= viridis::viridis(4, option = "D")) +
      geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.2) + 
      labs(y = expression(log[2] ~ "Mean-centered intensity " ~ "(± SEM)"), col = "Replicates") + 
      theme(axis.text.x = element_text(size=15, angle=90),
            axis.text.y = element_text(size=15)) +
      theme(legend.title = element_text(size=15, face="bold"), legend.text = element_text(size=13)) +
      theme(strip.background = element_rect(color="black", fill="lightgrey", linewidth=1, linetype="solid"),
            strip.text.x = element_text(size = 15, color = "black", face = "bold")) +
      facet_wrap(~rowname,  scales = "free_y")
  }
  if (type == "contrast") {
    df <- rowData(subset, use.names = FALSE) %>% data.frame() %>% 
      select(name, ends_with("_diff"), ends_with("_CI.L"), 
             ends_with("_CI.R")) %>% gather(var, val, 
                                            -name) %>% mutate(contrast = gsub("_diff|_CI.L|_CI.R", 
                                                                              "", var), var = gsub(".*_", "", 
                                                                                                   var)) %>% spread(var, val)
    df$name <- parse_factor(df$name, levels = proteins)
    suffix <- get_suffix(df$contrast)
    if (length(suffix)) {
      df$contrast <- delete_suffix(df$contrast)
    }
    p <- ggplot(df, aes(contrast, diff)) + geom_hline(yintercept = 0) + 
      geom_col(colour = "black", fill = "lightgrey", width = 0.6) + 
      geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.2) + 
      labs(x = suffix, y = expression(log[2] ~ "Fold change" ~ 
                                        "(± SEM)")) + facet_wrap(~name) + theme_DEP2()
  }
  if (plot) {
    return(p)
  }
  else {
    if (type == "centered") {
      df <- df %>% select(rowname, condition, mean, CI.L, 
                          CI.R)
      colnames(df) <- c("protein", "condition", 
                        "log2_intensity", "CI.L", "CI.R")
    }
    if (type == "contrast") {
      df <- df %>% select(name, contrast, diff, CI.L, CI.R) %>% 
        mutate(contrast = paste0(contrast, suffix))
      colnames(df) <- c("protein", "contrast", 
                        "log2_fold_change", "CI.L", "CI.R")
    }
    return(df)
  }
}


plot1 <- plot_single_Eli(dep, proteins = "TET3", type = "centered")
plot2 <- plot_single_Eli(dep, proteins = "MECP2", type = "centered")
plot3 <- plot_single_Eli(dep, proteins = "DNMT3B", type = "centered")
plot4 <- plot_single_Eli(dep, proteins = "NANOG", type = "centered")
plot5 <- plot_single_Eli(dep, proteins = "POU5F1", type = "centered")

ggarrange(plot1 + rremove("xlab"), plot2 + rremove("ylab") + rremove("xlab"), plot3 + rremove("ylab") + rremove("xlab"), 
          plot4, plot5 + rremove("ylab"), 
          labels = NULL, nrow = 3, ncol = 3, common.legend = TRUE, legend = "right", align = "hv")  # Exported 15 x 13 inches (PDF)

plot1 <- plot_single_Eli(dep, proteins = "FGFR1", type = "centered")
plot2 <- plot_single_Eli(dep, proteins = "PROM1", type = "centered")
plot3 <- plot_single_Eli(dep, proteins = "PAX6", type = "centered")
plot4 <- plot_single_Eli(dep, proteins = "SOX2", type = "centered")
plot5 <- plot_single_Eli(dep, proteins = "OTX2", type = "centered")
plot6 <- plot_single_Eli(dep, proteins = "ID4", type = "centered")
plot7 <- plot_single_Eli(dep, proteins = "PAX3", type = "centered")
plot8 <- plot_single_Eli(dep, proteins = "DLL1", type = "centered")
plot9 <- plot_single_Eli(dep, proteins = "HES1", type = "centered")

panel.plot <- ggarrange(plot1 + rremove("xlab"), plot2 + rremove("ylab") + rremove("xlab"), plot3 + rremove("ylab") + rremove("xlab"), 
                        plot4 + rremove("xlab"), plot5 + rremove("ylab") + rremove("xlab"), plot6 + rremove("ylab") + rremove("xlab"), 
                        plot7, plot8 + rremove("ylab"), plot9 + rremove("ylab"), nrow = 3, ncol = 3, common.legend = TRUE, legend = "right", align = "hv")
annotate_figure(panel.plot, top = text_grob("Developmental", face = "bold", size = 18)) 


plot1 <- plot_single_Eli(dep, proteins = "NEUROG2", type = "centered")
plot2 <- plot_single_Eli(dep, proteins = "NEUROG1", type = "centered")
plot3 <- plot_single_Eli(dep, proteins = "TUBB3", type = "centered")
plot4 <- plot_single_Eli(dep, proteins = "MAP2", type = "centered")
plot5 <- plot_single_Eli(dep, proteins = "DCX", type = "centered")
plot6 <- plot_single_Eli(dep, proteins = "SCRT1", type = "centered")
plot7 <- plot_single_Eli(dep, proteins = "NNAT.1", type = "centered") # For NNAT, check A. Rao's PNAS 2016
plot8 <- plot_single_Eli(dep, proteins = "HDAC1", type = "centered") 
plot9 <- plot_single_Eli(dep, proteins = "HDAC2", type = "centered")

panel.plot <- ggarrange(plot1 + rremove("xlab"), plot2 + rremove("ylab") + rremove("xlab"), plot3 + rremove("ylab") + rremove("xlab"), 
                        plot4 + rremove("xlab"), plot5 + rremove("ylab") + rremove("xlab"), plot6 + rremove("ylab") + rremove("xlab"), 
                        plot7, plot8 + rremove("ylab"), plot9 + rremove("ylab"), nrow = 3, ncol = 3, common.legend = TRUE, legend = "right", align = "hv")
annotate_figure(panel.plot, top = text_grob("Neuronal and others", face = "bold", size = 18)) 

plot1 <- plot_single_Eli(dep, proteins = "PCNA", type = "centered")
plot2 <- plot_single_Eli(dep, proteins = "MKI67", type = "centered")
plot3 <- plot_single_Eli(dep, proteins = "CDK5", type = "centered") # For CDK5, check V.K. Rao's Nucleic Acid Res 2019
plot4 <- plot_single_Eli(dep, proteins = "CDK10", type = "centered")
plot5 <- plot_single_Eli(dep, proteins = "CCNB1", type = "centered")
plot6 <- plot_single_Eli(dep, proteins = "CCND1", type = "centered")
plot7 <- plot_single_Eli(dep, proteins = "CCND3", type = "centered")
plot8 <- plot_single_Eli(dep, proteins = "CCNK", type = "centered")
plot9 <- plot_single_Eli(dep, proteins = "GMNN", type = "centered")

panel.plot <- ggarrange(plot1 + rremove("xlab"), plot2 + rremove("ylab") + rremove("xlab"), plot3 + rremove("ylab") + rremove("xlab"), 
                        plot4 + rremove("xlab"), plot5 + rremove("ylab") + rremove("xlab"), plot6 + rremove("ylab") + rremove("xlab"), 
                        plot7, plot8 + rremove("ylab"), plot9 + rremove("ylab"), nrow = 3, ncol = 3, common.legend = TRUE, legend = "right", align = "hv")
annotate_figure(panel.plot, top = text_grob("Proliferation", face = "bold", size = 18)) 



# Synaptogenesis 


plot1 <- plot_single_Eli(dep, proteins = "NRGN", type = "centered")
plot2 <- plot_single_Eli(dep, proteins = "NRXN2", type = "centered")
plot3 <- plot_single_Eli(dep, proteins = "NSG1", type = "centered") 
plot4 <- plot_single_Eli(dep, proteins = "VAMP7", type = "centered")
plot5 <- plot_single_Eli(dep, proteins = "CADM1", type = "centered")
plot6 <- plot_single_Eli(dep, proteins = "SLC17A6", type = "centered")
plot7 <- plot_single_Eli(dep, proteins = "STX17", type = "centered")
plot8 <- plot_single_Eli(dep, proteins = "SYN3", type = "centered")
plot9 <- plot_single_Eli(dep, proteins = "SYT4", type = "centered")

panel.plot <- ggarrange(plot1 + rremove("xlab"), plot2 + rremove("ylab") + rremove("xlab"), plot3 + rremove("ylab") + rremove("xlab"), 
                        plot4 + rremove("xlab"), plot5 + rremove("ylab") + rremove("xlab"), plot6 + rremove("ylab") + rremove("xlab"), 
                        plot7, plot8 + rremove("ylab"), plot9 + rremove("ylab"), nrow = 3, ncol = 3, common.legend = TRUE, legend = "right", align = "hv")
annotate_figure(panel.plot, top = text_grob("Synapses", face = "bold", size = 18)) 






# Correlation plot of significant proteins varying between 1M and 6M


df.dep <- get_df_wide(dep)
colnames(df.dep)

genotype_d4_d0 <- df.dep[df.dep$TET3KO_d4._vs_TET3KO_d0._significant == TRUE | df.dep$wt_d4._vs_wt_d0._significant == TRUE, ] 
genotype_d4_d1 <- df.dep[df.dep$TET3KO_d4._vs_TET3KO_d1._significant == TRUE | df.dep$wt_d4._vs_wt_d1._significant == TRUE, ] 
genotype_d4_d2 <- df.dep[df.dep$TET3KO_d4._vs_TET3KO_d2._significant == TRUE | df.dep$wt_d4._vs_wt_d2._significant == TRUE, ] 
genotype_d4_d3 <- df.dep[df.dep$TET3KO_d4._vs_TET3KO_d3._significant == TRUE | df.dep$wt_d4._vs_wt_d3._significant == TRUE, ] 

genotype_d1_d0 <- df.dep[df.dep$TET3KO_d1._vs_TET3KO_d0._significant == TRUE | df.dep$wt_d1._vs_wt_d0._significant == TRUE, ] 
genotype_d2_d0 <- df.dep[df.dep$TET3KO_d2._vs_TET3KO_d0._significant == TRUE | df.dep$wt_d2._vs_wt_d0._significant == TRUE, ] 
genotype_d3_d0 <- df.dep[df.dep$TET3KO_d3._vs_TET3KO_d0._significant == TRUE | df.dep$wt_d3._vs_wt_d0._significant == TRUE, ] 

genotype_d3_d2 <- df.dep[df.dep$TET3KO_d3._vs_TET3KO_d2._significant == TRUE | df.dep$wt_d3._vs_wt_d2._significant == TRUE, ]
genotype_d2_d1 <- df.dep[df.dep$TET3KO_d2._vs_TET3KO_d1._significant == TRUE | df.dep$wt_d2._vs_wt_d1._significant == TRUE, ]




dat <- genotype_d1_d0
A <- dat$wt_d3._vs_wt_d2._diff
B <- dat$TET3KO_d3._vs_TET3KO_d2._diff



plot <- ggplot(dat, aes(x = A, y = B, label = ifelse(((A > 20 & B > 1.5) | (A > 18 & B > 18) | 
                                                        (A < -19 & B < -1.5) |  (A < -16 & B < -16) |
                                                        (A > 1.5 & B < -1.5) | 
                                                        (A < -1.5 & B > 1.5)), name,"")))
plot + 
  theme_bw() +
  ggtitle(label="Correlation of significant proteins")+
  geom_vline(xintercept = 0, linewidth=1, linetype = "dotted") +
  geom_hline(yintercept = 0, linewidth=1, linetype = "dotted") +
  geom_vline(xintercept = 1.5, linewidth=1, linetype = "dotted") +
  geom_hline(yintercept = 1.5, linewidth=1, linetype = "dotted") +
  geom_vline(xintercept = -1.5, linewidth=1, linetype = "dotted") +
  geom_hline(yintercept = -1.5, linewidth=1, linetype = "dotted") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 16, face="bold"), legend.position = "bottom",
    legend.text = element_text(size = 14),
    axis.text=element_text(size=16),
    axis.title=element_text(size=16),
    axis.line = element_line(colour = "black")
  )+
  labs(x="wt, day3 vs. day2", y="TET3KO, day3 vs. day2") +
  geom_point(data=dat, aes(col = ifelse((A > 1.5 & B > 1.5) | (A < -1.5 & B < -1.5), "Convergent", 
                                        ifelse((A > 1.5 & B < -1.5) | (A < -1.5 & B > 1.5), "Divergent", "Other"))), size=3.5) +
  scale_colour_manual(name="Significance", breaks=c("Convergent", "Divergent", "Other"), values = c("seagreen1", "orchid", "grey30")) +
  ggrepel::geom_text_repel(box.padding = 0.3, nudge_x = 0.5, min.segment.length = 0, point.padding = unit(0.2, "cm"), max.overlaps = Inf,
                           size = 5) # Export 1200 x 1200

# ---

time_d3_d2 <- df.dep[df.dep$TET3KO_d3._vs_wt_d3._significant == TRUE | df.dep$TET3KO_d2._vs_wt_d2._significant == TRUE, ]
time_d4_d3 <- df.dep[df.dep$TET3KO_d4._vs_wt_d4._significant == TRUE | df.dep$TET3KO_d3._vs_wt_d3._significant == TRUE, ]
time_d2_d1 <- df.dep[df.dep$TET3KO_d2._vs_wt_d2._significant == TRUE | df.dep$TET3KO_d1._vs_wt_d1._significant == TRUE, ]
time_d1_d0 <- df.dep[df.dep$TET3KO_d1._vs_wt_d1._significant == TRUE | df.dep$TET3KO_d0._vs_wt_d0._significant == TRUE, ]

time_d3_d2_MeCP2 <- df.dep[df.dep$MeCP2KO_d3._vs_TET3KO_d3._significant == TRUE | df.dep$MeCP2KO_d2._vs_TET3KO_d2._significant == TRUE, ]
time_d4_d3_MeCP2 <- df.dep[df.dep$MeCP2KO_d4._vs_TET3KO_d4._significant == TRUE | df.dep$MeCP2KO_d3._vs_TET3KO_d3._significant == TRUE, ]
time_d2_d1_MeCP2 <- df.dep[df.dep$MeCP2KO_d2._vs_TET3KO_d2._significant == TRUE | df.dep$MeCP2KO_d1._vs_TET3KO_d1._significant == TRUE, ]
time_d1_d0_MeCP2 <- df.dep[df.dep$MeCP2KO_d1._vs_TET3KO_d1._significant == TRUE | df.dep$MeCP2KO_d0._vs_TET3KO_d0._significant == TRUE, ]



dat <- time_d1_d0_MeCP2
A <- dat$MeCP2KO_d0._vs_TET3KO_d0._diff
B <- dat$MeCP2KO_d1._vs_TET3KO_d1._diff

plot2 <- ggplot(dat, aes(x = A, y = B, label = ifelse(((A > 18 & B > 18) | 
                                                         (A < -19 & B < -1.5) |  (A < -17 & B < -17) |
                                                         (A > 1.5 & B < -1.5) | 
                                                         (A < -1.5 & B > 1.5)), name,"")))
plot2 + 
  theme_bw() +
  ggtitle(label="Correlation of significant proteins")+
  geom_vline(xintercept = 0, linewidth=1, linetype = "dotted") +
  geom_hline(yintercept = 0, linewidth=1, linetype = "dotted") +
  geom_vline(xintercept = 1.5, linewidth=1, linetype = "dotted") +
  geom_hline(yintercept = 1.5, linewidth=1, linetype = "dotted") +
  geom_vline(xintercept = -1.5, linewidth=1, linetype = "dotted") +
  geom_hline(yintercept = -1.5, linewidth=1, linetype = "dotted") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 16, face="bold"), legend.position = "bottom",
    legend.text = element_text(size = 14),
    axis.text=element_text(size=16),
    axis.title=element_text(size=16),
    axis.line = element_line(colour = "black")
  )+
  labs(x="MeCP2KO vs. TET3KO, day0", y="MeCP2KO vs. TET3KO, day1") +
  geom_point(data=dat, aes(col = ifelse((A > 1.5 & B > 1.5) | (A < -1.5 & B < -1.5), "Convergent", 
                                        ifelse((A > 1.5 & B < -1.5) | (A < -1.5 & B > 1.5), "Divergent", "Other"))), size=3.5) +
  scale_colour_manual(name="Significance", breaks=c("Convergent", "Divergent", "Other"), values = c("seagreen1", "orchid", "grey30")) +
  ggrepel::geom_text_repel(box.padding = 0.3, nudge_x = 0.5, min.segment.length = 0, point.padding = unit(0.2, "cm"), max.overlaps = Inf,
                           size = 5) # Export 1200 x 1200

#create a Summarized Experiment object for cluster in cmat
extract_cluster_se <-function(cmat, #matrix with cluster annotation
                               cluster#cluster number
                               )
  {
  cluster<- as.numeric(cluster)
  
  cmat_x <- cmat[cmat[,"clusterannot"] == cluster,,drop=FALSE]
  df_x <- as.data.frame(cmat_x)
  data_columns <- grep("d", colnames(df_x)) # get data column numbers
  se <- make_se_parse(df_x, data_columns)
  rm(cmat_x, df_x, data_columns)
  return(get_df_wide(se))
}

# IF DONE FOR THE FIRST TIME ONLY!!! Export following data.frames
# dep/xlsx is the one excel file i always import

xlsx <- get_df_wide(dep)
xlsx1 <- get_df_wide(data_norm)
write.xlsx(xlsx, "C:\\Users\\Elisa\\Documents\\Post-doc\\Tet3 Project\\Tet3-MeCP2 project\\iNGNs project\\Proteomics_DIA_Franzi23122022\\Proteomics_iNGNs_DIA_R_DataOutput_impute_1+Gaussian_with-MeCP2.xlsx"
           , sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE)
write.xlsx(xlsx1, "C:\\Users\\Elisa\\Documents\\Post-doc\\Tet3 Project\\Tet3-MeCP2 project\\iNGNs project\\Proteomics_DIA_Franzi23122022\\Proteomics_iNGNs_DIA_R_DataOutput_Pre-imputation_with-MeCP2.xlsx"
           , sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE)

# If done the first time, or if the tested contrasts are changed, do also the following. File generated from xlsx variable (see step above) should not be ever modified 

xlsx2 <- get_df_wide(dep)
xlsx2 <- xlsx2[, -grep("CI", names(xlsx2))]
write.xlsx(xlsx2, "C:\\Users\\Elisa\\Documents\\Post-doc\\Tet3 Project\\Tet3-MeCP2 project\\iNGNs project\\Proteomics_DIA_Franzi23122022\\Proteomics_iNGNs_DIA_R_DataOutput_analysis_with-MeCP2.xlsx",
           sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE)

