#!/usr/bin/env Rscript
# INSTALL THE FOLLOWING #
#pkgs=c("tidyverse", "igraph", "ggraph", "scriptName",
#       "readxl", "patchwork", "RColorBrewer",
#       "viridis", "parallel", "doParallel",
#       "furrr", "future"
#)
source("scripts/FUNCTIONS.R") # loads packages too

# Arguments
args = commandArgs(trailingOnly=T)
if (length(args) != 0) {
    stop (paste0"Usage:
            ", scriptName::current_cli_filename(),
          " --cores # --expr Expression_table_long_averaged_z.tsv",
          " --HVgenes high_var_genes.tsv"
          " --baits Baits.tsv --ETdir results/Main_analysis_edgeTable")
}

cores=args[(which(args == "--cores"))+1]

high_var_genes = args[(which(args == "--HVgenes"))+1] %>% 
    read_delim() %>% pull(gene_ID)

Exp_table_long_averaged_z_high_var = 
    args[(which(args == "--expr"))+1] %>% 
    read_delim() %>% 
    filter(gene_ID %in% high_var_genes)

Baits = args[(which(args == "--baits"))+1] %>% read_delim()
edgetable_dir = args[(which(args == "--ETdir"))+1]

# Wide z-score table for selected genes
z_score_wide <- Exp_table_long_averaged_z_high_var %>%
    select(gene_ID, tissue, z.score.TPM) %>%
    pivot_wider(names_from = tissue, values_from = z.score.TPM) %>%
    as.data.frame()

row.names(z_score_wide) <- z_score_wide$gene_ID
nreps = ncol(z_score_wide) - 1
# Gene correlation matrix
cor_matrix <- cor(t(z_score_wide[, -1]))

# Transform the lower triangle of the matrix into NA, since it is a duplication of the upper part
cor_matrix_upper_tri <- cor_matrix
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA

# Edge selection: select only statistically significant correlations
#' t-distribution approximation
#' For each correlation coefficient r, you approximate a t statistics.
#' The equation is t = r * ( (n-2) / (1 - r^2) )^0.5
plan( multisession, workers = cores )
print("")
print("Calculate edge table")
edge_table <- future_map(.progress = T, 1:nrow(cor_matrix_upper_tri), \(i) {
    
    r = cor_matrix_upper_tri[i,];
    as.data.frame(r) %>%
        mutate(from = row.names(cor_matrix),
               to = rownames(cor_matrix_upper_tri)[i]) %>%
        # remove the lower triangle
        filter(!is.na(r)) %>%
        # remove self-to-self correlations
        filter(from != to) %>%
        mutate(t = r*sqrt( (nreps-2) / (1-r^2) ) ) %>%
        mutate(p.value = case_when(
            t > 0 ~ pt(t, df = nreps-2, lower.tail = F),
            t <=0 ~ pt(t, df = nreps-2, lower.tail = T)
        )) %>%
        mutate(FDR = p.adjust(p.value, method = "fdr"),
               significant = ifelse(FDR < 0.01, T, F))
    
}) %>% list_rbind()

# edge_table <- cor_matrix_upper_tri %>%
# as.data.frame() %>%
# mutate(from = row.names(cor_matrix)) %>%
# pivot_longer(cols = !from, names_to = "to", values_to = "r") %>%
# filter(!is.na(r)) %>% # remove the lower triangle
# filter(from != to) %>% # remove self-to-self correlations
# mutate(t = r*sqrt( (nreps-2) / (1-r^2) ) ) %>%
# mutate(p.value = case_when(
# t > 0 ~ pt(t, df = nreps-2, lower.tail = F),
# t <=0 ~ pt(t, df = nreps-2, lower.tail = T)
# )) %>%
# mutate(FDR = p.adjust(p.value, method = "fdr"))

# Plot the distribution of r values


#### DECIDE r_cutoff ####
r_bins =
    edge_table %>%
    mutate(r_bins = round(r, digits = 5)) %>%
    filter(r_bins > 0) %>%
    group_by(r_bins, significant) %>%
    count

print("")
print("Finding lowest possible r cutoff")
plan( multisession, workers = cores )
bins_cumul =
    future_map(.progress=T, 1:nrow(r_bins), \(curRow) {
        r_value = r_bins$r_bins[curRow]
        data.frame(r_value = r_value,
                   n_edges = (r_bins %>% filter(r_bins >= r_value))$n %>% sum,
                   significant_edges = (r_bins %>% filter(r_bins >= r_value, significant))$n %>% sum
        ) %>%
            mutate(Prop_significant = significant_edges / n_edges)
    }) %>% list_rbind

r_cutoff = min((bins_cumul %>% filter(Prop_significant > .9))$r_value)

print(paste0("Best r threshold = ", r_cutoff))

edge_table %>%
    ggplot() +
    geom_histogram(aes(x = r,
                       fill = significant),
                   bins = 1000) +
    geom_vline(xintercept = r_cutoff,
               color = "black",
               linewidth = 0.5) +
    # check at which r value, the number of correlations starts to drop fast
    labs(fill = "FDR < 0.01",
         title = "Gene correlations") +
    scale_x_continuous(
        breaks = seq(-1, 1, .4)) +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.text = element_text(color = "black")
    )

ggsave("r_histogram.svg", height = 5.6, width = 8)

edge_table %>%
    filter(from %in% Baits$gene_ID |
               to %in% Baits$gene_ID) %>%
    ggplot() +
    geom_histogram(aes(x = r,
                       fill = significant),
                   bins = 100) +
    geom_vline(xintercept = r_cutoff,
               color = "black",
               linewidth = 0.85) +
    # check at which r value, the number of correlations starts to drop fast
    labs(fill = "FDR < 0.01",
         title = "Gene correlations - Bait genes") +
    scale_x_continuous(
        breaks = seq(-1, 1, .4)) +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.text = element_text(color = "black")
    )
ggsave("r_histogram_baitGenes.svg",
       height = 5.6, width = 8)

# co-expressed in the same direction
split_at <- 90000

if (dir.exists(edgetable_dir)) {
    for (i in list.files(path = edgetable_dir)) {
        paste0(edgetable_dir, "/", i) %>%
            file.remove()
    }
} else {
    dir.create(edgetable_dir)
}


filtered_edge_table <-
    edge_table %>%
    filter(r > r_cutoff) %>%
    relocate(from, to)

if (nrow(filtered_edge_table) < split_at) {
    fileName = paste0(edgetable_dir, "/edge_table_r",
                      r_cutoff, "_",
                      0, ".tsv")
    filtered_edge_table %>%
        mutate(group = 0) %>%
        write_delim(fileName,
                    quote = "none", append = F,
                    col_names = T, delim = "\t")
} else {
    
    tmp <-  filtered_edge_table %>%
        mutate(group = (row_number() - 1) %/% !! split_at) %>%
        group_split(group)
    for (x in seq_len(length(tmp))) {
        fileName = paste0(edgetable_dir, "/edge_table_r",
                          r_cutoff, "_",
                          x, ".tsv")
        write_delim(tmp[[x]], fileName,
                    quote = "none", append = F,
                    col_names = T, delim = "\t")
        
    }
}
