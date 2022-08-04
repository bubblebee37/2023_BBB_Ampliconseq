#tximport of kallisto data
library(tximport)
library(limma)
library(statmod)
library(readr)
library(edgeR)
library(ggplot2)
folder_1 <-
  "/home/kyungha/Downloads/BBB_phage_display/Data/211216_RNAseq/c_to_m/"
file_list_c <- list.files(path = folder_1, pattern = "*.tsv")
for (i in 1:length(file_list_c)) {
  assign(file_list_c[i],
         read.csv(paste(folder_1, file_list_c[i], sep = '')))
}
folder_2 <-
  "/home/kyungha/Downloads/BBB_phage_display/Data/211216_RNAseq/t_to_m/"
file_list_t <- list.files(path = folder_2, pattern = "*.tsv")
for (i in 1:length(file_list_t)) {
  assign(file_list_t[i],
         read.csv(paste(folder_2, file_list_t[i], sep = '')))
}
setwd("~/Downloads/BBB_phage_display/Data/211216_RNAseq")
mouse_trans <-
  read.table(
    "/home/kyungha/Downloads/BBB_phage_display/Data/211216_RNAseq/mus_musculus_transcripts_to_genes_v2.csv",
    sep = ",",
    header = TRUE
  )
head(mouse_trans)
setwd("~/Downloads/BBB_phage_display/Data/211216_RNAseq/c_to_m/")
txi_c <-
  tximport(file_list_c, type = "kallisto", tx2gene = mouse_trans, ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")
names(txi_c)
c_counts <- txi_c$counts
colnames(c_counts) <- c("Invivo_1","Invivo_2", "Invivo_3", "Chip2um_1", "Chip2um_2", "Chip2um_3")
head(c_counts)
setwd("~/Downloads/BBB_phage_display/Data/211216_RNAseq/t_to_m/")
txi_t <-
  tximport(file_list_t, type = "kallisto", tx2gene = mouse_trans, ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")
names(txi_t)
#Create count file for each sample
t_counts <- txi_t$counts
colnames(t_counts) <- c("Invivo_1","Invivo_2", "Invivo_3", "Transwell_1", "Transwell_2", "Transwell_3")
setwd("~/Downloads/BBB_phage_display/Data/211216_RNAseq")
write.table(c_counts, "Chip2um_count_RNAseq_220117.csv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(t_counts, "Transwell_count_RNAseq_220117.csv", sep = "\t", quote = FALSE, row.names = TRUE)

#Find DEG & add each file with CPM
setwd("~/Downloads/BBB_phage_display/Data/211216_RNAseq/c_to_m/")
col_c_c_m <- data.frame(condition = factor(rep(c("Invivo", "chip2um"), each = 3)))
rownames(col_c_c_m) <- colnames(c_counts)
head(col_c_c_m)
group_c <- c('Invivo', 'Invivo', 'Invivo', 'chip2um', 'chip2um', 'chip2um')
y_c <- DGEList(txi_c$counts, group = group_c)
# filtering using filterByExpr : delete gene with minimum count less than 4 (total 6 samples)
keep_c <- filterByExpr(y_c, min.count = 4)
y_c <- y_c[keep_c, ]
y_c <- calcNormFactors(y_c)
design_c <- model.matrix(~0+condition, data = col_c_c_m)
v_c <- voom(y_c, design_c, plot = T)
fit_c <- lmFit(v_c, design_c)
contr_c <- makeContrasts(conditionchip2um - conditionInvivo, levels = colnames(coef(fit_c)))
fit_c <- contrasts.fit(fit_c, contr_c)
fit_c <- eBayes(fit_c)
setwd("~/Downloads/BBB_phage_display/Data/211216_RNAseq/")
list_invivo_chip <- topTable(fit_c, coef=1, confint=TRUE, n=Inf, adjust="BH")
total_list_c <- list_invivo_chip
f_c_out <- 'DEG_from_invivo_chip'
DE_c_keep <- total_list_c$adj.P.Val < 0.05 & abs(total_list_c$logFC) >= 1
DE_c_list <- total_list_c[DE_c_keep,]
f_c_all <- paste(f_c_out, 'limma_all.txt', sep = '.')
write.table(total_list_c, f_c_all, sep = '\t')
f_c_DE <- paste(f_c_out, 'limma_DE.txt', sep = '.')
write.table(DE_c_list, f_c_DE, sep = '\t')
f_c_cpm <- paste(f_c_out, 'abundance.txt', sep = '.')
write.table(txi_c$abundance, f_c_cpm, sep = '\t')

setwd("~/Downloads/BBB_phage_display/Data/211216_RNAseq/t_to_m/")
col_c_t_m <- data.frame(condition = factor(rep(c("invivo", "transwell"), each = 3)))
rownames(col_c_t_m) <- colnames(t_counts)
#col_c_t_m$condition <- factor(col_c_t_m$condition, levels = c("Transwell", "Invivo"))
head(col_c_t_m)
group_t <- c('invivo', 'invivo', 'invivo', 'transwell', 'transwell', 'transwell')
y_t <- DGEList(txi_t$counts, group = group_t)
# filtering
keep_t <- filterByExpr(y_t, min.count = 4)
y_t <- y_t[keep_t, ]
y_t <- calcNormFactors(y_t)
design_t <- model.matrix(~0+condition, data = col_c_t_m)
v_t <- voom(y_t, design_t, plot = T)
fit_t <- lmFit(v_t, design_t)
contr_t <- makeContrasts(conditiontranswell - conditioninvivo, levels = colnames(coef(fit_t)))
fit_t <- contrasts.fit(fit_t, contr_t)
fit_t <- eBayes(fit_t)
setwd("~/Downloads/BBB_phage_display/Data/211216_RNAseq/")
list_invivo_well <- topTable(fit_t, coef=1, confint=TRUE, n=Inf, adjust="BH")
total_list_t <- list_invivo_well
f_t_out <- 'DEG_from_invivo_well'
DE_t_keep <- total_list_t$adj.P.Val < 0.05 & abs(total_list_t$logFC) >= 1
DE_t_list <- total_list_t[DE_t_keep,]
f_t_all <- paste(f_t_out, 'limma_all.txt', sep = '.')
write.table(total_list_t, f_t_all, sep = '\t')
f_t_DE <- paste(f_t_out, 'limma_DE.txt', sep = '.')
write.table(DE_t_list, f_t_DE, sep = '\t')
f_t_cpm <- paste(f_t_out, 'abundance.txt', sep = '.')
write.table(txi_t$abundance, f_t_cpm, sep = '\t')
