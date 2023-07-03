# Madhavan Proteostasis R-BHB Manuscript R File
# Base R version 4.2.2 
# Version Annotation: 2022_10_31 "Innocent and Trusting"
# MS/MS Datasets are from Christina King (Birgit Schilling) at Buck Institute

# Key of Colors for Manuscript
# 1/0 mM Ex Vivo: #44AA99
# 5/0 mM Ex Vivo: #117733
# Fraction 1: #332288
# Fraction 2: #88CCEE
# Fraction 3: #CC6677
# Fraction 4: #DDCC77
# All Fractions: #332288
# Primary Targets: #6699CC

# Set Working Directory
setwd('/Users/sidmadhavan/Documents/Mass Spectrometry/')
options(stringsAsFactors = F)











##################################################
# EX VIVO DATA
##################################################

# Load data from Christina, I abbreviated, alphabetized by UNIPROT and corrected
# where labeling of gene symbols was incorrect bc of Excel
library(readxl)

primtarg1 <- read_excel("2022_11_26.ExVivo.0.1mM.xlsx")
primtarg2 <- as.matrix(primtarg1)
primtarg2 <- as.data.frame(primtarg2)

sectarg1 <- read_excel("2022_11_26.ExVivo.0.5mM.xlsx")
sectarg2 <- as.matrix(sectarg1)
sectarg2 <- as.data.frame(sectarg2)

mastertarg1 <- read_excel("2022_11_26.ExVivo.REPLICATE.xlsx")
mastertarg2 <- as.matrix(mastertarg1)
mastertarg2 <- as.data.frame(mastertarg2)


# Add replicates to target dfs, clean off unneeded, add log2 and QValue(numeric)
primtarg3 <- mastertarg2[mastertarg2$C_UNIPROT %in% primtarg2$C_UNIPROT,]
primtarg3 <- primtarg3[,-c(1, 11:14)]
primtarg4 <- cbind(primtarg3, as.numeric(primtarg2$`log2(1/0)`),
                   as.numeric(primtarg2$`Q Value`))

sectarg3 <- mastertarg2[mastertarg2$C_UNIPROT %in% sectarg2$C_UNIPROT,]
sectarg3 <- sectarg3[,-c(1, 7:10)]
sectarg4 <- cbind(sectarg3, as.numeric(sectarg2$`log2(1/0)`),
                  as.numeric(sectarg2$`Q Value`))


# Clean extraneous details from Christina's UNIPROT IDs, name columns for ease
primgenes1 <- as.matrix(primtarg4$C_UNIPROT)
primgenes1 <- gsub(";...*","", primgenes1)
primtarg5 <- cbind(primtarg4,primgenes1)
colnames(primtarg5)[10] <- "log2FC"
colnames(primtarg5)[11] <- "QValue"
colnames(primtarg5)[12] <- "UNIPROT"

secgenes1 <- as.matrix(sectarg4$C_UNIPROT)
secgenes1 <- gsub(";...*","", secgenes1)
sectarg5 <- cbind(sectarg4,secgenes1)
colnames(sectarg5)[10] <- "log2FC"
colnames(sectarg5)[11] <- "QValue"
colnames(sectarg5)[12] <- "UNIPROT"


# Create background lists with total detectable proteome in fractions
masback1 <- as.matrix(mastertarg2$C_UNIPROT)
masback2 <- gsub(";...*","", masback1)
colnames(masback2)[1] <- "UNIPROT"


# Export directional bitr input (Uniprot Accessions) to 
# https://www.uniprot.org/id-mapping to convert
Uniprim <- primtarg5$UNIPROT[primtarg5$log2FC > 0]
write.csv(Uniprim, file="Uniprim.csv", row.names = F)
Unisec <- sectarg5$UNIPROT[sectarg5$log2FC > 0]
write.csv(Unisec, file="Unisec.csv", row.names = F)
Unimaster <- masback2
write.csv(Unimaster, file="Unimaster.csv", row.names = F)


# Converted from UniProtKB AC/ID to UniProt / Gene Name, see tsv file for raw
library(readr)

test1 <- read_tsv('Uniprot_Primary_Return.tsv', show_col_types = F)
test2 <- read_tsv('Uniprot_Secondary_Return.tsv', show_col_types = F)
test3 <- read_tsv('Uniprot_Master_Return.tsv', show_col_types = F)

# Convert all Gene Symbols to Entrez IDs for Visualization Analysis
# Create background denominator for analysis from whole genome
library("clusterProfiler")
require(DOSE)
library(org.Mm.eg.db)

primgenes2 <- bitr(test1$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
secgenes2 <- bitr(test2$To, fromType = "SYMBOL", toType = "ENTREZID", 
                  OrgDb = "org.Mm.eg.db")
masgenes2 <- bitr(test3$To, fromType = "SYMBOL", toType = "ENTREZID", 
                  OrgDb = "org.Mm.eg.db")
unigene <- as.data.frame(org.Mm.egGO)
mmgenome <- unique(sort(unigene$gene_id))


# ID extraction and ORA
primids <- primgenes2$ENTREZID
secids <- secgenes2$ENTREZID
master <- masgenes2$ENTREZID

# KEGG ORA
kegg.prim <- enrichKEGG(gene = primids, universe = mmgenome,
                        organism = "mmu", keyType = "kegg",
                        pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                        qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.prim <- setReadable(kegg.prim, org.Mm.eg.db, keyType = "ENTREZID")

kegg.sec <- enrichKEGG(gene = secids, universe = mmgenome,
                       organism = "mmu", keyType = "kegg",
                       pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                       qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.sec <- setReadable(kegg.sec, org.Mm.eg.db, keyType = "ENTREZID")


# KEGG Prism Figures with BRITE Annotations Exportation
# https://www.genome.jp/brite/query=01100&htext=br08901.keg&option=-a&node_
# proc=br08901_org&proc_enabled=mmu&panel=collapse
write.csv(readkegg.prim, file="KEGGPrim.csv", row.names = F)
write.csv(readkegg.sec, file="KEGGSec.csv", row.names = F)


# BP ORA
ego.bp.prim <- enrichGO(gene = primids, universe =  mmgenome,
                        OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                        pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                        qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                        maxGSSize = 5000)
ego.bp.sec <- enrichGO(gene = secids, universe =  mmgenome,
                       OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                       qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                       maxGSSize = 5000)

# Volcano Plots
library("EnhancedVolcano")

# 1 mM Targets
volp1 <- read_excel("2023_01_12.ExVivo.1mM.MASTER.xlsx")
volp2 <- as.matrix(volp1)
volp2 <- as.data.frame(volp2)
volp3 <- as.matrix(volp2$Group)
volp3 <- gsub(";...*","", volp3)
volp4 <- data.frame(as.numeric(volp2$log2FC), as.numeric(volp2$`Q Value`))
rownames(volp4) <- volp3
volp4 <- as.data.frame(volp4)
colnames(volp4) <- c("log2FC", "Q Value")

p5 <- EnhancedVolcano(volp4, lab = rownames(volp4), x = 'log2FC', y = 'Q Value',
                      pCutoff = 0.05, FCcutoff = 0.58, pointSize = 3.0, labSize = 0,
                      col = c("grey", "grey", "grey", "#44AA99"), colAlpha = 1,
                      xlim = c(-3, 5), xlab = "log2(1mM/0mM)", 
                      ylab = "-log(QValue)", 
                      gridlines.major = F, gridlines.minor = F)

# 5 mM Targets
vols1 <- read_excel("2023_01_12.ExVivo.5mM.MASTER.xlsx")
vols2 <- as.matrix(vols1)
vols2 <- as.data.frame(vols2)
vols3 <- as.matrix(vols2$Group)
vols3 <- gsub(";...*","", vols3)
vols4 <- data.frame(as.numeric(vols2$log2FC), as.numeric(vols2$`Q Value`))
rownames(vols4) <- vols3
vols4 <- as.data.frame(vols4)
colnames(vols4) <- c("log2FC", "Q Value")

p6 <- EnhancedVolcano(vols4, lab = rownames(vols4), x = 'log2FC', y = 'Q Value',
                      pCutoff = 0.05, FCcutoff = 0.58, pointSize = 3.0, labSize = 0,
                      col = c("grey", "grey", "grey", "#117733"), colAlpha = 1,
                      xlim = c(-5, 15), xlab = "log2(5mM/0mM)", 
                      ylab = "-log(QValue)", 
                      gridlines.major = F, gridlines.minor = F)

# Volcano Plots Exported Together
pdf(paste(Sys.Date(),"Volcano_ExVivo.pdf", sep = "_"), width = 15)
grid.arrange(p5, p6, ncol = 2)
dev.off()











##################################################
# IN VIVO DATA
##################################################

# Load data from Christina, I abbreviated, alphabetized by UNIPROT and corrected
# where labeling of gene symbols was incorrect bc of Excel
library(readxl)

f1targ1 <- read_excel("2023_01_04.InVivo.F1.xlsx")
f1targ2 <- as.matrix(f1targ1)
f1targ2 <- as.data.frame(f1targ2)
f1master1 <- read_excel("2023_01_04.InVivo.F1.REPLICATE.xlsx")
f1master2 <- as.matrix(f1master1)
f1master2 <- as.data.frame(f1master2)

f2targ1 <- read_excel("2023_01_23.InVivo.F2.xlsx")
f2targ2 <- as.matrix(f2targ1)
f2targ2 <- as.data.frame(f2targ2)
f2master1 <- read_excel("2023_01_23.InVivo.F2.REPLICATE.xlsx")
f2master2 <- as.matrix(f2master1)
f2master2 <- as.data.frame(f2master2)

f3targ1 <- read_excel("2023_01_23.InVivo.F3.xlsx")
f3targ2 <- as.matrix(f3targ1)
f3targ2 <- as.data.frame(f3targ2)
f3master1 <- read_excel("2023_01_23.InVivo.F3.REPLICATE.xlsx")
f3master2 <- as.matrix(f3master1)
f3master2 <- as.data.frame(f3master2)

f4targ1 <- read_excel("2023_01_04.InVivo.F4.xlsx")
f4targ2 <- as.matrix(f4targ1)
f4targ2 <- as.data.frame(f4targ2)
f4master1 <- read_excel("2023_01_04.InVivo.F4.REPLICATE.xlsx")
f4master2 <- as.matrix(f4master1)
f4master2 <- as.data.frame(f4master2)


# Add replicates to target dfs, clean off unneeded, add log2 and QValue(numeric)
f1targ3 <- f1master2[f1master2$C_UNIPROT %in% f1targ2$C_UNIPROT,]
f1targ3 <- f1targ3[,-1]
f1targ4 <- cbind(f1targ3, as.numeric(f1targ2$`log2(1/0)`),
                 as.numeric(f1targ2$`Q Value`))

f2targ3 <- f2master2[f2master2$C_UNIPROT %in% f2targ2$C_UNIPROT,]
f2targ3 <- f2targ3[,-1]
f2targ4 <- cbind(f2targ3, as.numeric(f2targ2$`log2(1/0)`),
                 as.numeric(f2targ2$`Q Value`))

f3targ3 <- f3master2[f3master2$C_UNIPROT %in% f3targ2$C_UNIPROT,]
f3targ3 <- f3targ3[,-1]
f3targ4 <- cbind(f3targ3, as.numeric(f3targ2$`log2(1/0)`),
                 as.numeric(f3targ2$`Q Value`))

f4targ3 <- f4master2[f4master2$C_UNIPROT %in% f4targ2$C_UNIPROT,]
f4targ3 <- f4targ3[,-1]
f4targ4 <- cbind(f4targ3, as.numeric(f4targ2$`log2(1/0)`),
                 as.numeric(f4targ2$`Q Value`))


# Clean extraneous details from Christina's UNIPROT IDs, name columns for ease
f1genes1 <- as.matrix(f1targ4$C_UNIPROT)
f1genes1 <- gsub(";...*","", f1genes1)
f1targ5 <- cbind(f1targ4,f1genes1)
colnames(f1targ5)[8] <- "log2FC"
colnames(f1targ5)[9] <- "QValue"
colnames(f1targ5)[10] <- "UNIPROT"

f2genes1 <- as.matrix(f2targ4$C_UNIPROT)
f2genes1 <- gsub(";...*","", f2genes1)
f2targ5 <- cbind(f2targ4,f2genes1)
colnames(f2targ5)[10] <- "log2FC"
colnames(f2targ5)[11] <- "QValue"
colnames(f2targ5)[12] <- "UNIPROT"

f3genes1 <- as.matrix(f3targ4$C_UNIPROT)
f3genes1 <- gsub(";...*","", f3genes1)
f3targ5 <- cbind(f3targ4,f3genes1)
colnames(f3targ5)[9] <- "log2FC"
colnames(f3targ5)[10] <- "QValue"
colnames(f3targ5)[11] <- "UNIPROT"

f4genes1 <- as.matrix(f4targ4$C_UNIPROT)
f4genes1 <- gsub(";...*","", f4genes1)
f4targ5 <- cbind(f4targ4,f4genes1)
colnames(f4targ5)[10] <- "log2FC"
colnames(f4targ5)[11] <- "QValue"
colnames(f4targ5)[12] <- "UNIPROT"

# Create background lists with total detectable proteome in fractions
f1back1 <- as.matrix(f1master2$C_UNIPROT)
f1back2 <- gsub(";...*","", f1back1)
colnames(f1back2)[1] <- "UNIPROT"

f2back1 <- as.matrix(f2master2$C_UNIPROT)
f2back2 <- gsub(";...*","", f2back1)
colnames(f2back2)[1] <- "UNIPROT"

f3back1 <- as.matrix(f3master2$C_UNIPROT)
f3back2 <- gsub(";...*","", f3back1)
colnames(f3back2)[1] <- "UNIPROT"

f4back1 <- as.matrix(f4master2$C_UNIPROT)
f4back2 <- gsub(";...*","", f4back1)
colnames(f4back2)[1] <- "UNIPROT"

# Export directional bitr input (Uniprot Accessions) to 
# https://www.uniprot.org/id-mapping to convert
Unif1.up <- f1targ5$UNIPROT[f1targ5$log2FC > 0]
write.csv(Unif1.up, file="Unif1.up.csv", row.names = F)
Unif1.dn <- f1targ5$UNIPROT[f1targ5$log2FC < 0]
write.csv(Unif1.dn, file="Unif1.dn.csv", row.names = F)
Unif1.bk <- f1back2
write.csv(Unif1.bk, file="Unif1.bk.csv", row.names = F)

Unif2.up <- f2targ5$UNIPROT[f2targ5$log2FC > 0]
write.csv(Unif2.up, file="Unif2.up.csv", row.names = F)
Unif2.dn <- f2targ5$UNIPROT[f2targ5$log2FC < 0]
write.csv(Unif2.dn, file="Unif2.dn.csv", row.names = F)
Unif2.bk <- f2back2
write.csv(Unif2.bk, file="Unif2.bk.csv", row.names = F)

Unif3.up <- f3targ5$UNIPROT[f3targ5$log2FC > 0]
write.csv(Unif3.up, file="Unif3.up.csv", row.names = F)
Unif3.dn <- f3targ5$UNIPROT[f3targ5$log2FC < 0]
write.csv(Unif3.dn, file="Unif3.dn.csv", row.names = F)
Unif3.bk <- f3back2
write.csv(Unif3.bk, file="Unif3.bk.csv", row.names = F)

Unif4.up <- f4targ5$UNIPROT[f4targ5$log2FC > 0]
write.csv(Unif4.up, file="Unif4.up.csv", row.names = F)
Unif4.dn <- f4targ5$UNIPROT[f4targ5$log2FC < 0]
write.csv(Unif4.dn, file="Unif4.dn.csv", row.names = F)
Unif4.bk <- f4back2
write.csv(Unif4.bk, file="Unif4.bk.csv", row.names = F)

# Converted from UniProtKB AC/ID to UniProt / Gene Name, see tsv file for raw
library(readr)

return1 <- read_tsv('Uniprot_f1up_Return.tsv', show_col_types = F)
return2 <- read_tsv('Uniprot_f1dn_Return.tsv', show_col_types = F)
returnbk1 <- read_tsv('Uniprot_f1bk_Return.tsv', show_col_types = F)

return3 <- read_tsv('Uniprot_f2up_Return.tsv', show_col_types = F)
return4 <- read_tsv('Uniprot_f2dn_Return.tsv', show_col_types = F)
returnbk2 <- read_tsv('Uniprot_f2bk_Return.tsv', show_col_types = F)

return5 <- read_tsv('Uniprot_f3up_Return.tsv', show_col_types = F)
return6 <- read_tsv('Uniprot_f3dn_Return.tsv', show_col_types = F)
returnbk3 <- read_tsv('Uniprot_f3bk_Return.tsv', show_col_types = F)

return7 <- read_tsv('Uniprot_f4up_Return.tsv', show_col_types = F)
return8 <- read_tsv('Uniprot_f4dn_Return.tsv', show_col_types = F)
returnbk4 <- read_tsv('Uniprot_f4bk_Return.tsv', show_col_types = F)


# Convert all Gene Symbols to Entrez IDs for Visualization Analysis
# Create background denominator for analysis from whole genome
library("clusterProfiler")
require(DOSE)
library(org.Mm.eg.db)

f1genes.up <- bitr(return1$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
f1genes.dn <- bitr(return2$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
f1genes.bk <- bitr(returnbk1$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")

f2genes.up <- bitr(return3$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
f2genes.dn <- bitr(return4$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
f2genes.bk <- bitr(returnbk2$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")

f3genes.up <- bitr(return5$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
f3genes.dn <- bitr(return6$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
f3genes.bk <- bitr(returnbk3$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")

f4genes.up <- bitr(return7$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
f4genes.dn <- bitr(return8$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
f4genes.bk <- bitr(returnbk4$To, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")

unigene <- as.data.frame(org.Mm.egGO)
mmgenome <- unique(sort(unigene$gene_id))

# ID extraction and ORA
f1up <- f1genes.up$ENTREZID
f1dn <- f1genes.dn$ENTREZID
f1bk <- f1genes.bk$ENTREZID

f2up <- f2genes.up$ENTREZID
f2dn <- f2genes.dn$ENTREZID
f2bk <- f2genes.bk$ENTREZID

f3up <- f3genes.up$ENTREZID
f3dn <- f3genes.dn$ENTREZID
f3bk <- f3genes.bk$ENTREZID

f4up <- f4genes.up$ENTREZID
f4dn <- f4genes.dn$ENTREZID
f4bk <- f4genes.bk$ENTREZID

totalbk <- c(f1bk, f2bk, f3bk, f4bk)



##########
# F1 Analysis
##########

# KEGG ORA
kegg.f1.up <- enrichKEGG(gene = f1up, universe = mmgenome,
                        organism = "mmu", keyType = "kegg",
                        pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                        qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.f1.up <- setReadable(kegg.f1.up, org.Mm.eg.db, keyType = "ENTREZID")

kegg.f1.dn <- enrichKEGG(gene = f1dn, universe = mmgenome,
                         organism = "mmu", keyType = "kegg",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.f1.dn <- setReadable(kegg.f1.dn, org.Mm.eg.db, keyType = "ENTREZID")

# BP ORA
ego.bp.f1.up <- enrichGO(gene = f1up, universe =  totalbk,
                        OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                        pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                        qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                        maxGSSize = 5000)
ego.bp.f1.dn <- enrichGO(gene = f1dn, universe =  totalbk,
                         OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                         maxGSSize = 5000)



##########
# F2 Analysis
##########

# KEGG ORA
kegg.f2.up <- enrichKEGG(gene = f2up, universe = mmgenome,
                         organism = "mmu", keyType = "kegg",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.f2.up <- setReadable(kegg.f2.up, org.Mm.eg.db, keyType = "ENTREZID")

kegg.f2.dn <- enrichKEGG(gene = f2dn, universe = mmgenome,
                         organism = "mmu", keyType = "kegg",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.f2.dn <- setReadable(kegg.f2.dn, org.Mm.eg.db, keyType = "ENTREZID")

# BP ORA
ego.bp.f2.up <- enrichGO(gene = f2up, universe =  totalbk,
                         OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                         maxGSSize = 5000)
ego.bp.f2.dn <- enrichGO(gene = f2dn, universe =  totalbk,
                         OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                         maxGSSize = 5000)



##########
# F3 Analysis
##########

# KEGG ORA
kegg.f3.up <- enrichKEGG(gene = f3up, universe = mmgenome,
                         organism = "mmu", keyType = "kegg",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.f3.up <- setReadable(kegg.f3.up, org.Mm.eg.db, keyType = "ENTREZID")

kegg.f3.dn <- enrichKEGG(gene = f3dn, universe = mmgenome,
                         organism = "mmu", keyType = "kegg",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.f3.dn <- setReadable(kegg.f3.dn, org.Mm.eg.db, keyType = "ENTREZID")

# BP ORA
ego.bp.f3.up <- enrichGO(gene = f3up, universe =  totalbk,
                         OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                         maxGSSize = 5000)
ego.bp.f3.dn <- enrichGO(gene = f3dn, universe =  totalbk,
                         OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                         maxGSSize = 5000)



##########
# F4 Analysis
##########

# KEGG ORA
kegg.f4.up <- enrichKEGG(gene = f4up, universe = mmgenome,
                         organism = "mmu", keyType = "kegg",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.f4.up <- setReadable(kegg.f4.up, org.Mm.eg.db, keyType = "ENTREZID")

kegg.f4.dn <- enrichKEGG(gene = f4dn, universe = mmgenome,
                         organism = "mmu", keyType = "kegg",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.f4.dn <- setReadable(kegg.f4.dn, org.Mm.eg.db, keyType = "ENTREZID")

# BP ORA
ego.bp.f4.up <- enrichGO(gene = f4up, universe =  totalbk,
                         OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                         maxGSSize = 5000)
ego.bp.f4.dn <- enrichGO(gene = f4dn, universe =  totalbk,
                         OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, readable = TRUE, minGSSize = 10, 
                         maxGSSize = 5000)



##########
# Volcano Plots
##########
library("EnhancedVolcano")

# Fraction 1
vol1f1 <- read_excel("2023_01_23.InVivo.F1.MASTER.xlsx")
vol1f2 <- as.matrix(vol1f1)
vol1f2 <- as.data.frame(vol1f2)
vol1f3 <- as.matrix(vol1f2$Group)
vol1f4 <- data.frame(as.numeric(vol1f2$log2FC), as.numeric(vol1f2$`Q Value`))
colnames(vol1f4) <- c("log2FC", "Q Value")

p1 <- EnhancedVolcano(vol1f4, lab = rownames(vol1f4), x = 'log2FC', y = 'Q Value',
                      pCutoff = 0.05, FCcutoff = 0.58, pointSize = 3.0, labSize = 0, 
                      col = c("grey", "grey", "grey", "#332288"), colAlpha = 1,
                      xlab = "log2(BH-BD/canola)", ylab = "-log(QValue)", xlim = c(-7, 7),
                      ylim = c(0, 4), gridlines.major = F, gridlines.minor = F)

# Fraction 2
vol2f1 <- read_excel("2023_01_23.InVivo.F2.MASTER.xlsx")
vol2f2 <- as.matrix(vol2f1)
vol2f2 <- as.data.frame(vol2f2)
vol2f3 <- as.matrix(vol2f2$Group)
vol2f4 <- data.frame(as.numeric(vol2f2$log2FC), as.numeric(vol2f2$`Q Value`))
colnames(vol2f4) <- c("log2FC", "Q Value")

p2 <- EnhancedVolcano(vol2f4, lab = rownames(vol2f4), x = 'log2FC', y = 'Q Value', 
                      pCutoff = 0.05, FCcutoff = 0.58, pointSize = 3.0, labSize = 0,
                      col = c("grey", "grey", "grey", "#88CCEE"), colAlpha = 1, 
                      xlab = "log2(BH-BD/canola)", ylab = "-log(QValue)", 
                      xlim = c(-7.5, 3.5), ylim = c(0, 100), 
                      gridlines.major = F, gridlines.minor = F)

# Fraction 3
vol3f1 <- read_excel("2023_01_23.InVivo.F3.MASTER.xlsx")
vol3f2 <- as.matrix(vol3f1)
vol3f2 <- as.data.frame(vol3f2)
vol3f3 <- as.matrix(vol3f2$Group)
vol3f4 <- data.frame(as.numeric(vol3f2$log2FC), as.numeric(vol3f2$`Q Value`))
colnames(vol3f4) <- c("log2FC", "Q Value")

p3 <- EnhancedVolcano(vol3f4, lab = rownames(vol3f4), x = 'log2FC', y = 'Q Value',
                      pCutoff = 0.05, FCcutoff = 0.58, pointSize = 3.0, labSize = 0,
                      col = c("grey", "grey", "grey", "#CC6677"), colAlpha = 1,
                      xlab = "log2(BH-BD/canola)", ylab = "-log(QValue)", 
                      xlim = c(-6.5, 1.5), ylim = c(0, 40), 
                      gridlines.major = F, gridlines.minor = F)

# Fraction 4
vol4f1 <- read_excel("2023_01_23.InVivo.F4.MASTER.xlsx")
vol4f2 <- as.matrix(vol4f1)
vol4f2 <- as.data.frame(vol4f2)
vol4f3 <- as.matrix(vol4f2$Group)
vol4f4 <- data.frame(as.numeric(vol4f2$log2FC), as.numeric(vol4f2$`Q Value`))
colnames(vol4f4) <- c("log2FC", "Q Value")

p4 <- EnhancedVolcano(vol4f4, lab = rownames(vol4f4), x = 'log2FC', y = 'Q Value',
                      pCutoff = 0.05, FCcutoff = 0.58, pointSize = 3.0, labSize = 0,
                      col = c("grey", "grey", "grey", "#DDCC77"), colAlpha = 1, 
                      xlab = "log2(BH-BD/canola)", ylab = "-log(QValue)", xlim = c(-5, 5),
                      ylim = c(0, 35), gridlines.major = F, gridlines.minor = F)

# All Fractions print together
pdf(paste(Sys.Date(),"Volcano_AllFractions.pdf", sep = "_"), width = 20)
grid.arrange(p1, p2, p3, p4, ncol = 4)
dev.off()
















##################################################
# COMBINED DATA
##################################################



##########
# Venn Diagrams and Primary R-BHB Target KEGG ORA
##########
library(ggvenn)

# Ex vivo target comparisons and visualization
venn1 <- list('1/0 mM' = Uniprim, '5/0 mM' = Unisec)

pdf(paste(Sys.Date(),"Venn_Diagram_ExVivo.pdf", sep = "_"), width = 12)
ggvenn(venn1, show_percentage = F, fill_color =  c("#44AA99", "#117733"), 
       text_size = 12, set_name_size = 10)
dev.off()

# Comparison of targets in In Vivo all four fractions, both directions
venn2 <- list('Fraction 1' = c(Unif1.up, Unif1.dn), 
              'Fraction 2' = c(Unif2.up, Unif2.dn),
              'Fraction 3' = c(Unif3.up, Unif3.dn), 
              'Fraction 4' = c(Unif4.up, Unif4.dn))

pdf(paste(Sys.Date(),"Venn_Diagram_InVivo.pdf", sep = "_"))
ggvenn(venn2, show_percentage = F, 
       fill_color =  c("#332288", "#88CCEE", "#CC6677", "#DDCC77"),
       text_size = 10, set_name_size = 8)
dev.off()

# In vivo and Ex vivo comparisons and visualization
venn3 <- list('BH-BD Targets' = c(Unif1.up, Unif2.up, Unif3.up, Unif4.up), 
              '1/0 mM' = Uniprim, '5/0 mM' = Unisec)

pdf(paste(Sys.Date(),"Venn_Diagram_ComparisonExIn.pdf", sep = "_"), width = 12)
ggvenn(venn3, show_percentage = F, 
       fill_color =  c("#332288", "#44AA99", "#117733"), 
       text_size = 10)
dev.off()

# Extract and export the lists of intersecting proteins from the Venn Diagram
# Exported to Uniprot like in previous section
vennlist1 <- intersect(c(Unif1.up, Unif2.up, Unif3.up, Unif4.up), Uniprim)
vennlist2 <- intersect(vennlist1, Unisec)
write.csv(vennlist1, file="Primary Binding Targets.csv", row.names = F)

# Return lists and convert using bitr for KEGG ORA
library(readr)

preturn <- read_tsv('Primary_Binding_Target_Return.tsv', show_col_types = F)
primvenn1 <- bitr(preturn$To, fromType = "SYMBOL", toType = "ENTREZID", 
                  OrgDb = "org.Mm.eg.db")

# Run KEGG ORA
primvenn2 <- primvenn1$ENTREZID

kegg.vprim <- enrichKEGG(gene = primvenn2, universe = mmgenome,
                         organism = "mmu", keyType = "kegg",
                         pAdjustMethod = "BH", pvalueCutoff  = 0.01, 
                         qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 5000)
readkegg.vprim <- setReadable(kegg.vprim, org.Mm.eg.db, keyType = "ENTREZID")



##########
# KEGG Barplots
##########
library(ggplot2)
library(readxl)
library(gridExtra)
library(stringr)

# Export KEGG Data and amend for only top 10 names and log2(qvalues)
write.csv(readkegg.f1.up, file="KEGGF1up.csv", row.names = F)
write.csv(readkegg.f2.dn, file="KEGGF2dn.csv", row.names = F)
write.csv(readkegg.f2.up, file="KEGGF2up.csv", row.names = F)
write.csv(readkegg.f3.dn, file="KEGGF3dn.csv", row.names = F)
write.csv(readkegg.f3.up, file="KEGGF3up.csv", row.names = F)
write.csv(readkegg.f4.dn, file="KEGGF4dn.csv", row.names = F)
write.csv(readkegg.vprim, file="KEGGVenn.csv", row.names = F)

# Return and name anew
df1d <- read_csv("KEGGF1up_v2.csv")
df1d <- as.data.frame(df1d)
colnames(df1d) <- c("mmu", "QValue")
df2d <- read_csv("KEGGF2dn_v2.csv")
df2d <- as.data.frame(df2d)
colnames(df2d) <- c("mmu", "QValue")
df2u <- read_csv("KEGGF2up_v2.csv")
df2u <- as.data.frame(df2u)
colnames(df2u) <- c("mmu", "QValue")
df3d <- read_csv("KEGGF3dn_v2.csv")
df3d <- as.data.frame(df3d)
colnames(df3d) <- c("mmu", "QValue")
df3u <- read_csv("KEGGF3up_v2.csv")
df3u <- as.data.frame(df3u)
colnames(df3u) <- c("mmu", "QValue")
df4d <- read_csv("KEGGF4dn_v2.csv")
df4d <- as.data.frame(df4d)
colnames(df4d) <- c("mmu", "QValue")
dfvn <- read_csv("KEGGVenn_v2.csv")
dfvn <- as.data.frame(dfvn)
colnames(dfvn) <- c("mmu", "QValue")

# Load figures
f1 <- ggplot(df1d, aes(x=mmu, y=QValue)) +
      geom_bar(stat = "identity", color = '#000000',
               fill = '#332288', width = 0.5) +
      coord_flip() +
      ylab("-log(QValue)") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                       limits=df1d$mmu) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.line = element_line(colour = '#000000'),
            axis.ticks = element_line(colour = '#000000'))
f2 <- ggplot(df2d, aes(x=mmu, y=QValue)) +
      geom_bar(stat = "identity", color = '#000000',
               fill = '#88CCEE', width = 0.5) +
      coord_flip() +
      ylab("-log(QValue)") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                       limits=df2d$mmu) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.line = element_line(colour = '#000000'),
            axis.ticks = element_line(colour = '#000000'))
f3 <- ggplot(df2u, aes(x=mmu, y=QValue)) +
      geom_bar(stat = "identity", color = '#000000',
               fill = '#88CCEE', width = 0.5) +
      coord_flip() +
      ylab("-log(QValue)") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                       limits=df2u$mmu) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.line = element_line(colour = '#000000'),
            axis.ticks = element_line(colour = '#000000'))
f4 <- ggplot(df3d, aes(x=mmu, y=QValue)) +
      geom_bar(stat = "identity", color = '#000000',
               fill = '#CC6677', width = 0.5) +
      coord_flip() +
      ylab("-log(QValue)") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                       limits=df3d$mmu) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.line = element_line(colour = '#000000'),
            axis.ticks = element_line(colour = '#000000'))
f5 <- ggplot(df3u, aes(x=mmu, y=QValue)) +
      geom_bar(stat = "identity", color = '#000000', 
               fill = '#CC6677', width = 0.5) +
      coord_flip() +
      ylab("-log(QValue)") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                       limits=df3u$mmu) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.line = element_line(colour = '#000000'),
            axis.ticks = element_line(colour = '#000000'))
f6 <- ggplot(df4d, aes(x=mmu, y=QValue)) +
      geom_bar(stat = "identity", color = '#000000', 
               fill = '#DDCC77', width = 0.5) +
      coord_flip() +
      ylab("-log(QValue)") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                       limits=df4d$mmu) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.line = element_line(colour = '#000000'),
            axis.ticks = element_line(colour = '#000000'))
f7 <- ggplot(dfvn, aes(x=mmu, y=QValue)) +
      geom_bar(stat = "identity", color = '#000000', 
               fill = '#6699CC', width = 0.5) +
      coord_flip() +
      ylab("-log(QValue)") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                       limits=dfvn$mmu) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.line = element_line(colour = '#000000'),
            axis.ticks = element_line(colour = '#000000'))

# Print to pdf
# Set 1 Main Figures
pdf(paste(Sys.Date(), "KEGG Dotplots Main", sep = "_"), width = 30, height = 8)
grid.arrange(f3, f5, f6, ncol = 3)
dev.off()
# Set 2 Supplementary Figures
pdf(paste(Sys.Date(), "KEGG Dotplots Supplementary", sep = "_"), width = 30, 
    height = 5)
grid.arrange(f1, f2, f4, ncol = 3)
dev.off()
# Set 3 Main Figure
pdf(paste(Sys.Date(), "KEGG Dotplots Venn Primary", sep = "_"), width = 9, 
    height = 7)
f7
dev.off()



##########
# GO BP Barplots
##########
library(ggplot2)
library(readxl)
library(gridExtra)
library(stringr)

# Export GO BP Data and amend for only top 10 names and log2(qvalues)
write.csv(ego.bp.prim, file="GOBPprim.csv", row.names = F)
write.csv(ego.bp.sec, file="GOBPsec.csv", row.names = F)

write.csv(ego.bp.f2.dn, file="GOBPF2dn.csv", row.names = F)
write.csv(ego.bp.f2.up, file="GOBPF2up.csv", row.names = F)
write.csv(ego.bp.f3.dn, file="GOBPF3dn.csv", row.names = F)
write.csv(ego.bp.f3.up, file="GOBPF3up.csv", row.names = F)
write.csv(ego.bp.f4.dn, file="GOBPF4dn.csv", row.names = F)

# Return and name anew
dtfprim <- read_csv("GOBPprim_v2.csv")
dtfprim <- as.data.frame(dtfprim)
colnames(dtfprim) <- c("go", "QValue")
dtfsec <- read_csv("GOBPsec_v2.csv")
dtfsec <- as.data.frame(dtfsec)
colnames(dtfsec) <- c("go", "QValue")

dtf2d <- read_csv("GOBPF2dn_v2.csv")
dtf2d <- as.data.frame(dtf2d)
colnames(dtf2d) <- c("go", "QValue")
dtf2u <- read_csv("GOBPF2up_v2.csv")
dtf2u <- as.data.frame(dtf2u)
colnames(dtf2u) <- c("go", "QValue")
dtf3d <- read_csv("GOBPF3dn_v2.csv")
dtf3d <- as.data.frame(dtf3d)
colnames(dtf3d) <- c("go", "QValue")
dtf3u <- read_csv("GOBPF3up_v2.csv")
dtf3u <- as.data.frame(dtf3u)
colnames(dtf3u) <- c("go", "QValue")
dtf4d <- read_csv("GOBPF4dn_v2.csv")
dtf4d <- as.data.frame(dtf4d)
colnames(dtf4d) <- c("go", "QValue")

# Load figures
fg1 <- ggplot(dtfprim, aes(x=go, y=QValue)) +
        geom_bar(stat = "identity", color = '#000000',
                 fill = '#44AA99', width = 0.5) +
        coord_flip() +
        ylab("-log(QValue)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                         limits=dtfprim$go) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.line = element_line(colour = '#000000'),
              axis.ticks = element_line(colour = '#000000'))
fg2 <- ggplot(dtfsec, aes(x=go, y=QValue)) +
        geom_bar(stat = "identity", color = '#000000',
                 fill = '#117733', width = 0.5) +
        coord_flip() +
        ylab("-log(QValue)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                         limits=dtfsec$go) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.line = element_line(colour = '#000000'),
              axis.ticks = element_line(colour = '#000000'))
fg3 <- ggplot(dtf2d, aes(x=go, y=QValue)) +
        geom_bar(stat = "identity", color = '#000000',
                 fill = '#88CCEE', width = 0.5) +
        coord_flip() +
        ylab("-log(QValue)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                         limits=dtf2d$go) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.line = element_line(colour = '#000000'),
              axis.ticks = element_line(colour = '#000000'))
fg4 <- ggplot(dtf2u, aes(x=go, y=QValue)) +
        geom_bar(stat = "identity", color = '#000000',
                 fill = '#88CCEE', width = 0.5) +
        coord_flip() +
        ylab("-log(QValue)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                         limits=dtf2u$go) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.line = element_line(colour = '#000000'),
              axis.ticks = element_line(colour = '#000000'))
fg5 <- ggplot(dtf3d, aes(x=go, y=QValue)) +
        geom_bar(stat = "identity", color = '#000000',
                 fill = '#CC6677', width = 0.5) +
        coord_flip() +
        ylab("-log(QValue)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                         limits=dtf3d$go) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.line = element_line(colour = '#000000'),
              axis.ticks = element_line(colour = '#000000'))
fg6 <- ggplot(dtf3u, aes(x=go, y=QValue)) +
        geom_bar(stat = "identity", color = '#000000', 
                 fill = '#CC6677', width = 0.5) +
        coord_flip() +
        ylab("-log(QValue)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                         limits=dtf3u$go) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.line = element_line(colour = '#000000'),
              axis.ticks = element_line(colour = '#000000'))
fg7 <- ggplot(dtf4d, aes(x=go, y=QValue)) +
        geom_bar(stat = "identity", color = '#000000', 
                 fill = '#DDCC77', width = 0.5) +
        coord_flip() +
        ylab("-log(QValue)") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30), 
                         limits=dtf4d$go) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.line = element_line(colour = '#000000'),
              axis.ticks = element_line(colour = '#000000'))

# Print to pdf
# Set 1 Supplementary Figures Ex Vivo
pdf(paste(Sys.Date(), "GO BP Dotplots Ex Vivo", sep = "_"), width = 15,
    height = 7)
grid.arrange(fg1, fg2, ncol = 2)
dev.off()
# Set 2 Supplementary Figures In Vivo
pdf(paste(Sys.Date(), "GO BP Dotplots In Vivo", sep = "_"), width = 30, 
    height = 14)
grid.arrange(fg3, fg4, fg5, fg6, fg7, ncol = 3)
dev.off()



##########
# PLS-DA
##########
library(mixOmics)

# Ex Vivo Total
ValExT <- read.csv("2023_05_16.ExVivo1.PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)
GrpExT <- read.csv("2023_05_16.ExVivo2.PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)

log.ValExT <- data.frame(sapply(ValExT[,c(2:12)], log2))
t.log.ValExT <- data.frame(t(log.ValExT))
GrpExT.F <- factor(GrpExT$Group)
plsda1 <- plsda(t.log.ValExT, GrpExT.F, ncomp = 11)

pdf(paste(Sys.Date(), "PLSDA Ex Vivo", sep = "_"), width = 5, height = 5.5)
plotIndiv(plsda1, ind.names = F, legend = T, ellipse = T, title = 'Ex Vivo',
          col.per.group = c("#999999", "#44AA99", "#117733"), 
          legend.position = "bottom")
dev.off()

# In Vivo Fraction 1
ValIn1 <- read.csv("2023_05_16.InVivoF1.1PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)
GrpIn1 <- read.csv("2023_05_16.InVivoF1.2PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)

log.ValIn1 <- data.frame(sapply(ValIn1[,c(2:6)], log2))
t.log.ValIn1 <- data.frame(t(log.ValIn1))
GrpIn1.F <- factor(GrpIn1$Group)
plsda2 <- plsda(t.log.ValIn1, GrpIn1.F, ncomp = 5)

pdf(paste(Sys.Date(), "PLSDA Fraction 1", sep = "_"), width = 5, height = 5)
plotIndiv(plsda2, ind.names = F, legend = T, ellipse = T, title = 'Fraction 1',
          col.per.group = c("#6699CC", "#999999"))
dev.off()

# In Vivo Fraction 2
ValIn2 <- read.csv("2023_05_16.InVivoF2.1PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)
GrpIn2 <- read.csv("2023_05_16.InVivoF2.2PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)

log.ValIn2 <- data.frame(sapply(ValIn2[,c(2:8)], log2))
t.log.ValIn2 <- data.frame(t(log.ValIn2))
GrpIn2.F <- factor(GrpIn2$Group)
plsda3 <- plsda(t.log.ValIn2, GrpIn2.F, ncomp = 7)

pdf(paste(Sys.Date(), "PLSDA Fraction 2", sep = "_"), width = 5, height = 5)
plotIndiv(plsda3, ind.names = F, legend = T, ellipse = T, title = 'Fraction 2',
          col.per.group = c("#6699CC", "#999999"))
dev.off()

# In Vivo Fraction 3
ValIn3 <- read.csv("2023_05_16.InVivoF3.1PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)
GrpIn3 <- read.csv("2023_05_16.InVivoF3.2PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)

log.ValIn3 <- data.frame(sapply(ValIn3[,c(2:7)], log2))
t.log.ValIn3 <- data.frame(t(log.ValIn3))
GrpIn3.F <- factor(GrpIn3$Group)
plsda4 <- plsda(t.log.ValIn3, GrpIn3.F, ncomp = 6)

pdf(paste(Sys.Date(), "PLSDA Fraction 3", sep = "_"), width = 5, height = 5)
plotIndiv(plsda4, ind.names = F, legend = T, ellipse = T, title = 'Fraction 3',
          col.per.group = c("#6699CC", "#999999"))
dev.off()

# In Vivo Fraction 4
ValIn4 <- read.csv("2023_05_16.InVivoF4.1PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)
GrpIn4 <- read.csv("2023_05_16.InVivoF4.2PLSDA.csv", sep = ",", dec = ".",
                   header = T, na = "NA", stringsAsFactors = F)

log.ValIn4 <- data.frame(sapply(ValIn4[,c(2:8)], log2))
t.log.ValIn4 <- data.frame(t(log.ValIn4))
GrpIn4.F <- factor(GrpIn4$Group)
plsda5 <- plsda(t.log.ValIn4, GrpIn4.F, ncomp = 7)

pdf(paste(Sys.Date(), "PLSDA Fraction 4", sep = "_"), width = 5, height = 5)
plotIndiv(plsda5, ind.names = F, legend = T, ellipse = T, title = 'Fraction 4',
          col.per.group = c("#6699CC", "#999999"))
dev.off()

