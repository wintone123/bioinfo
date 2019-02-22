library(ballgown)
library(genefilter)
library(dplyr)

pheno_data = read.csv("geuvadis_phenodata.csv") 
bg_chr = ballgown(dataDir = "ballgown", samplePattern = "SRR", pData = pheno_data) 
bg_chr_filt = subset(bg_chr, "rowVars(texpr(bg_chr)) > 1", genomesubset = TRUE)
results_transcripts = stattest(bg_chr_filt, feature = "transcript", covariate = "condition", 
                               adjustvars = c("repeat."), getFC = TRUE, meas = "FPKM")
results_transcripts = data.frame(geneNames = ballgown::geneNames(bg_chr_filt),
                                 geneIDs = ballgown::geneIDs(bg_chr_filt), results_transcripts)
write.csv(results_transcripts, "results_transcripts.csv")

results_genes = stattest(bg_chr_filt, feature = "gene", covariate = "condition",
                         adjustvars = c("repeat."), getFC = TRUE, meas = "FPKM")
write.csv(results_genes, "results_genes.csv")

transcript_FPKM_unfilt = texpr(bg_chr,"FPKM")
transcript_FPKM_unfilt = data.frame(indexes(bg_chr)$t2g, transcript_FPKM_unfilt)
transcript_FPKM_unfilt = data.frame(geneName = (ballgown::geneNames(bg_chr)), transcript_FPKM_unfilt)
write.csv(transcript_FPKM_unfilt, "transcript_FPKM_unfilt.csv")

# indices <- match(results_genes$id, texpr(bg, 'all')$gene_id)
# gene_names_for_result <- texpr(bg, 'all')$gene_name[indices]
# results_genes <- data.frame(geneNames=gene_names_for_result, results_genes)