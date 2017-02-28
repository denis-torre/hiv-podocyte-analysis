#################################################################
#################################################################
############### HIV Podocyte Analysis Pipeline ##################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
source('/Users/denis/Documents/Projects/scripts/Support.R')

##### 2. Other libraries #####

#######################################################
#######################################################
########## S3. Remove Batch Effects
#######################################################
#######################################################

#############################################
########## 1. Remove batch effects
#############################################

remove_batch_effects <- function(vstDataframe, annotationDataframe, outfile) {

	# Load libraries
	library(sva)
	
	# Create design matrix
	treatmentDesign <- model.matrix(~treatment, data=annotationDataframe)

	# Get batch
	patient_id <- annotationDataframe$patient_id

	# Correct data
	vstDataframeCorrected <- ComBat(dat=vstDataframe, batch=patient_id, mod=treatmentDesign, par.prior=TRUE, prior.plots=FALSE)
	vstDataframeCorrected <- cbind(gene_symbol=rownames(vstDataframeCorrected), vstDataframeCorrected)

	# Save file
	write.table(vstDataframeCorrected, file=outfile, sep='\t', quote=FALSE, row.names=FALSE)
}

#######################################################
#######################################################
########## S2. Run PCA
#######################################################
#######################################################

#############################################
########## 1. Run PCA
#############################################

run_pca <- function(infiles, outfile) {
	
	# Read data
	expressionDataframe <- read.table(infiles[[1]], sep='\t', header=TRUE, row.names='gene_symbol')
	sampleAnnotationDataframe <- read.table(infiles[[2]], sep='\t', header=TRUE, row.names='sample_id')

	# Remove columns
	expressionDataframe <- expressionDataframe[,setdiff(colnames(expressionDataframe), c('NK1', 'NK2', 'B8N', 'B10C'))]

	# Run PCA
	pca <- runPCA(expressionDataframe)

	# Get PCs
	PCs <- c('PC1', 'PC2', 'PC3')

	# Get PCA dataframe
	pcaDataframe <- as.data.frame(pca$x[,PCs])

	# Add variance labels
	pcaDataframe <- rbind(varLabels=pca$varLabels[colnames(pcaDataframe)], pcaDataframe)

	# Merge with sample annotations
	mergedPcaDataframe <- merge(pcaDataframe, sampleAnnotationDataframe, by='row.names', all.x=TRUE)
	colnames(mergedPcaDataframe)[1] <- 'sample_id'

	# Save
	write.table(mergedPcaDataframe, file=outfile, sep='\t', quote=FALSE, row.names=FALSE)
}

#######################################################
#######################################################
########## S4. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Characteristic Direction
#############################################

run_characteristic_direction <- function(expressionDataframe, outfile) {

	# Read data
	source('pipeline/scripts/support/nipals.R')
	source('pipeline/scripts/support/chdir.R')

	# Get timepoints
	timepoints <- unique(sapply(colnames(expressionDataframe), function(x) strsplit(x, '.', fixed=TRUE)[[1]][1]))
	    
	# Get column names grouped by timepoint
	timepointColnames <- sapply(timepoints, function(x) colnames(expressionDataframe)[grepl(x, colnames(expressionDataframe))])
	    
	# Split and filter data
	timepointData <- sapply(timepointColnames, function(x) {
	    resultDataframe <- round(expressionDataframe[,x], digits=2);
	    geneVar <- apply(resultDataframe, 1, var);
	    resultDataframe <- resultDataframe[names(geneVar)[geneVar > 0],];
	    return(resultDataframe);
	}, simplify=FALSE)
	    
	# Get common genes
	commonGenes <- Reduce(intersect, sapply(timepointData, function(x) rownames(x)))
	    
	# Filter data
	timepointData <- sapply(timepointData, function(x) x[commonGenes,])
	    
	# Run CD
	cdList <- sapply(timepoints[-1], function(x) chdir(timepointData[['h0']], timepointData[[x]], commonGenes)[commonGenes,])
	
	# Fix dimension names
	cdList <- cdList[,c('h6', 'h12', 'h24', 'h48')]
	cdList <- cbind(gene_symbol=rownames(cdList), cdList)
	    
	# Save
	write.table(cdList, file=outfile, sep='\t', quote=FALSE, row.names=FALSE)
}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################