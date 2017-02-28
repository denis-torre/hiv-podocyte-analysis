#################################################################
#################################################################
############### HIV Podocyte Analysis Pipeline ##################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, rpy2, os
import pandas as pd
import rpy2.robjects as robjects
import pandas.rpy.common as com

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineHivPodocyteAnalysis as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
rawDataFile = 'rawdata.dir/new_norm_geneData.txt'

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-hiv-podocyte-analysis.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Process Data
#######################################################
#######################################################

#############################################
########## 1. Get rawcount table
#############################################


@follows(mkdir('f1-expression.dir'))

@files(rawDataFile,
	   'f1-expression.dir/hiv_podocyte-rawcounts.txt')

def getRawcountTable(infile, outfile):

	# Read table
	expressionDataframe = pd.read_table(infile)

	# Define columns to keep
	columnsToKeep = ['Symbol', 'B10C', 'B11C', 'B11G', 'B11N', 'B12G', 'B12N', 'B1C', 'B1G', 'B1N', 'B2C', 'B2G', 'B2N', 'B3C', 'B3G', 'B3N', 'B8N', 'B9C', 'B9G', 'B9N', 'NK1', 'NK2']

	# Get subset and rename column
	rawcountDataframe = expressionDataframe[columnsToKeep].set_index('Symbol').dropna()

	# Convert to integer
	rawcountDataframe = rawcountDataframe.astype('int')

	# Write file
	rawcountDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S2. Annotation
#######################################################
#######################################################

#############################################
########## 1. Get sample annotations
#############################################

@follows(mkdir('f2-annotation.dir'))

@transform(getRawcountTable,
		   regex(r'.*/(.*)-rawcounts.txt'),
		   r'f2-annotation.dir/\1-annotation.txt')

def getAnnotationTable(infile, outfile):

	# Read infile
	rawcountDataframe = pd.read_table(infile, index_col='gene_symbol')

	# Get sample names
	sampleNames = rawcountDataframe.columns.tolist()

	# Define attribute dict
	annotationDict = {x:{} for x in sampleNames}

	# Get treatment dict
	treatmentDict = {'N': 'hiv_infected', 'G': 'gfp_infected', 'C': 'control', '1': 'control', '2': 'control'}

	# Loop through sample names
	for sampleName in sampleNames:

	    # Get attributes
	    annotationDict[sampleName]['patient'] = sampleName[:-1]
	    annotationDict[sampleName]['treatment'] = treatmentDict[sampleName[-1]]

	# Convert to dataframe
	annotationDataframe = pd.DataFrame(annotationDict).T

	# Save file
	annotationDataframe.to_csv(outfile, sep='\t', index_label='sample_name')

#######################################################
#######################################################
########## S3. Normalize Expression
#######################################################
#######################################################

#############################################
########## 1. Variance Stabilizing Transform
#############################################

@follows(mkdir('f3-normalized_expression.dir'))

@transform(getRawcountTable,
		   regex(r'.*/(.*)-rawcounts.txt'),
		   r'f3-normalized_expression.dir/\1-vst.txt')

def runVst(infile, outfile):

	# Read expression dataframe
	rawcountDataframe = pd.read_table(infile, index_col='gene_symbol').drop(['NK1', 'NK2'], axis=1)

	# Run function
	vstMatrix = r.runVST(com.convert_to_r_dataframe(rawcountDataframe))

	# Convert to dataframe
	vstDataframe = com.convert_robj(vstMatrix)

	# Write file
	vstDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S4. Batch Effect Removal
#######################################################
#######################################################

#############################################
########## 1. Run ComBat
#############################################

@follows(mkdir('f4-combat.dir'))

@transform(runVst,
		   regex(r'.*/(.*)-vst.txt'),
		   add_inputs(getAnnotationTable),
		   r'f4-combat.dir/\1-combat.txt')

def runComBat(infiles, outfile):

	# Split infiles
	vstFile, annotationFile = infiles

	# Read expression dataframe
	vstDataframe = pd.read_table(vstFile, index_col='gene_symbol').drop(['B8N', 'B10C'], axis=1)

	# Read annotation dataframe
	annotationDataframe = pd.read_table(annotationFile, index_col='sample_name')

	# Get common samples
	annotationDataframe = annotationDataframe.loc[vstDataframe.columns]

	# Run function
	combatMatrix = r.runComBat(com.convert_to_r_dataframe(vstDataframe), com.convert_to_r_dataframe(annotationDataframe), covariateFormula='~treatment', batchColumn='patient')

	# Convert to dataframe
	combatDataframe = com.convert_robj(combatMatrix)

	# Write file
	combatDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S5. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Characteristic Direction
#############################################

@follows(mkdir('f5-characteristic_direction.dir'))

@transform(runComBat,
		   regex(r'.*/(.*).txt'),
		   add_inputs(getAnnotationTable),
		   r'f5-characteristic_direction.dir/\1_cd.txt')

def runCharacteristicDirection(infiles, outfile):

	# Split infiles
	expressionFile, annotationFile = infiles

	# Read expression file
	expressionDataframe = pd.read_table(expressionFile, index_col='gene_symbol')

	# Read annotation file
	annotationDataframe = pd.read_table(annotationFile, index_col='sample_name')

	# Get common samples
	annotationDataframe = annotationDataframe.loc[expressionDataframe.columns]

	# Get treatment dict
	treatmentDict = {x:annotationDataframe.index[annotationDataframe['treatment'] == x].tolist() for x in set(annotationDataframe['treatment'])}

	# Get comparisons
	comparisons = [['control', 'gfp_infected'], ['control', 'hiv_infected'], ['gfp_infected', 'hiv_infected']]

	# Initialize empty dataframe
	resultDataframe = pd.DataFrame()

	# Loop through comparisons
	for comparison in comparisons:
	    
		# Get columns
		controlColumns = treatmentDict[comparison[0]]
		experimentColumns = treatmentDict[comparison[1]]

		# Run characteristic direction
		cdResults = r.runCharacteristicDirection(com.convert_to_r_dataframe(expressionDataframe), experimentColumns, controlColumns)

		# Convert to dataframe
		cdDataframe = com.convert_robj(cdResults).reset_index()

		# Get comparison string
		comparisonString = '_v_'.join(comparison)

		# Add timepoint column
		cdDataframe['comparison'] = comparisonString

		# Append
		resultDataframe = pd.concat([resultDataframe, cdDataframe])

	# Pivot
	resultDataframeCast = resultDataframe.pivot(index='index', columns='comparison', values='CD')

	# Save
	resultDataframeCast.to_csv(outfile, sep='\t', index_label='gene_symbol')


#######################################################
#######################################################
########## S6. Enrichment Analysis
#######################################################
#######################################################

#############################################
########## 1. Submit Genesets 
#############################################

@follows(mkdir('f6-enrichr.dir'))

@transform(runCharacteristicDirection,
		   regex(r'.*/(.*)cd.txt'),
		   r'f6-enrichr.dir/\1enrichr_links.txt')

def submitEnrichrGenesets(infile, outfile):

	# Read infile
	cdDataframe = pd.read_table(infile, index_col='gene_symbol').fillna(0)

	# Initialize link dataframe
	resultDataframe = pd.DataFrame()

	# Loop through timepoints
	for comparison in cdDataframe.columns:

	    # Get Enrichr links
	    enrichrLinkDataframe = S.uploadToEnrichr(cdDataframe, comparison)

	    # Add comparison label
	    enrichrLinkDataframe['comparison'] = comparison

	    # Concatenate
	    resultDataframe = pd.concat([resultDataframe, enrichrLinkDataframe])

	# Save data
	resultDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Get Enrichment results
#############################################

@transform(submitEnrichrGenesets,
		   suffix('links.txt'),
		   'results.txt')

def getEnrichrResults(infile, outfile):

	# Read infile
	enrichrLinkDataframe = pd.read_table(infile, index_col=['geneset','comparison'])

	# Initialize result dataframe
	resultDataframe = pd.DataFrame()

	# Set libraries
	libraries = ['ChEA_2016', 'KEGG_2016', 'GO_Biological_Process_2015', 'GO_Cellular_Component_2015', 'GO_Molecular_Function_2015', 'VirusMINT']

	# Loop through timepoints, genesets and libraries
	for geneset in enrichrLinkDataframe.index.levels[0]:
	    for comparison in enrichrLinkDataframe.index.levels[1]:
	        for library in libraries:

	            # Get enrichment results
	            enrichmentResultDataframe = S.getEnrichmentResults(enrichrLinkDataframe.loc[(geneset, comparison), 'userListId'], library)

	            # Add labels
	            enrichmentResultDataframe['comparison'] = comparison
	            enrichmentResultDataframe['geneset'] = geneset
	            enrichmentResultDataframe['library'] = library

	            # Concatenate
	            resultDataframe = pd.concat([resultDataframe, enrichmentResultDataframe])

    # Write file
	resultDataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S7. L1000CDS2 Analysis
#######################################################
#######################################################

#############################################
########## 1. Submit analysis
#############################################

@follows(mkdir('f7-l1000cds2.dir'), getEnrichrResults)

@transform(runCharacteristicDirection,
		   regex(r'.*/(.*)cd.txt'),
		   r'f7-l1000cds2.dir/\1l1000cds2_links.txt')

def runL1000CDS2(infile, outfile):

	# Read infile
	cdDataframe = pd.read_table(infile, index_col='gene_symbol').fillna(0)

	# Initialize dataframes
	linkDataframe = pd.DataFrame()
	signatureDataframe = pd.DataFrame()

	# Loop through timepoints
	for comparison in cdDataframe.columns:

	    # Run L1000CDS2
	    resultDict = S.getL1000CDS2Results(cdDataframe, comparison)

	    # Add comparison labels
	    resultDict['links']['comparison'] = comparison
	    resultDict['signatures']['comparison'] = comparison

	    # Append dataframes
	    linkDataframe = pd.concat([linkDataframe, resultDict['links']])
	    signatureDataframe = pd.concat([signatureDataframe, resultDict['signatures']])

	# Write files
	linkDataframe.to_csv(outfile, sep='\t', index=False)
	signatureDataframe.to_csv(outfile.replace('links', 'signatures'), sep='\t', index=False)


#############################################
########## 2. Get annotation table
#############################################

@transform(rawDataFile,
		   regex(r'.*/(.*)vst.txt'),
		   r'f1-expression_data.dir/\1sample_annotation.txt')

def getAnnotationData(infile, outfile):

	# Read data
	vstDataframe = pd.read_table(infile).set_index('gene_symbol').drop(['NK1', 'NK2'], axis=1)

	# Get treatment dict
	treatmentDict = {'C': 'Control', 'G': 'GFP-infected', 'N': 'HIV-infected'}

	# Create annotation dataframe
	annotationDataframe = pd.DataFrame([[x, x[:-1], treatmentDict[x[-1]]] for x in vstDataframe.columns], columns=['sample_id','patient_id','treatment'])

	# Save file
	annotationDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Remove batch effects
#############################################

@transform(getAnnotationData,
		   suffix('.txt'),
		   '_corrected.txt')

def removeBatchEffects(infile, outfile):

	# Read data
	vstDataframe = pd.read_table(infile).set_index('gene_symbol').drop(['NK1', 'NK2', 'B8N', 'B10C'], axis=1)

	# Create annotation dataframe
	annotationDataframe = pd.DataFrame([[x, x[:-1], x[-1]] for x in vstDataframe.columns], columns=['sample_id','patient_id','treatment'])

	# Run function
	r.remove_batch_effects(com.convert_to_r_dataframe(vstDataframe), com.convert_to_r_dataframe(annotationDataframe), outfile)

#######################################################
#######################################################
########## S3. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
