#################################################################
#################################################################
############### HIV Podocyte Analysis Support ###################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
import sys, requests, json
import pandas as pd

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
import Support as S 

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

#######################################################
#######################################################
########## S1. Enrichr Functions
#######################################################
#######################################################

#############################################
########## 1. Add gene lists 
#############################################

def addGeneLists(geneset):
	ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
	genes_str = '\n'.join(geneset)
	payload = {
	    'list': (None, genes_str),
	}
	response = requests.post(ENRICHR_URL, files=payload)
	if not response.ok:
	    raise Exception('Error analyzing gene list')
	data = json.loads(response.text)
	return data

#############################################
########## 2. Get enrichment results
#############################################

def getEnrichmentResults(user_list_id, gene_set_library='GO_Biological_Process_2015'):
	ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
	query_string = '?userListId=%s&backgroundType=%s'
	response = requests.get(
	    ENRICHR_URL + query_string % (user_list_id, gene_set_library)
	 )
	if not response.ok:
	    raise Exception('Error fetching enrichment results')

	data = json.loads(response.text)
	resultDataframe = pd.DataFrame(data[gene_set_library], columns=['rank', 'term_name', 'pvalue', 'zscore', 'combined_score', 'overlapping_genes', 'FDR', 'old_pvalue', 'old_FDR'])
	resultDataframe = resultDataframe.loc[:,['rank','term_name','zscore','combined_score','FDR']]
	return resultDataframe

#######################################################
#######################################################
########## S2. L1000CS2 Analysis
#######################################################
#######################################################

#############################################
########## 1. Run L1000CDS2
#############################################

def runL1000CDS2(differentialExpressionDataframe, column, aggravate=False):
    # Set data
    data = {"genes": differentialExpressionDataframe.index.tolist(), "vals":differentialExpressionDataframe[column].tolist()}
    data['genes'] = [x.upper() for x in data['genes']]
    
    # Set configuration
    config = {"aggravate":aggravate, "searchMethod":"CD", "share":True, "combination":True, "db-version":"latest"}
    payload = {"data":data,"config":config}
    headers = {'content-type':'application/json'}
    
    # Perform request
    r = requests.post('http://amp.pharm.mssm.edu/L1000CDS2/query',data=json.dumps(payload),headers=headers)
    resCD= r.json()
    
    # Add URL
    resCD['resultUrl'] = 'http://amp.pharm.mssm.edu/L1000CDS2/#/result/' + resCD['shareId']
    
    # Return result
    return resCD


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

