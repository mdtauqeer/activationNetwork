# activationNetwork
The following files are required for the analysis of the activation network:

################ download files from mendeley

Files 1–4 are large and have been uploaded to the Mendeley Data repository. Please download them before starting the analysis from: https://data.mendeley.com/datasets/6w3bfjp4x6/1

(1) gene_expression_all_genome2.Rdata

This file contains processed microarray expression data (genome2 chip) retrieved from the ArrayExpress database.

(2) gene_expression_all_avgConditionWise_genome2.Rdata

In this file, samples from different studies belonging to the same condition are averaged to minimize bias from highly studied conditions.

(3) gene_co_expression_avgConditionWise_correlation.Rdata

Gene co-expression values were calculated using the gene_expression_all_avgConditionWise data through a correlation test.

(4) metabolite_KO_output_GEM9.txt

This file contains the flux profile for an in silico metabolite knockout.

################ download files form other sources

(5)sgadata_costanzo2009_stringentCutoff_101120.txt

This file is sourced from Costanzo et al. (2009) and contains genetic interactions with a stringent cutoff threshold. You can download it from: https://boonelab.ccbr.utoronto.ca/supplement/costanzo2009/

(6) sce00001.txt

This file contains the metabolic map of Saccharomyces cerevisiae from the KEGG database. You can download it from: https://www.genome.jp/brite/sce00001+YDR451C

################  these files are available from GitHub

(7) KO_Input_GEM9_Enzyme.txt

This file lists the enzyme knockout input for simulating enzyme knockouts, specifying which reactions should be deleted by the enzyme knockout.

(8) Enzyme_KO_output_GEM9.txt

This file contains the flux profile for the in silico enzyme knockout.

(9) KO_Input_GEM9_Metabolite.txt

This file lists the metabolite knockout input for simulating metabolite knockouts, specifying which reactions should be deleted by the metabolite knockout.

(10) KO_Input_GEM9_Gene.txt

This file lists the gene knockout input for simulating gene knockouts, specifying which reactions should be deleted by the gene knockout.

(11) Gene_KO_output_GEM9.txt

This file contains the flux profile for the in silico gene knockout.

(12) model_met.csv

This file contains a list of metabolites and their details from the activation network and Yeast9 metabolic model.

(13) activationNetwork.csv

This file contains the activation network constructed in this study.

(14) ec_details.csv

This file contains a list of enzymes and their details from the activation network and Yeast9 metabolic model.

(15) act_deg.csv

This file contains the degree of activators in the activation network.

(16) ec_class.csv

This file lists all metabolic enzymes along with their classes.

(17) actNetcode_MSB.R

This is the script used for the analysis.

