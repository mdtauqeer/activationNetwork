# activationNetwork
Below files are required in the analysis of activation network.

Files 1-4 are huge files and they are uploaded to Mendeley Data repository. Download those files before analysis from https://data.mendeley.com/preview/6w3bfjp4x6?a=2841696b-6823-4e96-849a-b9993fe4ebf7

(1) gene_expression_all_genome2.Rdata
We have retrived microarray data (genome2 chip) from arrayexpress database. This file contains the processed expression data.
(2) gene_expression_all_avgConditionWise_genome2.Rdata
In this file, if samples from different studies belong to same condition then we make average. This is to avoid any bias for highly studied conditions. 
(3) gene_co_expression_avgConditionWise_correlation.Rdata
Gene expression were calculated using gene_expression_all_avgConditionWise data. We use correlation test.
(4) #metabolite_KO_output_GEM9.txt
flux profile for in silico metabolite knockout.


(5)sgadata_costanzo2009_stringentCutoff_101120.txt
This file is taken from costanzo et al 2009. This is genetic interactions with stringent cutoff threshold.
Download this file from https://boonelab.ccbr.utoronto.ca/supplement/costanzo2009/
(6) sce00001.txt
Metabolic map of Saccharomyces cerevisiae from KEGG database. download from here: https://www.genome.jp/brite/sce00001+YDR451C

(7) KO_Input_GEM9_Enzyme.txt
This is enzyme knockout input file for simulating enzyme KO. Which reactions are to be deleted by an enzyme KO are listed here.
(8) Enzyme_KO_output_GEM9.txt
flux profile for in silico enzyme knockout.
(9) KO_Input_GEM9_Metabolite.txt
This is metabolite knockout input file for simulating metabolite KO. Which reactions are to be deleted by a metabolite KO are listed here.
(10) KO_Input_GEM9_Gene.txt
This is gene knockout input file for simulating gene KO. Which reactions are to be deleted by a gene KO are listed here.
(11) Gene_KO_output_GEM9.txt
flux profile for in silico gene knockout.
(12) model_met.csv
List of metabolite and their details from activation network and Yeast9 metabolic model.
(13) activationNetwork.csv
Activation network constructed in this study.
(14) ec_details.csv
List of enzymes and their details from activation network and Yeast9 metabolic model.
(15) act_deg.csv
degree of activators in the activation network
(16) ec_class.csv
All the metabolic enzymes with their class.
(17) actNetcode_MSB.R
script for the analysis.
