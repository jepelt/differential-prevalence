--- Data processing ---
Raw data was processed with data_tses_16s_030625.R and data_tses_sg_030625.R scripts. The raw data files for 16S data can be downloaded from here https://zenodo.org/records/1146764

The processed data are in the files data_tses_16s_030625.rds and data_tses_sg_030625.rds as lists of TreeSummarizedExperiment objects.

Further data processing was done with the data_090625.R script. (data_090625.R was for creating an additional file and data_10...R were used to split the large data object to smaller files.)


--- Running BDPE (working name blogrrall) ---
The stan code for BDPE is in the file stan_blogrrall_120625.stan. The script run_blogrrall_120625.R was used to run BDPE in the Computers of CSC – IT Center for Science Ltd.  The results from running BDPE are in the file res_blogrrall_120625.rds.


--- Running the other methods ---
Scripts run_[METHOD_NAME]….R LDM was run in CSC and reults for LDM are in the file res_ldm_260625.rds.


--- Summarizing results and creating figures ---
The results and the figures 3-5 in the Results section of the paper can be computed/created with the scripts fig_[FIGURE_NUMBER]….R  


