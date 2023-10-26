# surveillance_simulation_metagenomics

The repository contains all data and codes used to generate the results presented in the publication “Modelling the effectiveness of surveillance based on metagenomics in detecting, monitoring, and forecasting antimicrobial resistance in livestock production under economic constraints” published in Scientific reports.

The data is available as excel spreadsheets.
The codes are available as R scripts. In the study, we used R version 4.3.1 and R studio version 1.4.1717.

The EXCEL spreadsheet “reads_seqdepth_vetII_samples” contains count data from metagenomic analysis of pig faeces that is used to define the probability distribution describing the concentration of antimicrobial resistance genes in pigs.

The R script “libraries” loads R packages used in the simulation, plotting and analyses.

The R codes “fig_1_6_Tet_efflux_samp5_20_50_100” and “fig_2_6_Aph_samp5_20_50_100” are simulating surveillance data and plotting the results equivalent to fig 1, 2 and 3 presented in the publication.

The R codes “plot_supl_fig_1_2_tet” and “plot_supl_fig_1_2_aph” are simulating surveillance data and plotting the results equivalent to fig 1, 2 and 3 presented in the publication.

The R scripts “Tete_dbest_samp5”, …., “Tete_dbest_samp100” and “aph_dbest_samp5”, …,  “aph_dbest_samp100” are simulating surveillance data and analyse the simulated data using the DBEST package for tetracycline resistance and gentamicin resistance respectively. Because we are using 1,001 iterations in each scenario for each number of samples per pool (5, 20, 50 and 100) the running time for these scripts are several hours on a standard laptop, the line defining the number of iterations is set to 2 iterations (n_series<-2). Below that line, there is an inactivated line (#n_series <- 1001). The output from the DBEST analysis using 1.001 iterations are available in the EXCEL spreadsheets “EStart5_Tetra” …. “EStart100_Tetra” and “ESTARTAPH5” … “ESTARTAPH100”.

Equivalent simulation for surveillance based on phenotype characterisation of resistance using 20 samples per month is performed in the scripts “Ecoli_tet_dbest_20” and “Efaecium_gent_dbest_20” and the output from the DBEST analysis is available in the spreadsheets  “E_coli_feno_df_detection”, “E_coli_feno_df_negative”, “E_facium_feno_df_detection” and “E_facium_feno_df_negative”.

The R scripts “fig_3_tet_ppv_npv” and “fig_4_aph_ppv_npv” are generating the plots presented in figures 3 and 4 in the article.

The EXCEL spreadsheet “SIR_n_farms” is utilising the estimated SIR model presented in the paper, to generate the number of infected farms per month.

The simulation of spread of an emerging ARG and surveillance data are done in the R scrips “emerge_BLATEM_5” … “emerge_BLATEM_100”, and results of the DBEST analysis is stored in the EXCEL spreadsheets “blatem_time_to_detect_5s” … “blatem_time_to_detect_100s”. The results of the data analysis from the simulated surveillance of an emerging ARG is plotted in the R script “repository_fig_5_emerging”.
