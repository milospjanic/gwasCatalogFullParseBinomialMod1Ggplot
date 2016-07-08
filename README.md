# gwasCatalogFullParseBinomialMod1Ggplot

This script will connect to http://www.genome.gov, download GWAS Catalog. It will download and add CardiogramPlusC4D GWAS SNPs to the GWASCatalog. Next, it will convert downloaded GWAS Catalog file to a bed format with columns chr;position;position+1;proxy_gene;phenotype, and then create separate bed files for each unique GWAS Catalog category using 5th column as identifier. All special caracters will be removed and substituted with underscore. It will then perform binomial statistics on your input bed file and such parsed collection of GWAS Catalog files using modified binomial as described in https://github.com/milospjanic/gwasanalytics. Outputs are table with binomial -log10pvalues and fold changes in file output.table.txt and graphical representation in file output.pdf. Screen output will consist of series of intermediary steps and additional parameters. Examples of outputs shown below.

#Textual output in output.table.txt
<pre>
GWAS Catalog Phenotype Total SNPs Overlap Fold change Fraction of hg19 Peak coverage
Parkinson's_disease 76 2 2.63158 1.31801e-05 41348
Fat_distribution_(HIV) 14 1 7.14286 6.53425e-06 20499
Breast_size 31 1 3.22581 6.51831e-06 20449
Periodontitis_(CDC-AAP) 26 1 3.84615 6.66176e-06 20899
Keloid 4 1 25 6.50238e-06 20399
CardiogramPlusC4D 52 1 1.92308 6.50238e-06 20399
Response_to_alcohol_consumption_(flushing_response) 3 1 33.3333 6.66176e-06 20899
Multiple_sclerosis_(severity) 12 1 8.33333 6.55019e-06 20549
Tanning 15 2 13.3333 1.32119e-05 41448
Systemic_lupus_erythematosus 108 3 2.77778 1.94912e-05 61147
Obesity-related_traits 833 5 0.60024 3.38507e-05 106195
Migraine 70 2 2.85714 1.31482e-05 41248
Free_thyroxine_concentration 3 1 33.3333 6.56613e-06 20599
Hyperactive-impulsive_symptoms 13 1 7.69231 6.51831e-06 20449
Crohn's_disease_(need_for_surgery) 2 1 50 7.42678e-06 23299
Systemic_lupus_erythematosus_and_Systemic_sclerosis 20 1 5 6.48644e-06 20349
Bulimia_nervosa 23 1 4.34783 6.66176e-06 20899
Iris_characteristics 5 1 20 6.50238e-06 20399
Stroke_(pediatric) 3 1 33.3333 6.51831e-06 20449
Cystic_fibrosis_severity 6 1 16.6667 6.53425e-06 20499
Multiple_sclerosis 162 1 0.617284 6.4705e-06 20299
Serum_dimethylarginine_levels_(symmetric) 30 1 3.33333 1.26704e-05 39749
Alcohol_dependence_(age_at_onset) 26 2 7.69231 1.30844e-05 41048
...
</pre>
#Graphical output in output.pdf
![Screenshot](https://github.com/milospjanic/gwasCatalogFullParseBinomialMod1Ggplot/blob/master/example.png)

#Screen output
<pre>
./gwasCatalogFullScanBinomialGgplot.sh ARNT.chipseq.cut.10000 
--2016-06-16 12:32:09--  http://www.genome.gov/admin/gwascatalog.txt
Resolving www.genome.gov (www.genome.gov)... 156.40.242.24
Connecting to www.genome.gov (www.genome.gov)|156.40.242.24|:80... connected.
HTTP request sent, awaiting response... 301 Moved Permanently
Location: https://www.genome.gov/admin/gwascatalog.txt [following]
--2016-06-16 12:32:09--  https://www.genome.gov/admin/gwascatalog.txt
Connecting to www.genome.gov (www.genome.gov)|156.40.242.24|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 10407265 (9.9M) [text/plain]
Saving to: `gwascatalog.txt'

100%[=======================================================================================================================================================>] 10,407,265  1.47M/s   in 7.1s    

2016-06-16 12:32:17 (1.40 MB/s) - `gwascatalog.txt' saved [10407265/10407265]

--2016-06-16 12:32:18--  https://stanfordmedicine.box.com/shared/static/pqxkuzwgv8bhl8ne05a3ohlmzgwbir28.bed
Resolving stanfordmedicine.box.com (stanfordmedicine.box.com)... 74.112.184.85, 107.152.24.197, 74.112.185.182
Connecting to stanfordmedicine.box.com (stanfordmedicine.box.com)|74.112.184.85|:443... connected.
HTTP request sent, awaiting response... 302 Found
Location: https://stanfordmedicine.app.box.com/shared/static/pqxkuzwgv8bhl8ne05a3ohlmzgwbir28.bed [following]
--2016-06-16 12:32:18--  https://stanfordmedicine.app.box.com/shared/static/pqxkuzwgv8bhl8ne05a3ohlmzgwbir28.bed
Resolving stanfordmedicine.app.box.com (stanfordmedicine.app.box.com)... 74.112.185.87, 74.112.184.87, 107.152.24.199
Connecting to stanfordmedicine.app.box.com (stanfordmedicine.app.box.com)|74.112.185.87|:443... connected.
HTTP request sent, awaiting response... 302 Found
Location: https://dl.boxcloud.com/d/1/68Jy8knXO2QKrBbbRXd88q2tjBgSZ-O1IHUO8OMH7MwUf963W573e-JFFoHJ7EAJDHUlZRmb-X5ltqHgFJGdPlmsPUdEAFW0e4hGculzAFt-rk7108ylZovCOs-FO8ii9nXylV6UIxzHX9-w2-qg37FMwEuDVyH_Sg8aqEWI7MX4GjTDg5nj9JRzNumSDRjwL1zf8nULuiLBEE5c_eeZReW2Pq3BnI6G44Tu2pD0l9tfXjF9fFTRaANo5phCkSU_2cAc2aGsXeTaI5zBiHjbclc2C2QHKsiSByzg13FXKBcnqbge2qzLRnjfyTQZjR6guZbX3x0tiMIYRtAs4xesRZ71R5PqFMiJ4qFM1c_ZEZJ82QqOEhgzT-pN6shMVfe9Quxm0b3rO_vr3gbOkpFD3qnNn8bkmznaGs2s2IDoR7qLHIjLU4mRpkV9_SFzEiOj4827LOtPNKgV0tSPiOu6PpudRQlfkYtnAP5WnHNBK9-mHkBQ79weScffqf8oNbsAPd_rqLAtyy9vcBZVpev08QOy_hhnzTCdyxVUrsPtgUc8Vowh_m1iK-AlKQFDdDn8JSNH8bqpEqHQcaNEp6UA-haJ_hd6hpKrftZSENiXaE7cNlwgoOalny2DHr4pVpuNTYPNSN6V-WJSana1mQSCv0A9FgW2TsnN979wHJpmA9XPnPs4MEKbOjWaK4kcAKuldOOAvdehX-942AUmm9w6banYWMu79DTumvtrzh6N-CHMHpsgkY0WGGMiDuhtuRrgtv4T5cMF8fHM79MFDXLMpdBy2gKCxKP3hD3ifoeoy5vAlSGoqzivsrfzXtPLpc5juvzZFhi_iTHJNeyVNBEkzwMEILaF_u4CiT0cAcBrbt917LpurTE0IrDlypOJql6YvgscokLxqI2K30N64nd7VxY5Spp1g8m1ZnCTbIyF9HbFQtjn5Hi6tO0-B559hCzHEz4mY_5zNif5B_2dtuzC-2W9mP9PvDQVwW3hkxyddMKd0Oo8xd3cRlf7Rt6dBbEQAMTcpizNbZr9LKaB6-OMrvdr8SarBw2CWk58hXQp3JY_iJSWG5Xln4uCo5x_1hhWJOv5MiSRm6dP_BRThrVcTsgwtZ7ErOHwGh7eR0cGCDAH2U5F9XWauaB0mp0fnbwk6xMnxXRQ5Mpesem6BL21bxHUFFuP-eFQIvolyb3pnvCr3Njxo8K8/download [following]
--2016-06-16 12:32:18--  https://dl.boxcloud.com/d/1/68Jy8knXO2QKrBbbRXd88q2tjBgSZ-O1IHUO8OMH7MwUf963W573e-JFFoHJ7EAJDHUlZRmb-X5ltqHgFJGdPlmsPUdEAFW0e4hGculzAFt-rk7108ylZovCOs-FO8ii9nXylV6UIxzHX9-w2-qg37FMwEuDVyH_Sg8aqEWI7MX4GjTDg5nj9JRzNumSDRjwL1zf8nULuiLBEE5c_eeZReW2Pq3BnI6G44Tu2pD0l9tfXjF9fFTRaANo5phCkSU_2cAc2aGsXeTaI5zBiHjbclc2C2QHKsiSByzg13FXKBcnqbge2qzLRnjfyTQZjR6guZbX3x0tiMIYRtAs4xesRZ71R5PqFMiJ4qFM1c_ZEZJ82QqOEhgzT-pN6shMVfe9Quxm0b3rO_vr3gbOkpFD3qnNn8bkmznaGs2s2IDoR7qLHIjLU4mRpkV9_SFzEiOj4827LOtPNKgV0tSPiOu6PpudRQlfkYtnAP5WnHNBK9-mHkBQ79weScffqf8oNbsAPd_rqLAtyy9vcBZVpev08QOy_hhnzTCdyxVUrsPtgUc8Vowh_m1iK-AlKQFDdDn8JSNH8bqpEqHQcaNEp6UA-haJ_hd6hpKrftZSENiXaE7cNlwgoOalny2DHr4pVpuNTYPNSN6V-WJSana1mQSCv0A9FgW2TsnN979wHJpmA9XPnPs4MEKbOjWaK4kcAKuldOOAvdehX-942AUmm9w6banYWMu79DTumvtrzh6N-CHMHpsgkY0WGGMiDuhtuRrgtv4T5cMF8fHM79MFDXLMpdBy2gKCxKP3hD3ifoeoy5vAlSGoqzivsrfzXtPLpc5juvzZFhi_iTHJNeyVNBEkzwMEILaF_u4CiT0cAcBrbt917LpurTE0IrDlypOJql6YvgscokLxqI2K30N64nd7VxY5Spp1g8m1ZnCTbIyF9HbFQtjn5Hi6tO0-B559hCzHEz4mY_5zNif5B_2dtuzC-2W9mP9PvDQVwW3hkxyddMKd0Oo8xd3cRlf7Rt6dBbEQAMTcpizNbZr9LKaB6-OMrvdr8SarBw2CWk58hXQp3JY_iJSWG5Xln4uCo5x_1hhWJOv5MiSRm6dP_BRThrVcTsgwtZ7ErOHwGh7eR0cGCDAH2U5F9XWauaB0mp0fnbwk6xMnxXRQ5Mpesem6BL21bxHUFFuP-eFQIvolyb3pnvCr3Njxo8K8/download
Resolving dl.boxcloud.com (dl.boxcloud.com)... 107.152.24.200, 74.112.185.96, 74.112.184.96
Connecting to dl.boxcloud.com (dl.boxcloud.com)|107.152.24.200|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 1784 (1.7K) [application/octet-stream]
Saving to: `pqxkuzwgv8bhl8ne05a3ohlmzgwbir28.bed'

100%[=======================================================================================================================================================>] 1,784       --.-K/s   in 0s      

2016-06-16 12:32:18 (46.3 MB/s) - `pqxkuzwgv8bhl8ne05a3ohlmzgwbir28.bed' saved [1784/1784]

Gwas Catalog number of SNP-phenotype associations:
18950 GwasCatalog.bed
***Finished parsing Gwas Catalog SNP-phenotype associations per category***
***Converting to GWAS Phenotype specific bed files***
***Finished creating Gwas Catalog phenotype specific bed files***
***Overlapping Phenotype SNPs with input bed***
***Finished overlapping GWAS Catalog Phenotype specific bed with input bed***
***Downloading hg19 chromosome sizes***
--2016-06-16 12:32:55--  https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
Resolving genome.ucsc.edu (genome.ucsc.edu)... 128.114.119.133, 128.114.119.136, 128.114.119.132, ...
Connecting to genome.ucsc.edu (genome.ucsc.edu)|128.114.119.133|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 1971 (1.9K) [text/plain]
Saving to: `hg19.chrom.sizes'

100%[=======================================================================================================================================================>] 1,971       --.-K/s   in 0s      

2016-06-16 12:32:55 (92.6 MB/s) - `hg19.chrom.sizes' saved [1971/1971]

Human Genome size version hg19: 3137161264
1249 ARNT.chipseq.cut.10000
***Selecting FDR 10% overlaps***
   9 ./GWASCatalogPhenotype_Height.txt.cut.sort.uniq.chrXY.overlap
   6 ./GWASCatalogPhenotype_IgG_glycosylation.txt.cut.sort.uniq.chrXY.overlap
   5 ./GWASCatalogPhenotype_Obesity-related_traits.txt.cut.sort.uniq.chrXY.overlap
   5 ./GWASCatalogPhenotype_Crohn's_disease.txt.cut.sort.uniq.chrXY.overlap
   4 ./GWASCatalogPhenotype_Ulcerative_colitis.txt.cut.sort.uniq.chrXY.overlap
***Start calculating coverage***
***Overlapping Phenotype SNPs with input bed to calculate coverage***
***Creating R script***
GWAS Catalog Phenotype: Parkinson's_disease Total SNPs: 76 Overlap: 2 Fold change: 2.63158 Fraction of hg19 1.31801e-05 Peak coverage: 41348
GWAS Catalog Phenotype: Fat_distribution_(HIV) Total SNPs: 14 Overlap: 1 Fold change: 7.14286 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Breast_size Total SNPs: 31 Overlap: 1 Fold change: 3.22581 Fraction of hg19 6.51831e-06 Peak coverage: 20449
GWAS Catalog Phenotype: Periodontitis_(CDC-AAP) Total SNPs: 26 Overlap: 1 Fold change: 3.84615 Fraction of hg19 6.66176e-06 Peak coverage: 20899
GWAS Catalog Phenotype: Keloid Total SNPs: 4 Overlap: 1 Fold change: 25 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: CardiogramPlusC4D Total SNPs: 52 Overlap: 1 Fold change: 1.92308 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: Response_to_alcohol_consumption_(flushing_response) Total SNPs: 3 Overlap: 1 Fold change: 33.3333 Fraction of hg19 6.66176e-06 Peak coverage: 20899
GWAS Catalog Phenotype: Multiple_sclerosis_(severity) Total SNPs: 12 Overlap: 1 Fold change: 8.33333 Fraction of hg19 6.55019e-06 Peak coverage: 20549
GWAS Catalog Phenotype: Tanning Total SNPs: 15 Overlap: 2 Fold change: 13.3333 Fraction of hg19 1.32119e-05 Peak coverage: 41448
GWAS Catalog Phenotype: Systemic_lupus_erythematosus Total SNPs: 108 Overlap: 3 Fold change: 2.77778 Fraction of hg19 1.94912e-05 Peak coverage: 61147
GWAS Catalog Phenotype: Obesity-related_traits Total SNPs: 833 Overlap: 5 Fold change: 0.60024 Fraction of hg19 3.38507e-05 Peak coverage: 106195
GWAS Catalog Phenotype: Migraine Total SNPs: 70 Overlap: 2 Fold change: 2.85714 Fraction of hg19 1.31482e-05 Peak coverage: 41248
GWAS Catalog Phenotype: Free_thyroxine_concentration Total SNPs: 3 Overlap: 1 Fold change: 33.3333 Fraction of hg19 6.56613e-06 Peak coverage: 20599
GWAS Catalog Phenotype: Hyperactive-impulsive_symptoms Total SNPs: 13 Overlap: 1 Fold change: 7.69231 Fraction of hg19 6.51831e-06 Peak coverage: 20449
GWAS Catalog Phenotype: Crohn's_disease_(need_for_surgery) Total SNPs: 2 Overlap: 1 Fold change: 50 Fraction of hg19 7.42678e-06 Peak coverage: 23299
GWAS Catalog Phenotype: Systemic_lupus_erythematosus_and_Systemic_sclerosis Total SNPs: 20 Overlap: 1 Fold change: 5 Fraction of hg19 6.48644e-06 Peak coverage: 20349
GWAS Catalog Phenotype: Bulimia_nervosa Total SNPs: 23 Overlap: 1 Fold change: 4.34783 Fraction of hg19 6.66176e-06 Peak coverage: 20899
GWAS Catalog Phenotype: Iris_characteristics Total SNPs: 5 Overlap: 1 Fold change: 20 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: Stroke_(pediatric) Total SNPs: 3 Overlap: 1 Fold change: 33.3333 Fraction of hg19 6.51831e-06 Peak coverage: 20449
GWAS Catalog Phenotype: Cystic_fibrosis_severity Total SNPs: 6 Overlap: 1 Fold change: 16.6667 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Multiple_sclerosis Total SNPs: 162 Overlap: 1 Fold change: 0.617284 Fraction of hg19 6.4705e-06 Peak coverage: 20299
GWAS Catalog Phenotype: Serum_dimethylarginine_levels_(symmetric) Total SNPs: 30 Overlap: 1 Fold change: 3.33333 Fraction of hg19 1.26704e-05 Peak coverage: 39749
GWAS Catalog Phenotype: Alcohol_dependence_(age_at_onset) Total SNPs: 26 Overlap: 2 Fold change: 7.69231 Fraction of hg19 1.30844e-05 Peak coverage: 41048
GWAS Catalog Phenotype: Upper_aerodigestive_tract_cancers Total SNPs: 5 Overlap: 1 Fold change: 20 Fraction of hg19 6.70957e-06 Peak coverage: 21049
GWAS Catalog Phenotype: Narcolepsy_with_cataplexy Total SNPs: 4 Overlap: 1 Fold change: 25 Fraction of hg19 6.55019e-06 Peak coverage: 20549
GWAS Catalog Phenotype: Type_2_diabetes Total SNPs: 192 Overlap: 2 Fold change: 1.04167 Fraction of hg19 1.30048e-05 Peak coverage: 40798
GWAS Catalog Phenotype: Response_to_tamoxifen_in_breast_cancer Total SNPs: 1 Overlap: 1 Fold change: 100 Fraction of hg19 6.55019e-06 Peak coverage: 20549
GWAS Catalog Phenotype: Large_artery_stroke Total SNPs: 9 Overlap: 1 Fold change: 11.1111 Fraction of hg19 6.62988e-06 Peak coverage: 20799
GWAS Catalog Phenotype: Amyotrophic_lateral_sclerosis_(sporadic) Total SNPs: 95 Overlap: 2 Fold change: 2.10526 Fraction of hg19 1.83918e-05 Peak coverage: 57698
GWAS Catalog Phenotype: Crohn's_disease Total SNPs: 180 Overlap: 5 Fold change: 2.77778 Fraction of hg19 3.8122e-05 Peak coverage: 119595
GWAS Catalog Phenotype: Plasma_plasminogen_activator_levels Total SNPs: 5 Overlap: 1 Fold change: 20 Fraction of hg19 6.598e-06 Peak coverage: 20699
GWAS Catalog Phenotype: Post-traumatic_stress_disorder_(asjusted_for_relatedness) Total SNPs: 15 Overlap: 1 Fold change: 6.66667 Fraction of hg19 6.48644e-06 Peak coverage: 20349
GWAS Catalog Phenotype: Nicotine_dependence Total SNPs: 4 Overlap: 1 Fold change: 25 Fraction of hg19 6.58207e-06 Peak coverage: 20649
GWAS Catalog Phenotype: Response_to_mTOR_inhibitor_(rapamycin) Total SNPs: 5 Overlap: 1 Fold change: 20 Fraction of hg19 6.48644e-06 Peak coverage: 20349
GWAS Catalog Phenotype: Dialysis-related_mortality Total SNPs: 26 Overlap: 1 Fold change: 3.84615 Fraction of hg19 1.18257e-05 Peak coverage: 37099
GWAS Catalog Phenotype: Pulmonary_function_(interaction) Total SNPs: 43 Overlap: 1 Fold change: 2.32558 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Refractive_error Total SNPs: 29 Overlap: 1 Fold change: 3.44828 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: IgG_glycosylation Total SNPs: 395 Overlap: 6 Fold change: 1.51899 Fraction of hg19 3.91896e-05 Peak coverage: 122944
GWAS Catalog Phenotype: Pulmonary_function_decline Total SNPs: 39 Overlap: 1 Fold change: 2.5641 Fraction of hg19 6.48644e-06 Peak coverage: 20349
GWAS Catalog Phenotype: Testicular_germ_cell_tumor Total SNPs: 21 Overlap: 1 Fold change: 4.7619 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Mean_platelet_volume Total SNPs: 53 Overlap: 1 Fold change: 1.88679 Fraction of hg19 6.55019e-06 Peak coverage: 20549
GWAS Catalog Phenotype: Cannabis_dependence Total SNPs: 10 Overlap: 1 Fold change: 10 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Non-melanoma_skin_cancer Total SNPs: 4 Overlap: 1 Fold change: 25 Fraction of hg19 6.64582e-06 Peak coverage: 20849
GWAS Catalog Phenotype: Airflow_obstruction_ Total SNPs: 30 Overlap: 2 Fold change: 6.66667 Fraction of hg19 1.31482e-05 Peak coverage: 41248
GWAS Catalog Phenotype: Pancreatic_cancer Total SNPs: 39 Overlap: 1 Fold change: 2.5641 Fraction of hg19 6.48644e-06 Peak coverage: 20349
GWAS Catalog Phenotype: Breast_cancer Total SNPs: 114 Overlap: 1 Fold change: 0.877193 Fraction of hg19 6.72551e-06 Peak coverage: 21099
GWAS Catalog Phenotype: Interstitial_lung_disease_ Total SNPs: 18 Overlap: 1 Fold change: 5.55556 Fraction of hg19 6.62988e-06 Peak coverage: 20799
GWAS Catalog Phenotype: Blood_metabolite_levels Total SNPs: 183 Overlap: 3 Fold change: 1.63934 Fraction of hg19 2.55636e-05 Peak coverage: 80197
GWAS Catalog Phenotype: Bone_mineral_density Total SNPs: 112 Overlap: 1 Fold change: 0.892857 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Response_to_serotonin_reuptake_inhibitors_in_major_depressive_disorder_(plasma_drug_and_metabolite_levels) Total SNPs: 12 Overlap: 1 Fold change: 8.33333 Fraction of hg19 1.65433e-05 Peak coverage: 51899
GWAS Catalog Phenotype: Platelet_counts Total SNPs: 76 Overlap: 1 Fold change: 1.31579 Fraction of hg19 6.64582e-06 Peak coverage: 20849
GWAS Catalog Phenotype: Type_1_diabetes Total SNPs: 78 Overlap: 2 Fold change: 2.5641 Fraction of hg19 1.84715e-05 Peak coverage: 57948
GWAS Catalog Phenotype: Heart_rate Total SNPs: 37 Overlap: 1 Fold change: 2.7027 Fraction of hg19 1.1491e-05 Peak coverage: 36049
GWAS Catalog Phenotype: Coronary_heart_disease Total SNPs: 130 Overlap: 4 Fold change: 3.07692 Fraction of hg19 2.60414e-05 Peak coverage: 81696
GWAS Catalog Phenotype: Metabolite_levels__(X-11787) Total SNPs: 29 Overlap: 1 Fold change: 3.44828 Fraction of hg19 6.85301e-06 Peak coverage: 21499
GWAS Catalog Phenotype: Response_to_mTOR_inhibitor_(everolimus) Total SNPs: 8 Overlap: 1 Fold change: 12.5 Fraction of hg19 9.08433e-06 Peak coverage: 28499
GWAS Catalog Phenotype: Response_to_taxane_treatment_(placlitaxel) Total SNPs: 13 Overlap: 1 Fold change: 7.69231 Fraction of hg19 1.26704e-05 Peak coverage: 39749
GWAS Catalog Phenotype: Hematological_parameters Total SNPs: 11 Overlap: 1 Fold change: 9.09091 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Neutrophil_count Total SNPs: 22 Overlap: 1 Fold change: 4.54545 Fraction of hg19 9.16402e-06 Peak coverage: 28749
GWAS Catalog Phenotype: Intracranial_aneurysm Total SNPs: 13 Overlap: 3 Fold change: 23.0769 Fraction of hg19 1.95071e-05 Peak coverage: 61197
GWAS Catalog Phenotype: Non-obstructive_azoospermia Total SNPs: 5 Overlap: 1 Fold change: 20 Fraction of hg19 6.48644e-06 Peak coverage: 20349
GWAS Catalog Phenotype: Response_to_antipsychotic_therapy_(extrapyramidal_side_effects) Total SNPs: 18 Overlap: 1 Fold change: 5.55556 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: Dental_caries Total SNPs: 63 Overlap: 1 Fold change: 1.5873 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Stroke_(ischemic) Total SNPs: 8 Overlap: 1 Fold change: 12.5 Fraction of hg19 6.62988e-06 Peak coverage: 20799
GWAS Catalog Phenotype: Response_to_anti-TNF_treatment_in_rheumatoid_arthritis Total SNPs: 4 Overlap: 1 Fold change: 25 Fraction of hg19 6.51831e-06 Peak coverage: 20449
GWAS Catalog Phenotype: Leprosy Total SNPs: 9 Overlap: 1 Fold change: 11.1111 Fraction of hg19 6.56613e-06 Peak coverage: 20599
GWAS Catalog Phenotype: Height Total SNPs: 460 Overlap: 9 Fold change: 1.95652 Fraction of hg19 6.26493e-05 Peak coverage: 196541
GWAS Catalog Phenotype: Ulcerative_colitis Total SNPs: 118 Overlap: 4 Fold change: 3.38983 Fraction of hg19 2.62645e-05 Peak coverage: 82396
GWAS Catalog Phenotype: Sudden_cardiac_arrest Total SNPs: 50 Overlap: 2 Fold change: 4 Fraction of hg19 1.31641e-05 Peak coverage: 41298
GWAS Catalog Phenotype: Prostate_cancer Total SNPs: 107 Overlap: 4 Fold change: 3.73832 Fraction of hg19 2.68383e-05 Peak coverage: 84196
GWAS Catalog Phenotype: Response_to_anti-retroviral_therapy_(ddI-d4T)_in_HIV-1_infection_(Grade_3_peripheral_neuropathy) Total SNPs: 16 Overlap: 2 Fold change: 12.5 Fraction of hg19 1.30685e-05 Peak coverage: 40998
GWAS Catalog Phenotype: Endometriosis Total SNPs: 25 Overlap: 1 Fold change: 4 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: Bilirubin_levels Total SNPs: 38 Overlap: 1 Fold change: 2.63158 Fraction of hg19 6.56613e-06 Peak coverage: 20599
GWAS Catalog Phenotype: Schizophrenia Total SNPs: 138 Overlap: 1 Fold change: 0.724638 Fraction of hg19 6.64582e-06 Peak coverage: 20849
GWAS Catalog Phenotype: Migraine_with_aura Total SNPs: 23 Overlap: 1 Fold change: 4.34783 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: Butyrylcholinesterase_levels Total SNPs: 3 Overlap: 1 Fold change: 33.3333 Fraction of hg19 6.66176e-06 Peak coverage: 20899
GWAS Catalog Phenotype: Self-reported_allergy Total SNPs: 35 Overlap: 1 Fold change: 2.85714 Fraction of hg19 6.61394e-06 Peak coverage: 20749
GWAS Catalog Phenotype: Myocardial_infarction_(early_onset) Total SNPs: 9 Overlap: 1 Fold change: 11.1111 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: Chronic_lymphocytic_leukemia Total SNPs: 38 Overlap: 1 Fold change: 2.63158 Fraction of hg19 6.64582e-06 Peak coverage: 20849
GWAS Catalog Phenotype: Freckles Total SNPs: 3 Overlap: 1 Fold change: 33.3333 Fraction of hg19 6.64582e-06 Peak coverage: 20849
GWAS Catalog Phenotype: Coronary_artery_calcification Total SNPs: 55 Overlap: 1 Fold change: 1.81818 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: QT_interval Total SNPs: 76 Overlap: 2 Fold change: 2.63158 Fraction of hg19 2.49104e-05 Peak coverage: 78148
GWAS Catalog Phenotype: Myasthenia_gravis_ Total SNPs: 5 Overlap: 1 Fold change: 20 Fraction of hg19 6.48644e-06 Peak coverage: 20349
GWAS Catalog Phenotype: HIV-1_susceptibility Total SNPs: 16 Overlap: 1 Fold change: 6.25 Fraction of hg19 6.58207e-06 Peak coverage: 20649
GWAS Catalog Phenotype: Venous_thromboembolism_(gene_x_gene_interaction) Total SNPs: 37 Overlap: 1 Fold change: 2.7027 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Cognitive_performance Total SNPs: 118 Overlap: 1 Fold change: 0.847458 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Plasminogen_activator_inhibitor_type_1_levels_(PAI-1) Total SNPs: 4 Overlap: 1 Fold change: 25 Fraction of hg19 6.598e-06 Peak coverage: 20699
GWAS Catalog Phenotype: Bipolar_disorder Total SNPs: 117 Overlap: 2 Fold change: 1.7094 Fraction of hg19 1.3196e-05 Peak coverage: 41398
GWAS Catalog Phenotype: Inflammatory_bowel_disease Total SNPs: 117 Overlap: 3 Fold change: 2.5641 Fraction of hg19 1.97781e-05 Peak coverage: 62047
GWAS Catalog Phenotype: Metabolic_traits Total SNPs: 48 Overlap: 1 Fold change: 2.08333 Fraction of hg19 6.53425e-06 Peak coverage: 20499
GWAS Catalog Phenotype: Menopause_(age_at_onset) Total SNPs: 36 Overlap: 2 Fold change: 5.55556 Fraction of hg19 1.31641e-05 Peak coverage: 41298
GWAS Catalog Phenotype: Smooth-surface_caries_ Total SNPs: 42 Overlap: 1 Fold change: 2.38095 Fraction of hg19 6.51831e-06 Peak coverage: 20449
GWAS Catalog Phenotype: Immune_response_to_smallpox_vaccine_(IL-6) Total SNPs: 46 Overlap: 1 Fold change: 2.17391 Fraction of hg19 6.50238e-06 Peak coverage: 20399
GWAS Catalog Phenotype: Primary_biliary_cirrhosis Total SNPs: 31 Overlap: 1 Fold change: 3.22581 Fraction of hg19 6.56613e-06 Peak coverage: 20599
GWAS Catalog Phenotype: Serum_IgA_levels Total SNPs: 1 Overlap: 1 Fold change: 100 Fraction of hg19 6.56613e-06 Peak coverage: 20599
GWAS Catalog Phenotype: Prostate_cancer_(gene_x_gene_interaction) Total SNPs: 51 Overlap: 1 Fold change: 1.96078 Fraction of hg19 1.17141e-05 Peak coverage: 36749
GWAS Catalog Phenotype: Insulin_resistance-response Total SNPs: 12 Overlap: 1 Fold change: 8.33333 Fraction of hg19 6.64582e-06 Peak coverage: 20849
GWAS Catalog Phenotype: Body_mass_index Total SNPs: 151 Overlap: 1 Fold change: 0.662252 Fraction of hg19 6.61394e-06 Peak coverage: 20749
GWAS Catalog Phenotype: Bladder_cancer_(smoking_interaction) Total SNPs: 5 Overlap: 1 Fold change: 20 Fraction of hg19 6.61394e-06 Peak coverage: 20749
GWAS Catalog Phenotype: Blood_pressure Total SNPs: 90 Overlap: 1 Fold change: 1.11111 Fraction of hg19 6.66176e-06 Peak coverage: 20899
Finishing R script
Writing ggplot2 code for plot in output.pdf
[1] "***Final table - top 25 categories - for plot - output.pdf***"
                                                                                                     LogP
Height                                                                                           44.82941
IgG_glycosylation                                                                                31.64193
Crohn's_disease                                                                                  29.75900
Intracranial_aneurysm                                                                            26.87840
Prostate_cancer                                                                                  26.64892
Ulcerative_colitis                                                                               26.33884
Coronary_heart_disease                                                                           25.98106
Obesity-related_traits                                                                           22.67013
Systemic_lupus_erythematosus                                                                     20.31205
Inflammatory_bowel_disease                                                                       20.02612
Blood_metabolite_levels                                                                          17.90739
Tanning                                                                                          17.81500
Response_to_anti-retroviral_therapy_(ddI-d4T)_in_HIV-1_infection_(Grade_3_peripheral_neuropathy) 17.70330
Alcohol_dependence_(age_at_onset)                                                                16.70467
Airflow_obstruction_                                                                             16.40347
Menopause_(age_at_onset)                                                                         16.03076
Sudden_cardiac_arrest                                                                            15.36597
Migraine                                                                                         14.68989
Parkinson's_disease                                                                              14.51951
Type_1_diabetes                                                                                  13.79260
Bipolar_disorder                                                                                 13.65009
Amyotrophic_lateral_sclerosis_(sporadic)                                                         13.40490
QT_interval                                                                                      13.24722
Type_2_diabetes                                                                                  12.68623
Response_to_tamoxifen_in_breast_cancer                                                           11.93602
                                                                                                        FC
Height                                                                                             1.95652
IgG_glycosylation                                                                                  1.51899
Crohn's_disease                                                                                    2.77778
Intracranial_aneurysm                                                                             23.07690
Prostate_cancer                                                                                    3.73832
Ulcerative_colitis                                                                                 3.38983
Coronary_heart_disease                                                                             3.07692
Obesity-related_traits                                                                             0.60024
Systemic_lupus_erythematosus                                                                       2.77778
Inflammatory_bowel_disease                                                                         2.56410
Blood_metabolite_levels                                                                            1.63934
Tanning                                                                                           13.33330
Response_to_anti-retroviral_therapy_(ddI-d4T)_in_HIV-1_infection_(Grade_3_peripheral_neuropathy)  12.50000
Alcohol_dependence_(age_at_onset)                                                                  7.69231
Airflow_obstruction_                                                                               6.66667
Menopause_(age_at_onset)                                                                           5.55556
Sudden_cardiac_arrest                                                                              4.00000
Migraine                                                                                           2.85714
Parkinson's_disease                                                                                2.63158
Type_1_diabetes                                                                                    2.56410
Bipolar_disorder                                                                                   1.70940
Amyotrophic_lateral_sclerosis_(sporadic)                                                           2.10526
QT_interval                                                                                        2.63158
Type_2_diabetes                                                                                    1.04167
Response_to_tamoxifen_in_breast_cancer                                                           100.00000
                                                                                                 Phenotype SNPs
Height                                                                                                      460
IgG_glycosylation                                                                                           395
Crohn's_disease                                                                                             180
Intracranial_aneurysm                                                                                        13
Prostate_cancer                                                                                             107
Ulcerative_colitis                                                                                          118
Coronary_heart_disease                                                                                      130
Obesity-related_traits                                                                                      833
Systemic_lupus_erythematosus                                                                                108
Inflammatory_bowel_disease                                                                                  117
Blood_metabolite_levels                                                                                     183
Tanning                                                                                                      15
Response_to_anti-retroviral_therapy_(ddI-d4T)_in_HIV-1_infection_(Grade_3_peripheral_neuropathy)             16
Alcohol_dependence_(age_at_onset)                                                                            26
Airflow_obstruction_                                                                                         30
Menopause_(age_at_onset)                                                                                     36
Sudden_cardiac_arrest                                                                                        50
Migraine                                                                                                     70
Parkinson's_disease                                                                                          76
Type_1_diabetes                                                                                              78
Bipolar_disorder                                                                                            117
Amyotrophic_lateral_sclerosis_(sporadic)                                                                     95
QT_interval                                                                                                  76
Type_2_diabetes                                                                                             192
Response_to_tamoxifen_in_breast_cancer                                                                        1
                                                                                                             Category
Height                                                                                                          Other
IgG_glycosylation                                                                                               Other
Crohn's_disease                                                                                                 Other
Intracranial_aneurysm                                                                                           Other
Prostate_cancer                                                                                                Cancer
Ulcerative_colitis                                                                                              Other
Coronary_heart_disease                                                                                 Cardiovascular
Obesity-related_traits                                                                                          Other
Systemic_lupus_erythematosus                                                                     Chronic Inflammatory
Inflammatory_bowel_disease                                                                                      Other
Blood_metabolite_levels                                                                                         Other
Tanning                                                                                                         Other
Response_to_anti-retroviral_therapy_(ddI-d4T)_in_HIV-1_infection_(Grade_3_peripheral_neuropathy)                Other
Alcohol_dependence_(age_at_onset)                                                                               Other
Airflow_obstruction_                                                                                            Other
Menopause_(age_at_onset)                                                                                        Other
Sudden_cardiac_arrest                                                                                           Other
Migraine                                                                                                        Other
Parkinson's_disease                                                                                             Other
Type_1_diabetes                                                                                                 Other
Bipolar_disorder                                                                                                Other
Amyotrophic_lateral_sclerosis_(sporadic)                                                                        Other
QT_interval                                                                                                     Other
Type_2_diabetes                                                                                                 Other
Response_to_tamoxifen_in_breast_cancer                                                                         Cancer
Loading required package: methods
Loading required package: grid
Loading required package: quadprog
Loading required package: proto
Scale for 'colour' is already present. Adding another scale for 'colour', which will replace the existing scale.
null device 
          1 

</pre>
