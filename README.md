# Standardised Minimum 3D Distance (SM3DD)

Standardised Minimum 3D Distance (SM3DD)
Source code and binary executable for "SM3DD with Segmented PCA: A Comprehensive Method for Interpreting 3D Spatial Transcriptomics".
 
# Installation and Compilation Instructions
Running the standalone binary executable requires MATLAB Runtime to be installed (https://au.mathworks.com/help/compiler/install-the-matlab-runtime.html). Currently the source code/binary executable is only the Windows version. Check the lab's GitHub for updates and additions (https://github.com/clinicalomx/Standardised-Minimum-3D-Distance-SM3DD).
 
File structure requires a folder named "SM3DD" to be in the same folder as the script/binary executable. The 4 CSV files need to be in the "SM3DD" folder, along with another folder named "Data". The probe coordinate files (..._tx_file(s) for NanoString assays) and at least one probe count matrix file (..._exprMat_file(s) for NanoString assays) should be put the "Data" folder (any folder structure within the "Data" folder is permissible).
 
Use Sample Groups.csv to: (i) list the samples in each group to be compared. (ii) name the sample groups (iii) set the paired test binary value (1 if paired test)
 
Use Samples to FOVs.csv to assign sample names to each Field of View on each slide. The script skips calculating SM3DDs for FOVs with sample names that are not included in Sample Groups.csv
 
Use Platform.csv to assign which assay the data comes from and define the data structure within the probe coordinate and probe count matrix files, along with the first 6 letters of any control probe type names. Additional lines can be added for custom assays etc. The "number of sections" options defines the segmentation level for the FROM-probe-to-TO-probe matrix. Particularly for whole transcriptome assays, the RAM requirements for the FROM-probe-to-TO-probe matrix become quite large. Segmenting reduces RAM requirements without noticably impacting performance.
 
Use PCA options.csv to indicate whether you want FOVs with missing data to be included in the segmented PCA assays. TO-probes are excluded if there is any missing data, so datasets with one or more sparse FOVs may suffer from their inclusion.
 
SM3DD will terminate with and generate an error message to the same folder as the script if the files are not as expected, or if Sample Groups.csv includes any samples that are not listed in Samples to FOVs.csv SM3DD will generate an error message without terminating if FOVs include no count data.
 
The script checks for what SM3DDs have already been calculated etc, so if additional comparisons are run that include FOVs for which the SM3DDs have already been calculated, then it doesn't recalculate them. Similarly, if the script is interrupted for any reason, then it will, when restarted, continue on from where it left off.
 
To aid in estimating completion time, the script generates and updates FOVCompletionTimes.csv, which lists the time taken to calculate the M3DDs for each FOV (in seconds).
 
The script outputs DirectionalLog10P values, which are positive if "The FROM and TO probes are closer in Sample Group 2". The script outputs a Composite PCA differences matrix (one row per FROM-probe, PCA dimensions in columns), ready for pathway level analysis.
 
Disk space required for the M3DD data can be considerable. As a guide, while our NanoString 1K assay, with 116 FOVs and ~100 millin transcripts in total only required ~362 Gb, a 48-FOV Nanostring WTA assay required ~18 Tb.
 
The M3DD calculations are parallelised and the PCA function is internally parallelised, so core availability will impact run time.
 
# Citation
SM3DD with Segmented PCA: A Comprehensive Method for Interpreting 3D Spatial Transcriptomics, NAR Genomics and Bioinformatics (2026), in press.
 
Installation includes fdr_bky.m David Groppe (2025). Two-stage Benjamini, Krieger, & Yekutieli FDR procedure (https://www.mathworks.com/matlabcentral/fileexchange/27423-two-stage-benjamini-krieger-yekutieli-fdr-procedure), MATLAB Central File Exchange. Retrieved June 20, 2025.
 
# License
Licensed under the Academic Free License version 3.0
