# TRAP-RNAscope-Analysis
Code and data related to Figures 4 and S5 in Ryan et al., 2024

Fluorescent in situ hybridization data (RNAscope) was quantified using Qupath and exported as csv files.

RNAscopeAnalysis_D1_D2 is code used to analyze csv files generated from tissue stained for D1R (Drd1a), D2R (Drd2a), and TRAP-tdTomato.
RNAscopeAnalysis_D1_pDyn is code used to analyze csv files generated from tissue stained for D1R (Drd1a), prodynorphin (pDyn), and TRAP-tdTomato.

RNAdata is a struct that contains all the raw and processed data for D1R (Drd1a) and D2R (Drd2a) expression, as well as TRAP-tdTomato expression.
Summary is a struct that contains averaged values for all slices that are used to generate summary figures 4L-M and S5D.

pDynRNAdata is a struct that contains all the raw and processed data for D1R (Drd1a) and prodynorphin (pDyn) expression, as well as TRAP-tdTomato expression.
pDynSummary is a struct that contains averaged values for all slices that are used to generate summary figures 4N and S5E.

