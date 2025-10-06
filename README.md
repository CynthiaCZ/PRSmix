# PRSmix
 Multi-Trait Polygenic Risk Score and Proteomic Analysis of COPD Exacerbations

## Study Schematic 
![Workflow][workflow_chart]


### PRS Analysis (PRS_analysis/)
<ul>
  <li> Implements construction and testing of multi-trait PRS (PRSmulti) across multiple cohorts, including COPDGene, ECLIPSE, MGBB, UK Biobank, and All of Us. </li>
  <ul>
    <li> PRSmix_COPDGene_NHW_train.Rmd: training PRSmulti with PRSmix+ in COPDGene non-Hispanic White participants </li> 
    <li> *_test.Rmd: testing association of PRSmulti with logistic regression for COPD status and negative binomial for exacerbation counts in testing cohorts </li> 
    <li> PRSmix_meta_analysis.Rmd: meta-analysis of PRSmulti results across cohorts and generation of tables and figures </li> 
  </ul>
</ul>

### Proteomic Analysis (proteomic_analysis/)
<ul>
  <li> Validate genetically predicted protein levels (via S-PrediXcan) with measured plasma protein levels.
    <ul>
    <li> proteomics_COPDGene.Rmd: analysis in COPDGene </li> 
    <li> proteomics_UKB.Rmd: analysis in UK Biobank </li> 
  </ul>
  <li> PRSmix_meta_analysis.Rmd: concordance check between the two cohorts and table generation </li> 
</ul>

### Clinical Outcome Definition (define_COPD_exac_ICD.py)
<ul>
  <li> Defines COPD cases and exacerbation events with ICD codes.
</ul>

<!-- MARKDOWN LINKS & IMAGES -->
[workflow_chart]: ./workflow_chart.png

