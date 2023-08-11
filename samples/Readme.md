# Samples Included in the Joint-Call of Variant Discovery
- Non-Pf samples (n=4): Used for analysis outside the Posseleff project and not
listed in the meta information tables.
- MalariaGen Samples (n=7,007)
- AFRIMS_UMB Samples (n=964): This includes both previously published and new
samples.

The file `joint_call_samples_meta.tsv` provides the list of samples, along with
their meta information, included in the joint-call.

## MalariaGEN Pf6 Samples Excluded from Our Joint-Call

The Pf6 meta file lists 7,113 samples. However, 106 of these (comprising 91 OCE
and 15 WAF samples) were not incorporated into our joint-call. We faced
challenges downloading the raw reads for these samples using the ENA accession
number provided in the Pf6 metadata file, as of 2021.

The list of these 106 samples can be found in the file:
`samples_in_pf6_but_not_joint_calls.txt`.

# Analyzable Samples 
Some samples within the Joint-call set were excluded from downstream IBD-based
analyses for various reasons. These include:
- Samples exhibiting high genotype missingness.
- Polyclonal samples lacking strong evidence of a dominant clone.
- Duplicate samples that weren't previously excluded from the joint-call due to
reasons like resequencing caused by quality issues.

The file `analyzable_samples_meta.tsv` lists the samples, with their meta
information, considered in the analyzable set. Note: For the `posseleff_paper`,
only samples from ESEA (Eastern Southeast Asia) and WAF (West Africa)
populations were utilized to evaluate the impact of positive selection on
IBD-based analyses.
