# Command 
```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_empirical/05_ibd_ne_ifm.nf \
        -profile sge \
        --vcf '/local/chib/toconnor_grp/bing/posseleff_empirical/runs/r231118/input/*.vcf.gz' \
        --peak_validate_meth ihs \
        --ibdne_no_diploid_convertion true \
        -resume
```
# input
see details in `input/input_prep_notes.txt`

# version
This run was done after the follow update of ibdutils 

commit with message "allow calc peak impact and provide arg to only remove high impact peaks"



# plots

see `posseff` repo (internal use only)
