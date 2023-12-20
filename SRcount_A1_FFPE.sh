# SUBMIT JOB
srun --nodes=1 --ntasks=1 --cpus-per-task=24 --mem=128G --pty /bin/bash

# PATH
fastq='/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/fastq/A1'
ref='/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/refdata/refdata-gex-GRCh38-2020-A'
cytaimage='/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/CytAssist_Image/A1.tif'
#image='/mnt/rdisk/duydao/scRNAseq/Spatial_transcriptomics/GeneSmart/datasets/Microscope_Image/A1.tiff'
manual_align_json="/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/CytAssist_Image/V42L05-352-A1.json"
probe='/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/refdata/Probe_set/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv'

# Run command FFPE
spaceranger count --id="HTI_spaceranger_count_A1" \
                 --description="Patient A1 with lung carcinoma (FFPE)" \
                 --transcriptome=$ref \
                 --probe-set=$probe \
                 --fastqs=$fastq \
                 --cytaimage=$cytaimage \
                 --loupe-alignment=$manual_align_json \
                 --slide=V42L05-352 \
                 --area=A1 \
                 --localcores=24 \
                 --localmem=128

--reorient-images=true \
## notes
## error: The argument '--reorient-images <true|false>' cannot be used with '--loupe-alignment <PATH>'
