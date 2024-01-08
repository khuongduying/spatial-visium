# SUBMIT JOB
srun --nodes=1 --ntasks=1 --cpus-per-task=24 --mem=128G --pty /bin/bash

# PATH
fastq='/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/fastq/10X'
ref='/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/refdata/refdata-gex-GRCh38-2020-A'
cytaimage='/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/CytAssist_Image/10X_crop2.tif'
#image='/mnt/rdisk/duydao/scRNAseq/Spatial_transcriptomics/GeneSmart/datasets/Microscope_Image/A1.tiff'
manual_align_json="/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/CytAssist_Image/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_alignment_file.json"
probe='/mnt/rdisk/duydao/SPATIAL/GeneSmart/datasets/refdata/Probe_set/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_probe_set.csv'

# Run command FFPE
spaceranger count --id="Sample_spaceranger_count_10X" \
                 --description="Patient with lung squamous cell (FFPE)" \
                 --transcriptome=$ref \
                 --probe-set=$probe \
                 --fastqs=$fastq \
                 --cytaimage=$cytaimage \
                 --slide=V42A20-354 \
                 --reorient-images=true \
                 --area=D1 \
                 --localcores=24 \
                 --localmem=128

--loupe-alignment=$manual_align_json \
                 --reorient-images=true \