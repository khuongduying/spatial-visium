# PATH
fastq='/mnt/d4t/SPATIAL/GeneSmart/datasets/fastq/D1'
ref='/mnt/d4t/SPATIAL/GeneSmart/datasets/refdata/refdata-gex-GRCh38-2020-A'
cytaimage='/mnt/d4t/SPATIAL/GeneSmart/datasets/CytAssist_Image/D1.tif'
#image='/mnt/rdisk/duydao/scRNAseq/Spatial_transcriptomics/GeneSmart/datasets/Microscope_Image/D1.tiff'
#manual_align_json="/mnt/d4t/SPATIAL/GeneSmart/datasets/CytAssist_Image/V42L05-352-A1.json"
probe='/mnt/d4t/SPATIAL/GeneSmart/datasets/refdata/Probe_set/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv'

# Run command FFPE (CytAssist)
spaceranger count --id="HTI_spaceranger_count_D1" \
                 --description="Patient D1 with lung carcinoma (FFPE)" \
                 --transcriptome=$ref \
                 --probe-set=$probe \
                 --fastqs=$fastq \
                 --cytaimage=$cytaimage \
                 --slide=V42L05-352 \
                 --area=D1 \
                 --reorient-images=true \
                 --localcores=16 \
                 --localmem=50