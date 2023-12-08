# PATH
fastq='/mnt/d4t/SPATIAL/GeneSmart/datasets/fastq/10X'
ref='/mnt/d4t/SPATIAL/GeneSmart/datasets/refdata/refdata-gex-GRCh38-2020-A'
cytaimage='/mnt/d4t/SPATIAL/GeneSmart/datasets/CytAssist_Image/10X.tif'
#image='/mnt/rdisk/duydao/scRNAseq/Spatial_transcriptomics/GeneSmart/datasets/Microscope_Image/A1.tiff'
manual_align_json="/mnt/d4t/SPATIAL/GeneSmart/datasets/CytAssist_Image/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_alignment_file.json"
probe='/mnt/d4t/SPATIAL/GeneSmart/datasets/refdata/Probe_set/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_probe_set.csv'

# Run command FFPE
spaceranger count --id="Sample_spaceranger_count_10X" \
                 --description="Patient with lung squamous cell (FFPE)" \
                 --transcriptome=$ref \
                 --probe-set=$probe \
                 --fastqs=$fastq \
                 --cytaimage=$cytaimage \
                 --loupe-alignment=$manual_align_json \
                 --slide=V42A20-354 \
                 --area=D1 \
                 --localcores=16 \
                 --localmem=50