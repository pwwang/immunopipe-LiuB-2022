library(Seurat)
library(DropletUtils)
library(dplyr)

DATADIR = "prepared-data"
RNA_FILE = file.path(DATADIR, "GSE179994_all.Tcell.rawCounts.rds.gz")

counts <- readRDS(gzfile(RNA_FILE))
sobj <- CreateSeuratObject(counts = counts)
sobj@meta.data <- sobj@meta.data %>%
    mutate(
        Sample = sub('\\.[^.]+$', '', colnames(sobj)),
        Sample = case_when(
            Sample == 'P1.ut' ~ 'P1_pre',
            Sample == 'P2.ut' ~ 'P2_pre',
            Sample == 'P3.ut' ~ 'P3_pre',
            Sample == 'P4.ut' ~ 'P4_pre',
            Sample == 'P5.ut' ~ 'P5_pre',
            Sample == 'P6.ut' ~ 'P6_pre',
            Sample == 'P7.ut' ~ 'P7_pre',
            Sample == 'P8.ut' ~ 'P8_pre',
            Sample == 'P9.ut' ~ 'P9_pre',
            Sample == 'P10.ut' ~ 'P10_pre',
            Sample == 'P11.ut' ~ 'P11_pre',
            Sample == 'P12.ut' ~ 'P12_pre',
            Sample == 'P13.ut' ~ 'P13_pre',
            Sample == 'P14.ut' ~ 'P14_pre',
            Sample == 'P15.ut' ~ 'P15_pre',
            Sample == 'P16.ut' ~ 'P16_pre',
            Sample == 'P17.ut' ~ 'P17_pre',
            Sample == 'P18.ut' ~ 'P18_pre',
            Sample == 'P19.ut' ~ 'P19_pre',
            Sample == 'P20.ut' ~ 'P20_pre',
            Sample == 'P21.ut' ~ 'P21_pre',
            Sample == 'P22.ut' ~ 'P22_pre',
            Sample == 'P23.ut' ~ 'P23_pre',
            Sample == 'P24.ut' ~ 'P24_pre',
            Sample == 'P25.ut' ~ 'P25_pre',
            Sample == 'P26.ut' ~ 'P26_pre',
            Sample == 'P27.ut' ~ 'P27_pre',
            Sample == 'P28.ut' ~ 'P28_pre',
            Sample == 'P29.ut' ~ 'P29_pre',
            Sample == 'P30.ut' ~ 'P30_pre',
            Sample == 'P33.ut' ~ 'P33_pre',
            Sample == 'P34.ut' ~ 'P34_pre',
            Sample == 'P35.ut' ~ 'P35_pre',
            Sample == 'P1.tr.1' ~ 'P1_post_1',
            Sample == 'P1.tr.2' ~ 'P1_post_2',
            Sample == 'P1.tr.3' ~ 'P1_post_3',
            Sample == 'P10.tr.1' ~ 'P10_post_1',
            Sample == 'P13.tr.1' ~ 'P13_post_1',
            Sample == 'P13.tr.2' ~ 'P13_post_2',
            Sample == 'P19.tr.1' ~ 'P19_post_1',
            Sample == 'P29.tr.1' ~ 'P29_post_1',
            Sample == 'P30.tr.1' ~ 'P30_post_1',
            Sample == 'P33.tr.1' ~ 'P33_post_1',
            Sample == 'P35.tr.1' ~ 'P35_post_1',
            Sample == 'P36.tr.1' ~ 'P36_post_1',
            Sample == 'P37.tr.1' ~ 'P37_post_1',
            Sample == 'P38.tr.1' ~ 'P38_post_1'
        )
    )

for (sam in unique(sobj$Sample)) {
    print(sam)
    dir.create(file.path(DATADIR, sam), showWarnings = FALSE)
    obj <- subset(sobj, subset = Sample == sam)
    cts <- GetAssayData(obj, layer = 'counts')
    barcodes = sub('^.+\\.', '', colnames(cts))
    write10xCounts(file.path(DATADIR, sam, 'rna'), cts, barcodes = barcodes, overwrite = TRUE)
}