library(dplyr)
library(tidyr)

DATADIR = "prepared-data"
TCR_FILE = file.path(DATADIR, "GSE179994_all.scTCR.tsv.gz")

tcr <- read.table(
    gzfile(TCR_FILE),
    sep = '\t',
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE)

tcr <- tcr %>%
    mutate(
        barcode = sub('^.+\\.', '', CellName),
        is_cell = 'true',
        contig_id = paste(`Identifier(Alpha1)`, `Identifier(Beta1)`, sep = '!!!'),
        high_confidence = 'true',
        length = paste(`Length(Alpha1)`, `Length(Beta1)`, sep = '!!!'),
        chain = paste(substring(`V_gene(Alpha1)`, 1, 3), substring(`V_gene(Beta1)`, 1, 3), sep = '!!!'),
        v_gene = paste(`V_gene(Alpha1)`, `V_gene(Beta1)`, sep = '!!!'),
        d_gene = '',
        j_gene = paste(`J_gene(Alpha1)`, `J_gene(Beta1)`, sep = '!!!'),
        c_gene = '',
        full_length = 'true',
        productive = 'true',
        fwr1 = '',
        fwr1_nt = '',
        cdr1 = '',
        cdr1_nt = '',
        fwr2 = '',
        fwr2_nt = '',
        cdr2 = '',
        cdr2_nt = '',
        fwr3 = '',
        fwr3_nt = '',
        cdr3 = paste(`CDR3(Alpha1)`, `CDR3(Beta1)`, sep = '!!!'),
        cdr3_nt = paste(`CDR3_nt(Alpha1)`, `CDR3_nt(Beta1)`, sep = '!!!'),
        fwr4 = '',
        fwr4_nt = '',
        reads = paste(`nRead(Alpha1)`, `nRead(Beta1)`, sep = '!!!'),
        umis = paste(`nUMI(Alpha1)`, `nUMI(Beta1)`, sep = '!!!'),
        raw_clonotype_id = '',
        raw_consensus_id = '',
        exact_subclonotype_id = 1,
        Sample = sub('\\.[^.]+$', '', CellName),
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
    ) %>%
    separate_rows(
        contig_id,
        length,
        chain,
        v_gene,
        j_gene,
        cdr3,
        cdr3_nt,
        reads,
        umis,
        sep = '!!!'
    ) %>%
    select(
        barcode,
        is_cell,
        contig_id,
        high_confidence,
        length,
        chain,
        v_gene,
        d_gene,
        j_gene,
        c_gene,
        full_length,
        productive,
        fwr1,
        fwr1_nt,
        cdr1,
        cdr1_nt,
        fwr2,
        fwr2_nt,
        cdr2,
        cdr2_nt,
        fwr3,
        fwr3_nt,
        cdr3,
        cdr3_nt,
        fwr4,
        fwr4_nt,
        reads,
        umis,
        raw_clonotype_id,
        raw_consensus_id,
        exact_subclonotype_id,
        Sample
    ) %>%
    group_by(Sample) %>%
    group_walk(
        function(.x, .y) {
            sam = .y$Sample[1]
            print(sam)
            if (is.na(sam)) { return() }
            odir = file.path(DATADIR, sam, 'tcr')
            dir.create(odir, recursive = TRUE, showWarnings = FALSE)
            write.table(
                .x,
                file = file.path(odir, 'filtered_contig_annotations.csv'),
                sep = ',',
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE
            )
        }
    )