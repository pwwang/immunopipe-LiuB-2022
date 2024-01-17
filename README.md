# immunopipe-LiuB-2022

Reanalysis of the scRNA-seq and scTCR-seq data from [Liu, B., et al. 2022](https://www.nature.com/articles/s43018-021-00292-8) using [immunopipe](https://github.com/pwwang/immunopipe).

> [Liu, B., Hu, X., Feng, K., Gao, R., Xue, Z., Zhang, S., Zhang, Y., Corse, E., Hu, Y., Han, W., & Zhang, Z. (2022). Temporal single-cell tracing reveals clonal revival and expansion of precursor exhausted T cells during anti-PD-1 therapy in lung cancer. Nature Cancer, 3(1), 108-121.](https://www.nature.com/articles/s43018-021-00292-8)

## Data preparation

The data was downloaded from [GSE179994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179994), where 47 samples from 38 patients were sequenced with scRNA-seq and scTCR-seq.

See `prepare-data.sh` for details.

## Configuration

> [!NOTE]
> This is not a replication of the original paper, primarily due to the irreproducibility of the clustering results. This is a reanalysis of the data using [`immunopipe`](https://github.com/pwwang/immunopipe), showing the potential of the pipeline similar analyses listed in the paper.
>

The scRNA-seq data was prepared by following the instruction of the original paper. While selecting T cells, the TCR information was ignored (not mentioned in the paper). Instead, `CD3E`, `CD3D`, `CD3G` and `NKG7` were used to select T cells (and NK cells). The clusters were not annotated, so the clusters were named as `C1`, `C2`, etc. For further analysis, the clusters were selected by the expression of marker genes.

See details at `Immunopipe.config.toml`.

## Results/Reports

You can find the results in the `Immunopipe-output` directory.

The report can be found at [https://imp-liub-2022.pwwang.com/REPORTS](https://imp-liub-2022.pwwang.com/REPORTS).
