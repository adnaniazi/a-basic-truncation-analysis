Extract for each read the mean and median of the signal level (or any other factor used by Tombo/Albacore to normalise raw data) and plot the distribution of median per read as a box plot/violin plot comparing the three treatments (untreated, DMS, abasic)

Date 17 Aug 2016
Dataset: DMS lig


Recipe:
1. Will do it in R
2. Collect raw data from all reads one by ones 
3. Only use data which is actually basecalled from the base calling start position
4. pa-normalize the data, and median normalize the data
5. calculate mean and median for pa-normalized and median normalized data
6. collect the mean and median for both normalization

pa_norm_mean_bc08
pa_norm_median_bc08

pa_norm_mean_bc09
pa_norm_median_bc09

pa_norm_mean_bc10
pa_norm_median_bc10

md_norm_mean_bc08
md_norm_median_bc08

md_norm_mean_bc09
md_norm_median_bc09

md_norm_mean_bc10
md_norm_median_bc10