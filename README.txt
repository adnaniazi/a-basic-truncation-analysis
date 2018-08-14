# A basic sites truncation analysis
Date: 13 Aug 2018
Dataset: DMS-lig (Barcode 10 primarily)
Purpose: To asses if the a-basic sites causes artificial truncation of reads

## Thoughts
1. Get the unclassified data
2. align it to the backbone
3. fish out the primary aligned reads
4. get the channel wise start and ends of the reads for both barcode 10 and unclassified data
5. filter the start ends so that only those that are 5 seconds or less apart remain
6. search the read ids of the filtered reads in the SAM file and get the start POS argument
7. Annotate reads with the alignment start position.
8. Plot this information in bokeh as time segments with annotation
