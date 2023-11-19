# Scripts for paired single molecule accessibiity analysis 

 * generate_methCalls.R - Summarize methylation calls per read, giving each row the read name, the chromosomal and base pair position of the CpG or GpC call

 * filter_methCalls.R - Takes in per chromosome fst file, in a parallel per chromosome fashion: define filter intervals (eg ctcf peak motifs) and associate read1 and read2 to respective to these intervals through binary search. Reads to intervals exceeding defined window size are removed. Final reads bearing the same read name are merged. 

 * filter_methCalls2_CpGGpC.R - Similar to filter_methCalls.R but use this script to match and filter CpG methy marks as one anchor and GpC meth marks on the next anchor 

 * filter_methCalls2_control.R - Similar to filter_methCalls.R but use this script to match and filter using control regions as anchors 
