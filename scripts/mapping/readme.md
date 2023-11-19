
# Mapping for 3DRAM-seq data 

##### Scripts are modified from TAURUS Lee et al., Nat Methods 2019 

* General steps
1. Clone scripts to desired project folder and create conda environment
2. Install Bismark, Bowtie1, PICARD, python and pysam (python module), Juicer 
3. Modify paths in map_3DRAM.py. NOTE: TAURUS_loc should specify the overarching directory that contains the ram/ and map_3DRAM.py
4. Index genome with Bismark, provide path to reference in <genome_folder> (https://felixkrueger.github.io/Bismark/bismark/genome_preparation/) 
5. Adjust parameters as required and Run 
python map_3DRAM.py -r <genome_folder> -1 <G to A converted mate> -2 <C to T converted mate> -Trim1 <Trim off first Xbp of R1> -Trim2 <Trim off first Xbp of R2> -split1 <Use first Xbp of unmapped read for 1st split reads> -split2 <Use last Xbp of unmapped read for 2nd split reads> -Threads <number of threads> -cutoff <CpG coverage cutoff to be considered for analysis> -genome <genome to use, either mm10 or hg38>
6. This generates a bash script Run_TAURUS-MH_xx.sh in the same folder that must be executed in nohup
