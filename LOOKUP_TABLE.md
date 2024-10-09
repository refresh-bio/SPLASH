# SPLASH: lookup table

### General idea

This software is designed to query k-mers from extendors/compactors/fasta files against reference sequences.

The first step is to build an index for a given set of reference sequences.
This index is further used to perform queries.

### Building the index
To build an index one should use `build_lookup_table.py` script.
Let's use three simple example files:

`1.fa`:
```
>h1-cat_1
ACGACGACGACG
>h2-cat_2
ACGACGACGACT
>h3-cat_2
ACGACGACGACT
>h4-cat_3
ACGACTGCGACG
>h5-cat_4
TGATGATGATGA
>h6-cat_4
TGATGATGATGA
>h7-cat_4
TGATGATGATGA
>h8-poly
AAAAAAAACCCC
```
,
`2.fa`:
```
>h1-cat_1
ACGCCGACGACG
>h2-cat_2
TCGCCGACGACG
>h3-cat_2
TCGCCGACGACG
>h4-cat_2
TCGCCGACGACG
>h5-cat_3
ACGACTGCGACG
```
, and `3.fa`:
```
>h1-cat_3
ACGACTGCGACG
>h2-cat_3
ACGACTGCGACG
>h3-cat_4
TGATGATGATGA
>h4-cat_4
TGATGATGATGA
>h5-cat_4
TGATGATGATGA
```
Let's put paths to these files in `input.txt`:
```
1.fa
2.fa
3.fa
```

To build an index of 12-mers without any poly sequences of length at least 6 one may use:
```
./build_lookup_table.py --poly_ACGT_len 6 --kmer_len 12 input.txt
```
The resulting index is a single file: `lookup.slt` (slt extension stands for splash lookup table).

The default output file name may be overwritten with `--outname`.


The full usage of `build_lookup_table.py` is:
```
  -h, --help            show this help message and exit
  --transcriptomes TRANSCRIPTOMES
                        path to additional (optional) file where transcriptome input FASTA files are defined, the format is: per each line path to fasta file (default: )
  --outname OUTNAME     prefix of output file names (default: lookup.slt)
  --kmer_len KMER_LEN   k-mer length (default: 15)
  --bin_path BIN_PATH   path to a directory where kmc, kmc_tools (default: .)
  --n_threads N_THREADS
                        number of threads (default: 8)
  --poly_ACGT_len POLY_ACGT_LEN
                        all k-mers containing polyACGT of this length will be filtered out (0 means no filtering) (default: 0)
  --category_3_threshold CATEGORY_3_THRESHOLD
                        accept k-mer in category 3 if its present in a given file <=category_3_threshold times (default: 1)
  --precomputed_sbwt PRECOMPUTED_SBWT
                        path to precomputed sbwt index (if set lookup_table will use it instead of building own sbwt - must be build for exactly the same set of input files!) (default: )
  --dont_clean_up       if set then intermediate files will not be removed (default: False)
```

### Extending existing index
This is currently not supported and the whole index should be rebuilt.

### Querying the index
To query the index one should use `lookup_table` program in `query` mode.
It expects three parameters:

Input:
 - `input` - input file with sequences to be queried
 - `lookup` - path to a index being a result of running `build_lookup_table.py`
 
Output:
 - `output` - path to the output text file

It also provides additional configuration:
 - `--input_fmt <format>`  - input format, one of: `fasta`, `extendors`, `compactors` (default: fasta)
 - `--report_fmt <string>` - format of the detailed report, one of: plain (base format), ids (just ids), empty (do not report, useful when stats enabled with stats_fmt) (default: plain)
 - `--stats_fmt <string>`  - for each query report stats (#k-mers per category) one of: empty (don't print stats) or with_stats (print stats) (default: empty)
 - `--output_fmt <string>` - output format, one of: txt (one line for each query), extendors (only for extendors input, adds additional columns to the extendors file), compactors (only for compactors input, adds additional columns to the compactors file) (default: txt)
 - `--kmer_skip <int>`     - for each query the next queried k-mer start on position <kmer_skip> + 1
 - `--truncate_paths`      - if set the files paths will be truncated to file names only

Let's say we have a file `4.fa`:
```
>
ANACGACGACGACGTACGCCGACGACG
>
ACGACGACGACT
>
AAAAAACGACTGCGACG
>
TGATGATGATGA
```
that we want to use to perform queries.
To do this one may use the following command:
```
./lookup_table query 4.fa lookup.slt o.txt
```
The resulting file `out.txt`:
```
N, N, 1: ("1.fa","h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, 1: ("2.fa","h1-cat_1")
2: ("1.fa",2)
P, U, U, U, U, 3: ["1.fa": ["h4-cat_3"], "2.fa": ["h5-cat_3"]]
4
```
Each line corresponds to a single input sequence.

Categories of k-mers are defined as as below:
1. k-mer is unique to a single fasta and also appears once in it (report the type and the fasta where it is unique)
2. k-mer is unique to a single fasta appears multiple times in it (report the type and the fasta where it appears and report # times it appears)
3. k-mer is unique in that it appears once in one fasta but also appears in other fastas (report the type and the header where it appears)
4. Appears >1 in >1 fasta
 
For each k-mer we report its state using the following convention:
 - `1: ("<file_path>","<header>")` - k-mer is in category 1 and its file and header is reported
 - `2: ("<file_path>",<cnt>)` - k-mer is in category 2 and its file and count is reported
 - `3: ["<file_1>": ["<header_in_file1>"], "<file_2>": ["header_in_file2"], ...]` - k-mer is in category 3 and for each file it occurs only once report the file and header
 - `4` - k-mer is in category 4
 - `U` - unknown, means that a given k-mers was not observed in any of the references used to build an index
 - `P` - means k-mer contains polyACGT sequence of length at least equal to the one used to build an index
 - `N` - means this k-mer contains at least one non-ACGT symbol, note that a single `N` in a query sequence may produce up to `k` invalid k-mers
 
The states are comma separated. The number of states per each sequence is len(sequence) - k + 1.

 ### Extension

It is possible to increase the acceptance into the category three threshold using `--category_3_threshold` parameter in `build_lookup_table.py` (default 1). In such a case a list of headers for a given k-mer and given file in category three may be longer. In other words, some of k-mers instead of being classified into category 4 will be classified into category 3.

### Using transcriptomes files
To use transcriptomes fasta file one should use `--transcriptomes` switch. Each transcriptome file is handled a little differently than other files.

**Warning** The transcriptome files should be probably only contained in the file pointed by `--transcriptomes` switch and not in the regular input file.
If it is in both it will be read twice and treated once as transcriptome and once as non-transcriptome.
If one needs to only use transcriptome files the regular input file may be empty, but still must be provided to the build script.

**The rules to assign k-mers to categories are as follows:**

1. If k-mer exists in only single variant of single gene (so basicaly in single header) and

    a. it exists only in a single file it is stored in category 1 with gene variant name, or

    b. it exists in at least one other file, it is stored in category 3 with gene variant name.
2. If k-mer exists only in a single gene, but in multiple variants and

    a. it exists only in a single file it is stored in category 1 with gene name, or

    b. it exists in at least one other file, it is stored in category 3 with gene name.
3. If k-mer exists in multiple genes and it exists only in a single file it is stored in category 2 with transcriptome file name and number of its occurences.

#### Extracting gene name from fasta headers
There are currently two ways of extracting gene name.
If there is `(` character in the header The parentheses method is used in the opposite case the dot/dash method is used.
#### The parentheses method (like in GRCh38_latest_rna.fna)
 1. If there is only one `(something)` in the header this is a gene name
 2. In the opposite case
  2a. If there is a `transcript variant` in the header I remove the line starting with this text (`transcript variant`) to the end <- because for some headers there were `()` after transcript variant which I believe are not gene names
  2b. I read the content of the last `(something)` in the (optionally trimmed in 2a) line and this is a gene name
#### The dot/dash method (like in CAT_liftoff_genes.fa)
 1. Split the header by whitespace and take the first **part** for the next steps.
 2. If there is `.` in the **part** the gene name is everything from the beginning to the `.`
 3. In the opposite case if there is `-` in the **part** the gene name is everything from the beginning to the `-`
 4. In the opposite case the gene name is the **part**
 

I'm not sure if this parsing is fine and if this is consistent for other transcriptome files.
Please feel free to post issues in that matter.

**Transcriptome handling is currently experimental and may change, also if there are any performance/memory usage-related issues please post an issue with details**

#### Example for transcriptomes
Let's say we want to extend the previous example of a transcriptome file.
Let's assume we have a file `genes.fa` with content
```
>gene (ABC), mRNA
ACGTTTACGACGT
>gene (ABC), variant 1, mRNA
ACGTTTACGACGC
>gene(ABD), variant 2, mRNA
ACGACTGCGACG
>gene(ABD), variant 3, mRNA
ACACACACACAC
>gene(AAAC), variant 1, mRNA
ACACACACACAC
>gene(AAAC), variant 2, mRNA
ACACACACACAC
```
And we have a file `transcriptomes.txt` with the following content:
```
genes.fa
```
We use the following command to build a lookup table:
```
./build_lookup_table.py --poly_ACGT_len 6 --kmer_len 12 --transcriptomes transcriptomes.txt input.txt
```
And then query the transcriptome file itself:
```
./lookup_table query genes.fa lookup.slt o.txt
```
The content of `o.txt` will be:
```
3: ["genes.fa": ["ABC"]], 1: ("genes.fa", "gene (ABC), mRNA")
3: ["genes.fa": ["ABC"]], 1: ("genes.fa", "gene (ABC), variant 1, mRNA")
3: ["1.fa": ["h4-cat_3"], "2.fa": ["h5-cat_3"], "genes.fa": ["gene(ABD), variant 2, mRNA"]]
3: ["genes.fa": ["AAAC", "gene(ABD), variant 3, mRNA"]]
3: ["genes.fa": ["AAAC", "gene(ABD), variant 3, mRNA"]]
3: ["genes.fa": ["AAAC", "gene(ABD), variant 3, mRNA"]]

```
### Supported input formats for a query
In the above example, the fasta file was used for the query, but one may also use compactors (`--input_fmt compactors`) file or main SPLASH output file `result.after_correction.scores.tsv` (`--input_fmt extendors`).
In the case of compactors, each compactor is treated as a separate sequence and results in one line in the output.
In the case of extendors, for each line in the input file (`result.after_correction.scores.tsv`) the anchor is taken and there are as many sequences built as there are targets. 
For example, if there are three targets stored for each anchor there will be up to three lines in the output file for these targets (maybe less if a given anchor has fewer targets).

### Formatting the output
The result of a query may be extended using `--stats_fmt with_stats`.
In this case for each query after the standard result there will also be a summary of k-mers to category assignments.
The format is: `<categor>: <number of k-mers in query in this category>`
For example when the query is:
`../bin/lookup_table query --stats_fmt with_stats 4.fa lookup.slt o.txt`
The output is:
```
N, N, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, 1: ("2.fa", "h1-cat_1")	1: 2, 2: 0, 3: 0, 4: 0, U: 12, P: 0, N: 2
2: ("1.fa", 2)	1: 0, 2: 1, 3: 0, 4: 0, U: 0, P: 0, N: 0
P, U, U, U, U, 3: ["1.fa": ["h4-cat_3"], "2.fa": ["h5-cat_3"]]	1: 0, 2: 0, 3: 1, 4: 0, U: 4, P: 1, N: 0
4	1: 0, 2: 0, 3: 0, 4: 1, U: 0, P: 0, N: 0
```
If only stats are of interest one may suppress the standard list with `--report_fmt empty`.

**Warning**: using `--report_fmt empty` with `--stats_fmt empty` (or not setting `--stats_fmt` at all) will result in output containing #queries empty lines.

#### Output in SPLASH TSV format
When the input format is `extendors` (`result.after_correction.scores.tsv`) it is possible to set `--output_fmt extendors`. In this case additional columns are added. 
For each `most_freq_target_<i>` column there will be a new column `most_freq_extendor_<i>_query_res` if `report_fmt` is `plain` or `ids` and `most_freq_extendor_<i>_query_stats` if `stats_fmt` is `with_stats`.

**Warning**: If both (`report_fmt` and `stats_fmt`) are `empty` the output will be exactly the same as input file.

#### Output in COMPACTORS TSV format
When the input format is `compactors` it is possible to set `--output_fmt compactors`. In this case additional columns are added. 
 - `query_res` if `report_fmt` is `plain` or `ids` 
 - `query_stats` if `stats_fmt` is `with_stats`

**Warning**: If both (`report_fmt` and `stats_fmt`) are `empty` the output will be exactly the same as input file.

#### Example for `--input_fmt extendors`
Lets assume we have following `result.after_correction.scores.tsv` SPLASH output:
```
anchor	pval_opt	effect_size_bin	pval_base	M	anch_uniqTargs	number_nonzero_samples	target_entropy	avg_no_homopolymer_targets	avg_hamming_distance_max_target	avg_hamming_distance_all_pairs	avg_edit_distance_max_target	avg_edit_distance_all_pairs	most_freq_target_1	cnt_most_freq_target_1	most_freq_target_2	cnt_most_freq_target_2	pval_opt_corrected
CAACGACGACGACGTACACGACGACGA	1.26129e-21	0.213304	0.0218488	2639	54	342	1.27155	0	6.25161	7.40047	5.41455	6.40716	CGTACGCCGACGACGTGCTGAATGGCT	1520	CGTACGCCGACGTCGTCCTGTCGGGGT	1049	6.10416177521e-19
AGAGATGCTGGTGGCAGGGCAACGACG	9.0965e-07	0.14895	1	1333	36	177	1.27046	0	0.414854	0.535504	0.414854	0.535504	ACGACTTGGTTGTATCAACTATGAAGA	804	ACGACCTGGTTGTATCAACTATGAAGC	489	6.77285923852e-05
AGATGCTGGTGGCAGGGCATGAGAAAA	5.99838e-09	0.205298	1	1397	31	184	1.22632	0	0.778812	0.997034	0.778812	0.994373	AACGACTGCGACGCAACTATGAAGAGC	856	GCAATGGTTGTATCAACTATGAAGCGT	504	5.04867265421e-07
AGGACAGCAATGGTTGTGATGATGATG	1.6011e-22	0.234194	1	2116	48	275	1.31055	0	6.14272	6.99537	5.26701	5.99881	AAGAGCTCGTCCGCATGGTGCTGAATG	1171	AAGCGTTTGTGAGGCATATCCTGTCGG	881	1.33768977104e-19
```
We can run the following command:
```
./lookup_table query --input_fmt extendors --output_fmt extendors --stats_fmt with_stats result.after_correction.scores.tsv lookup.slt result.after_correction.scores.queried.tsv 
```
 - `--input_fmt extendors` says that the input is in SPLASH output format
 - `--output_fmt extendors` says that we want this same file with additional columns as the output
 - `--stats_fmt with_stats` says we want to have also summary for each query

In such a case the resulting file `result.after_correction.scores.queried.tsv` is:
```
anchor	pval_opt	effect_size_bin	pval_base	M	anch_uniqTargs	number_nonzero_samples	target_entropy	avg_no_homopolymer_targets	avg_hamming_distance_max_target	avg_hamming_distance_all_pairs	avg_edit_distance_max_target	avg_edit_distance_all_pairs	most_freq_target_1	cnt_most_freq_target_1	most_freq_target_2	cnt_most_freq_target_2	pval_opt_corrected	most_freq_extendor_1_query_res	most_freq_extendor_1_query_stats	most_freq_extendor_2_query_res	most_freq_extendor_2_query_stats
CAACGACGACGACGTACACGACGACGA	1.26129e-21	0.213304	0.0218488	2639	54	342	1.27155	0	6.25161	7.40047	5.41455	6.40716	CGTACGCCGACGACGTGCTGAATGGCT	1520	CGTACGCCGACGTCGTCCTGTCGGGGT	1049	6.10416177521e-19	U, U, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, U, U, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, 1: ("2.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U	1: 3, 2: 0, 3: 0, 4: 0, U: 40, P: 0, N: 0	U, U, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, U, U, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 2, 2: 0, 3: 0, 4: 0, U: 41, P: 0, N: 0
AGAGATGCTGGTGGCAGGGCAACGACG	9.0965e-07	0.14895	1	1333	36	177	1.27046	0	0.414854	0.535504	0.414854	0.535504	ACGACTTGGTTGTATCAACTATGAAGA	804	ACGACCTGGTTGTATCAACTATGAAGC	489	6.77285923852e-05	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 2: ("1.fa", 2), U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 1, 3: 0, 4: 0, U: 42, P: 0, N: 0	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 43, P: 0, N: 0
AGATGCTGGTGGCAGGGCATGAGAAAA	5.99838e-09	0.205298	1	1397	31	184	1.22632	0	0.778812	0.997034	0.778812	0.994373	AACGACTGCGACGCAACTATGAAGAGC	856	GCAATGGTTGTATCAACTATGAAGCGT	504	5.04867265421e-07	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, P, P, P, P, P, P, P, U, U, U, U, 3: ["1.fa": ["h4-cat_3"], "2.fa": ["h5-cat_3"]], U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 1, 4: 0, U: 35, P: 7, N: 0	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 43, P: 0, N: 0
AGGACAGCAATGGTTGTGATGATGATG	1.6011e-22	0.234194	1	2116	48	275	1.31055	0	6.14272	6.99537	5.26701	5.99881	AAGAGCTCGTCCGCATGGTGCTGAATG	1171	AAGCGTTTGTGAGGCATATCCTGTCGG	881	1.33768977104e-19	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 4, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 1, U: 42, P: 0, N: 0	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 4, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 1, U: 42, P: 0, N: 0
```
If we don't define `--output_fmt extendors` i.e. we use:
```
./lookup_table query --input_fmt extendors --stats_fmt with_stats result.after_correction.scores.tsv lookup.slt result.after_correction.scores.queried.txt 
```
we will have each result for extendor (anchor+target) as a separate line:
```
U, U, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, U, U, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, 1: ("2.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U	1: 3, 2: 0, 3: 0, 4: 0, U: 40, P: 0, N: 0
U, U, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, U, U, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 2, 2: 0, 3: 0, 4: 0, U: 41, P: 0, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 2: ("1.fa", 2), U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 1, 3: 0, 4: 0, U: 42, P: 0, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 43, P: 0, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, P, P, P, P, P, P, P, U, U, U, U, 3: ["1.fa": ["h4-cat_3"], "2.fa": ["h5-cat_3"]], U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 1, 4: 0, U: 35, P: 7, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 43, P: 0, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 4, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 1, U: 42, P: 0, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 4, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 1, U: 42, P: 0, N: 0
```

The software will look for columns with names `anchor` and `most_freq_target_<number>` and build extendor for each pair of anchor + target.

#### Example for `--input_fmt compactors`
Lets say we have a compactors output file `compactors.queried.txt`:
```
anchor  compactor       support exact_support   extender_specificity    num_extended
AAAAACAAGACTGGGGCTGCTCCCATC     ACGACGACGACTANACGACGACGACGTACGCCGACGACGAAAAAACGACTGCGACGTGATGATGATGAAAAAACAAGACTG       570381  358915  0.944578        0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAAACAAGACTGGGGCTGCTCCCATCATTGATGTGGTGCGATCGGGCTACTACAAAGTTCTGGGAAGGGAAAGCTCCCAA       916     193     0.870968        0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAAACAAGACTGGGGCTGCTCCCATCATTGATGTGGTGCGATCGGGCTACTACAAGTTCTGGGAAAGGGAAAGCTCCCAA       446     298     0.001254        0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAAACAAGACTGGGGCTGCTCCCATCATACGACGTACGCCGACGAGCTACTACAAAGTTCTGGGAAAAGGGAAAGCTCCC       300     117     1.000000        0
AAAAACAAGACTGGGGCTGCTCCCATC     AAACGACGACGACTGGCTGCTCCCATCATTGATGTGGTGCGATCGGGCTACTACAAAAGTTCTGGGAAAGGGAAAGCTCCC       246     148     0.307692        0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAAACAAGACTGGGGCTGCTCCCATCATTGATGTGGTGCGATCGGCTACTACAAAGTTCTGGGAAAGGGAAAGCTCCCAA       165     80      0.001254        0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAATGATGATGATGACTGCTCCCATCATTGATGATGATGATGAGGCTACTACAACGACGACGACTAGGGAAAGCTCCCAA       105     49      0.001254        0
```
We can run the following command:
```
./lookup_table query --input_fmt compactors --output_fmt compactors --stats_fmt with_stats compactors.out.tsv lookup.slt compactors.queried.tsv
```
 - `--input_fmt compactors` says that the input is in COMPACTORS output format
 - `--output_fmt compactors` says that we want this same file with additional columns as the output
 - `--stats_fmt with_stats` says we want to have also summary for each query

In such a case the resulting file `compactors.queried.tsv` is:
```
anchor  compactor       support exact_support   extender_specificity    num_extended	query_res	query_stats
AAAAACAAGACTGGGGCTGCTCCCATC     ACGACGACGACTANACGACGACGACGTACGCCGACGACGAAAAAACGACTGCGACGTGATGATGATGAAAAAACAAGACTG       570381  358915  0.944578        0	2: ("1.fa", 2), U, N, N, N, N, N, N, N, N, N, N, N, N, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, 1: ("2.fa", "h1-cat_1"), U, U, U, U, U, P, P, P, P, P, P, P, U, U, U, U, 3: ["1.fa": ["h4-cat_3"], "2.fa": ["h5-cat_3"]], U, U, U, U, U, U, U, U, U, U, U, 4, U, U, U, U, P, P, P, P, P, P, P, U, U	1: 2, 2: 1, 3: 1, 4: 1, U: 39, P: 14, N: 12
AAAAACAAGACTGGGGCTGCTCCCATC     AAAAACAAGACTGGGGCTGCTCCCATCATTGATGTGGTGCGATCGGGCTACTACAAAGTTCTGGGAAGGGAAAGCTCCCAA       916     193     0.870968        0	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 70, P: 0, N: 0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAAACAAGACTGGGGCTGCTCCCATCATTGATGTGGTGCGATCGGGCTACTACAAGTTCTGGGAAAGGGAAAGCTCCCAA       446     298     0.001254        0	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 70, P: 0, N: 0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAAACAAGACTGGGGCTGCTCCCATCATACGACGTACGCCGACGAGCTACTACAAAGTTCTGGGAAAAGGGAAAGCTCCC       300     117     1.000000        0	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 70, P: 0, N: 0
AAAAACAAGACTGGGGCTGCTCCCATC     AAACGACGACGACTGGCTGCTCCCATCATTGATGTGGTGCGATCGGGCTACTACAAAAGTTCTGGGAAAGGGAAAGCTCCC       246     148     0.307692        0	U, U, 2: ("1.fa", 2), U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 1, 3: 0, 4: 0, U: 69, P: 0, N: 0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAAACAAGACTGGGGCTGCTCCCATCATTGATGTGGTGCGATCGGCTACTACAAAGTTCTGGGAAAGGGAAAGCTCCCAA       165     80      0.001254        0	U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 70, P: 0, N: 0
AAAAACAAGACTGGGGCTGCTCCCATC     AAAATGATGATGATGACTGCTCCCATCATTGATGATGATGATGAGGCTACTACAACGACGACGACTAGGGAAAGCTCCCAA       105     49      0.001254        0	U, U, U, U, 4, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 4, U, U, 4, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 2: ("1.fa", 2), U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 1, 3: 0, 4: 3, U: 66, P: 0, N: 0
```
If we don't define `--output_fmt compactors` i.e. we use:
```
./lookup_table query --input_fmt compactors --stats_fmt with_stats compactors.out.tsv lookup.slt compactors.queried.txt
```
we will have each result for compactor as a separate line:
```
2: ("1.fa", 2), U, N, N, N, N, N, N, N, N, N, N, N, N, 1: ("1.fa", "h1-cat_1"), U, U, U, U, U, U, U, U, U, U, U, U, 1: ("2.fa", "h1-cat_1"), U, U, U, U, U, P, P, P, P, P, P, P, U, U, U, U, 3: ["1.fa": ["h4-cat_3"], "2.fa": ["h5-cat_3"]], U, U, U, U, U, U, U, U, U, U, U, 4, U, U, U, U, P, P, P, P, P, P, P, U, U	1: 2, 2: 1, 3: 1, 4: 1, U: 39, P: 14, N: 12
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 70, P: 0, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 70, P: 0, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 70, P: 0, N: 0
U, U, 2: ("1.fa", 2), U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 1, 3: 0, 4: 0, U: 69, P: 0, N: 0
U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 0, 3: 0, 4: 0, U: 70, P: 0, N: 0
U, U, U, U, 4, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 4, U, U, 4, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, 2: ("1.fa", 2), U, U, U, U, U, U, U, U, U, U, U, U, U, U, U	1: 0, 2: 1, 3: 0, 4: 3, U: 66, P: 0, N: 0
```
### Canonical k-mers
Lexicographically smaller of the k-mer and its reverse complement is the so-called canonical k-mer.
In the index, canonical k-mers are used. 
When queries are performed all k-mers are also transformed to the canonical form.
