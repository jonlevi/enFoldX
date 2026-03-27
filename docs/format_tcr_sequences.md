### Format Full Length TCR sequences

If you already have full length TCR sequences, you can skip this tutorial. 

In order to run the EnFoldX pipeline, you will need a CSV that contains one row for every TCR-pMHC complex you wish to predict. However, usually, you will only have `V`,`J`, and `CDR3` calls from sequencing. 

We prefer to use the full length TCR sequences that come from the output of the useful tool [stitchr](https://jamieheather.github.io/stitchr/index.html), which reconstructs full-length TCR sequences from minimal input data. We will be using the command-line utility `thimble`, a wrapper included with `stitchr`, for batch processing of TCR alpha and beta chains. Before you begin, ensure you have installed `stitchr` and the necessary data (for the genome of interest) by following the instructions [here] (https://jamieheather.github.io/stitchr/installation.html#installation).

To start, you should have a CSV that has paired chain TCR information for each TCR you are interested in. These should include the 6 columns: 
`['TRAV',	'TRAJ',	'TRA_CDR3',	'TRBV',	'TRBJ',	'TRB_CDR3']`. Thimble expects the input file to have a specific format (see [here](https://github.com/JamieHeather/stitchr/blob/main/templates/input_template_TRA-TRB.tsv) for a template).

To run thimble, you can pass in the file with the TCR calls and the appropriate species. For example:
`thimble -i my_human_tcrs.tsv -o stitchr_out_human_tcrs.tsv -r ab -s HUMAN`. 

The output will include columns containing the full AA sequences of the TCRa and TCRb chains using 1-letter amino acid codes, which can then be used as the input for AF3 in the next step of the pipeline.

