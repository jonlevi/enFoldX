### Format Full Length TCR sequences

If you already have full length TCR sequences, you can skip this tutorial. 

In order to run the EnFoldX pipeline, you will need a CSV that contains one row for every TCR-pMHC complex you wish to predict. However, usually, you will only have `V`,`J`, and `CDR3` calls from sequencing or from publicly available TCR data.

We prefer to use the full length TCR sequences that come from the output of the useful tool [stitchr](https://jamieheather.github.io/stitchr/index.html), which reconstructs full-length TCR sequences from minimal input data. We will be using the command-line utility `thimble`, a wrapper included with `stitchr`, for batch processing of TCR alpha and beta chains. Before you begin, ensure you have installed `stitchr` and the necessary data (for the genome of interest) by following the instructions [here] (https://jamieheather.github.io/stitchr/installation.html#installation).

To start, you should have a CSV that has paired chain TCR information for each TCR you are interested in. These should include the 6 columns: 
`['TRAV',	'TRAJ', 'TRA_CDR3','TRBV', 'TRBJ', 'TRB_CDR3']`. Thimble expects the input file to have a specific format (see [here](https://github.com/JamieHeather/stitchr/blob/main/templates/input_template_TRA-TRB.tsv) for a template).

You can leave the other input columns blank, although you are welcome to include additional data if you have it. If you do have TRAC/TRBC calls, then you should use those instead of leaving that columns blank. Otherwise, stitchr will try to to infer the best one to use.

Here is some example code to put your dataframe in the appropriate format, starting from some pandas dataframe with the apporpriate data:
```
for c in ['TRAC','TRBC','TRA_leader','TRB_leader','Linker',
          'Link_order','TRA_5_prime_seq','TRA_3_prime_seq','TRB_5_prime_seq','TRB_3_prime_seq']:
    my_df[c] = None
order=[
    'TCR_name','TRAV','TRAJ','TRA_CDR3',
    'TRBV','TRBJ','TRB_CDR3','TRAC','TRBC','TRA_leader','TRB_leader','Linker',
          'Link_order','TRA_5_prime_seq','TRA_3_prime_seq','TRB_5_prime_seq','TRB_3_prime_seq']
my_df[order].to_csv('path/to/stitchr_input.tsv',sep='\t',index=False)
```

To run thimble, you can pass in the file with the TCR calls and the appropriate species. For example:
`thimble -i path/to/stitchr_input.tsv.tsv -o path/to/stitchr_out_human_tcrs.tsv -r ab -s HUMAN`. 

The output will include columns containing the full AA sequences of the TCRa and TCRb chains using 1-letter amino acid codes, which can then be used as the input for AF3 in the next step of the pipeline.

