Usage:   perl str_type.pl [options]

Options:
        [required]
         --input/-i         Sample STR abundance table（Supports simultaneous analysis of multiple fragments）
         --motif/-m         STR motif information table
                            this file declares each fragment's motif information, copy number information, etc. (must include all fragments appearing in input, the header is fixed)
                            Must contain 6 columns, the column names are:
                            target：fragment name (please keep it same as the Sample STR abundance table )
                            motif：The minimum constituent unit of STR fragment. It is a short sequence composed of capitalized ATCG which must be consistent with the motif in the input file.
                            ploid：ploidy number of species
                            homology：The homology number of this fragment in the genome, usually set to 1. If there are n homologs (there are n different positions in the genome have exactly the same sequence), then set n
                                    Note: copy number of final type = ploid * homology
                            noise_cutoff：During copy number analysis, the noise threshold (STRs with frequencies below the this threshold will be directly excluded), usually set to 0.6 * type_cutoff
                            type_cutoff：During copy number analysis, the typing threshold (STRs with frequencies higher than this threshold will be considered true), usually set to 0.5 * (1 / (ploid * homology))
                                        The STR between the noise threshold and the typing threshold will be corrected and typed by a series of algorithms.

         --format/-f        the format of Sample STR abundance table：row or col
                            row format：each row records the number of reads of STR on a certain fragment for every sample. If a sample contains no reads, it can be 0 or left blank
                            First column: fragment name
                            Second column: STR sequence type
                            The third column and after: the number of reads of this STR sequence contained in each sample
                            The file must include a table header. You can give the first two columns any names, but we recommend 'target' and 'str'. The following columns should be named by sample names.
                            example file： input_example_row.txt

                            col format：each line records the number of reads of STR on a certain fragment for a certain sample. If the sample contains no reads, it can be 0 or this line can be deleted, blank is not allowed. Each line must have 4 columns.
                            First column: fragment name
                            Second column: STR sequence type
                            Third column: sample name
                            Fourth column: number of reads
                            example file： input_example_col.txt

                            Matters needing attention：
                            （a）Fragment name: choose any name, but use letters, numbers, or underscores only, spaces is not allowed. The fragment name must be unique (the STR of the same target fragment must have the same fragment name), otherwise unexpected errors will occur.
                            （b）STR sequence type: Genesky uses a fixed format expressed as 'motif (n)', which means this STR is a sequence composed of n motif. motif is a short sequence of ATCG and must be capitalized.e.g. AGT(8)
         --output_dir/-o     output directory

        [optional]
        --min_reads          the minimum reads required for typing, default: 30
        --min_rc             the minimum relative copy number required for typing, ranges from 0 to 1. if the relative copy number of the sample is lower than this value, alleles are directly excluded, default: 0

         --help/-h           help doc
