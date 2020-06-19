# SSRSeq
A typing tool based on NGS SSR abundance table



<div><p style="margin:0;"><br /></p><p style="margin:0;">Program: str type， A typing tool based on NGS SSR abundance table</p><p style="margin:0;">Version: 2019-12-17</p><p style="margin:0;">Contact: 129 bin gan</p><p style="margin:0;"><br /></p><p style="margin:0;">Usage:&nbsp; &nbsp;perl str_type.pl [options]</p><p style="margin:0;"><br /></p><p style="margin:0;">Options:</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; [required]</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;--input/-i&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;Sample STR abundance table(Supports simultaneous analysis of multiple fragments)</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;--motif/-m&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;STR motif information table</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; this file declares each fragment's motif information, copy number information, etc. (must include all fragments appearing in input, the header is fixed)</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Must contain 6 columns, the column names are:</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; target：fragment name (please keep it same as the Sample STR abundance table )</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; motif：The minimum constituent unit of STR fragment. It is a short sequence composed of capitalized ATCG which must be consistent with the motif in the input file.</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; ploid：ploidy number of species</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; homology：The homology number of this fragment in the genome, usually set to 1. If there are n homologs (there are n different positions in the genome have exactly the same sequence), then set n</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Note: copy number of final type = ploid * homology</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; noise_cutoff：During copy number analysis, the noise threshold (STRs with frequencies below the this threshold will be directly excluded), usually set to 0.6 * type_cutoff</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; type_cutoff：During copy number analysis, the typing threshold (STRs with frequencies higher than this threshold will be considered true), usually set to 0.5 * (1 / (ploid * homology))</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; The STR between the noise threshold and the typing threshold will be corrected and typed by a series of algorithms.</p><p style="margin:0;"><br /></p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;--format/-f&nbsp; &nbsp; &nbsp; &nbsp; the format of Sample STR abundance table：row or col</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; row format：each row records the number of reads of STR on a certain fragment for every sample. If a sample contains no reads, it can be 0 or left blank</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; First column: fragment name</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Second column: STR sequence type</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; The third column and after: the number of reads of this STR sequence contained in each sample</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; The file must include a table header. You can give the first two columns any names, but we recommend 'target' and 'str'. The following columns should be named by sample names.</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; example file： input_example_row.txt</p><p style="margin:0;"><br /></p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; col format：each line records the number of reads of STR on a certain fragment for a certain sample. If the sample contains no reads, it can be 0 or this line can be deleted, blank is not allowed. Each line must have 4 columns.</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; First column: fragment name</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Second column: STR sequence type</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Third column: sample name</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Fourth column: number of reads</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; example file： input_example_col.txt</p><p style="margin:0;"><br /></p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Matters needing attention：</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; （a）Fragment name: choose any name, but use letters, numbers, or underscores only, spaces is not allowed. The fragment name must be unique (the STR of the same target fragment must have the same fragment name), otherwise unexpected errors will occur.</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; （b）STR sequence type: Genesky uses a fixed format expressed as 'motif (n)', which means this STR is a sequence composed of n motif. motif is a short sequence of ATCG and must be capitalized.e.g. AGT(8)</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;--output_dir/-o&nbsp; &nbsp; &nbsp;output directory</p><p style="margin:0;"><br /></p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; [optional]</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; --min_reads&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the minimum reads required for typing, default: 30</p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; --min_rc&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;the minimum relative copy number required for typing, ranges from 0 to 1. if the relative copy number of the sample is lower than this value, alleles are directly excluded, default: 0</p><p style="margin:0;"><br /></p><p style="margin:0;">&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;--help/-h&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;help doc</p><p style="margin:0;">&nbsp; &nbsp;&nbsp;</p><div><br /></div><br /></div>
 
    

Dependencies:
    perl packages

    File::Spec
    Getopt::Long
    Statistics::R
    Statistics::Descriptive
    Excel::Writer::XLSX

Example:
perl str_type.pl -i input_example_row.txt -m motif_example.txt -f row

