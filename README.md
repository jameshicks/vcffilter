vcffilter: Filter VCF files for exome sequencing experiments
=========

vcffilter allows you to filter VCF files to prioritize variants in sequencing studies.

*vcffilter requires python 2.7* 

Info field filters 
-----
the `--info_filter` flag allows you to filter based on the contents of the INFO field. To avoid problems with shells, the symbols >, >=, <, <=, ==, != are replaced with gt, gte, lt, lte, eq, and neq, respecitively. You can also filter strings with `contains` and `ncontains`.

Examples:
* `--info_filter CG gt 2`: Only return records with GERP scores greater than 2
* `--info_filter PH contains damaging`: Only return records predicted damaging by polyphen

Model-based filtering
-----
You can filter for variants matching a Mendelian genetic model with the `--model dom` and `--model rec` flags. `--model dom` only returns genotypes that match a dominant model, namely those where each individual has a minor allele. `--model rec` only returns genotypes that match a recessive model, that is those that match the follwing criteria: 1) everyone has two minor alleles and 2) everyone is homozygous.

Other filters
-----
* `--region 5 10000 20000` filters only variants within Chromosome 5, from 10000 to 20000bp 
* `--no-qc` doesn't add the filter (on by default) that the FILTER column equals 'PASS' 
* `--qual 25` only returns variants with QUAL field > 25  

Filter priority
-----
Filters are applied in the following order
1. Region
2. Minimum quality
3. Minimum call rate
4. FILTER = 'PASS'
5. Info filters, in order they were specified.
6. Mendelian Model filters
