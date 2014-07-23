VarDict
=======

VarDict is a sensitive variant caller for both single and paired sample variant calling from BAM files.
VarDict implements several novel features such as amplicon bias aware variant calling from targeted
sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability
than Java based variant callers.

Code by Zhongwu Lai 2014.

Requirements
------------

- Perl (uses /usr/bin/env perl)
- R (uses /usr/bin/env R)
- samtools (must be in path)

Quick start
-----------

Make sure the VarDict folder (scripts ``vardict.pl``, ``vardict``, ``testsomatic.R``, ``teststrandbias.R``, ``var2vcf_valid.pl`` and ``var2vcf_somatic.pl``) is in path before running the following commands.

- Running in single sample mode::

         AF_THR="0.01" # minimum allele frequency
         vardict -G /path/to/hg19.fa -f $AF_THR -N sample_name -b /path/to/my.bam -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | teststrandbias.R | var2vcf_valid.pl -N sample_name -E -f $AF_THR


- Paired variant calling::

         AF_THR="0.01" # minimum allele frequency
         vardict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | testsomatic.R | var2vcf_somatic.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR


Contributors
------------

- `Zhongwu Lai`_, AstraZeneca

.. _Zhongwu Lai: https://github.com/zhongwulai

License
-------

The code is freely available under the `MIT license`_.

.. _MIT license: http://www.opensource.org/licenses/mit-license.html

