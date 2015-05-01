VarDict
=======

VarDict is a sensitive variant caller for both single and paired sample variant calling from BAM files.
VarDict implements several novel features such as amplicon bias aware variant calling from targeted
sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability
than other Java based variant callers.

A Java based drop-in replacement for vardict.pl is being developed at https://github.com/AstraZeneca-NGS/VarDictJava

To enable amplicon aware variant calling (single sample mode only; not supported in paired variant calling),
please make sure the bed file has 8 columns with the 7th and 8th columns containing the insert interval 
(therefore subset of the 2nd and 3rd column interval).

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
         vardict -G /path/to/hg19.fa -f $AF_THR -N sample_name -b /path/to/my.bam -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | teststrandbias.R | var2vcf_valid.pl -N sample_name -E -f $AF_THR


- Paired variant calling::

         AF_THR="0.01" # minimum allele frequency
         vardict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | testsomatic.R | var2vcf_somatic.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR


Contributors
------------

- `Zhongwu Lai`_, AstraZeneca
.. _Zhongwu Lai: https://github.com/zhongwulai
- `Miika Ahdesmaki`_, AstraZeneca
.. _Miika Ahdesmaki: https://github.com/mjafin
- `Brad Chapman`_, Harvard School of Public Health
.. _Brad Chapman: https://github.com/chapmanb


License
-------

The code is freely available under the `MIT license`_.

.. _MIT license: http://www.opensource.org/licenses/mit-license.html
