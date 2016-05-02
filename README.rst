VarDict
=======

VarDict is a sensitive variant caller for both single and paired sample variant calling from BAM files.
VarDict implements several novel features such as amplicon bias aware variant calling from targeted
sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability
than many Java based variant callers.

A Java based drop-in replacement for vardict.pl is being developed at https://github.com/AstraZeneca-NGS/VarDictJava. 
The Java implementation is approximately 10 times faster than the original 
Perl implementation and does not depend
on samtools

To enable amplicon aware variant calling (single sample mode only; not supported in paired variant calling),
please make sure the bed file has 8 columns with the 7th and 8th columns containing the insert interval 
(therefore subset of the 2nd and 3rd column interval).

VarDict is fully integrated in e.g. bcbio-nextgen, see https://github.com/chapmanb/bcbio-nextgen

Please cite VarDict:

Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, and Dry JR.  VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Res. 2016, pii: gkw227.

The link to is article can be accessed through: http://nar.oxfordjournals.org/cgi/content/full/gkw227?ijkey=Tk8eKQcYwNlQRNU&keytype=ref

Coded by Zhongwu Lai 2014.

Requirements
------------

- Perl (uses /usr/bin/env perl)
- R (uses /usr/bin/env R)
- samtools (must be in path, not required if using the Java implementation in place of vardict.pl)

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
