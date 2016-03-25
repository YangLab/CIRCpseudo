# CIRCpseudo

CIRCpseudo is a pipeline to map back-splicing junction sequences for circRNA-derived pseudogenes detection. 
Using this pipeline, you could previously detect circRNA-derived pseudogenes in genome,
after maual check you could do further characterizations for them.

A schematic flow shows the pipeline
-----------------------------------

Features
--------

* Not specific to certain cell line/tissue/species
* Not specific to certain RNA-seq technology (length/sequencing platform)
* Effective and efficient to map junction reads

Usage: perl CIRCpseudo.pl [options]
-----

```bash
Required:
        --circ          CircRNA file (CIRCexplorer format)
        --ref           Reference annotation file (refFlat format)
        --genome        Reference genome file (fa format)
        --bwaidx        Bwa index of reference genome
        --output        Output file
Optional:
        --mismatch      max mismaches between fusion sequences and genome, defalt 4
        --fusionlen     fusion lenth of back-splice exon-exon junctions defalt 40
```
*Please add the CIRCpseudo directory to your $PATH first.


Results
-------

You should get result file by --output.

Requirements
------------

* CircRNA file in CIRCexplorer format
* [Perl] (https://www.perl.org/) v5.14.1
* [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) v0.12.9
* [BWA](http://bio-bwa.sourceforge.net/) v0.6.2
* [bedtools](https://github.com/arq5x/bedtools2) v2.19.0

Citation
--------

**[Dong R, Zhang XO, Zhang Y, Ma XK, Chen LL# and Yang L#. CircRNA-derived pseudogenes. Cell Res, 2016**

License
-------

Copyright (C) 2016 YangLab.
See the [LICENSE](https://github.com/YangLab/CIRCpseudo/master/LICENSE)
file for license rights and limitations (MIT).
