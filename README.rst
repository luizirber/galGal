======
galGal
======

Analysis scripts for galGal

Repository organization
=======================

 - Makefile
    Driver script

 - inputs/

   * chicken_transcripts

     RNASeq prepared by @likit

   * galGal3

     Gallus gallus reference genome, version 3. Publicly available at
     http://hgdownload.soe.ucsc.edu/goldenPath/galGal3/

   * galGal4

     Gallus gallus reference genome, version 4. Publicly available at
     http://hgdownload.soe.ucsc.edu/goldenPath/galGal4/

   * galGal5

     Gallus gallus reference genome, draft release of version 5. Not publicly available.

   * human

   * moleculo

   * uniprot

 - outputs/

   Output files created by the Makefile or auxiliary scripts.

   Guideline: No output file should be generated based on files outside of inputs/ or workdirs/

 - scripts/

   Auxiliary scripts.

 - notebooks/

   IPython notebooks containing the analysis.

 - workdirs/

   Job submission system outputs and results.

Dependencies
============

MSU HPCC
--------

  On MSU HPCC you can source hpcc.modules for all the dependencies.

License
=======

BSD licensed. See the bundled `LICENSE <https://github.com/luizirber/galGal/blob/master/LICENSE>`_ file for more details.
