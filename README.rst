===============================================================
OEIS_rec_sequences: Guessing recurrence sequences in the OEIS
===============================================================

Introduction
=============

This repository contains the data and notebooks accompanying our study of recurrence sequences in the On-Line Encyclopedia of Integer Sequences (OEIS). The files in this repository can be used to repeat (and validate) our results. The study is based on SageMath and Python. For guessing, the SageMath package `rec_sequences <https://github.com/PhilippNuspl/rec_sequences>`_ is used (and therefore should be installed).

Workflow 
========

The workflow is as follows (all files are in the folder ``sources``):

1. The notebook ``downloadOEIS.ipynb`` can be used to download all the sequences together with their terms from the OEIS.
2. The notebook ``numberOfTermsOEIS.ipynb`` can be used to investigate the number of terms given in the sequences which were downloaded.
3. The notebook ``scrapeWiki.ipynb`` can be used to download the indices of the sequences which were recognized to be C-finite from the `OEIS Wiki page <http://oeis.org/wiki/Index_to_OEIS:_Section_Rec>`_.
4. The file ``getSequences.py`` is used to guess recurrences for the sequences which were downloaded in step 1. It saves all the recurrences in ``pickle`` files and does an analysis on the orders/degrees of the sequences. The file should be executed using SageMath: ``sage getSequences.py``.

Results
=======

We provide the results for guessing linear recurrences with constant (called C-finite) as well as polynomial coefficients. The folder ``results`` contains files with the C-finite results. ``csv`` files (with delimiter ``;``) with the results for different confidence levels are provided (for instance, the recurrences in ``ens15.csv`` are more likely to be true compared to the recurrences in ``ens1.csv``). Due to their size, the ``csv`` files with the results for linear recurrences with polynomial coefficients are hosted on a different place. They can be downloaded here:

1. `ens1.csv <https://www3.risc.jku.at/people/pnuspl/CSV_OEIS/dfin/ens1.csv>`_
2. `ens3.csv <https://www3.risc.jku.at/people/pnuspl/CSV_OEIS/dfin/ens3.csv>`_
3. `ens5.csv <https://www3.risc.jku.at/people/pnuspl/CSV_OEIS/dfin/ens5.csv>`_
4. `ens8.csv <https://www3.risc.jku.at/people/pnuspl/CSV_OEIS/dfin/ens8.csv>`_
5. `ens10.csv <https://www3.risc.jku.at/people/pnuspl/CSV_OEIS/dfin/ens10.csv>`_
6. `ens15.csv <https://www3.risc.jku.at/people/pnuspl/CSV_OEIS/dfin/ens15.csv>`_

The files contain the following columns:

1. ``Id``: The identifier as given in the OEIS.
2. ``Recurrence``: The recurrence operator which annihilates some initial terms of the sequence. The shift operator is denoted by ``Sn``.
3. ``Initial values``: A list ``[a(0),a(1),...,a(r)]`` of enough initial values to uniquely determine the sequence in combination with the given recurrence.
4. ``Order``: The order of the recurrence operator (i.e., the maximal shift).
5. ``Degree``: The degree of the recurrence operator (i.e., the maximal degree of the polynomial coefficients).
6. ``Squarefree`` (in the case of C-finite sequences): ``True`` if the characteristic polynomial of the recurrence is squarefree and ``False`` otherwise.
7. ``Polynomial`` (in the case of C-finite sequences): ``True`` if the sequence is simply a polynomial sequence and ``False`` otherwise.