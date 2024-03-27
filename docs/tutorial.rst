.. _tutorial:

========
Tutorial
========

The tutorial demonstrates how to load an example Phenopacket cohort and perform genotype-phenotype analysis.

Set up analysis
^^^^^^^^^^^^^^^

`genophenocorr` needs HPO to do the analysis. Let's load the ontology:

.. doctest:: tutorial

  >>> import os
  >>> import hpotk
  >>> hpo = hpotk.load_minimal_ontology(os.path.join('docs', 'data', 'hp.toy.json'))

.. tip::

  Use the latest HPO which you can get at `http://purl.obolibrary.org/obo/hp.json`

TODO - move the code from `workflow` and the notebook here.

Prepare samples
^^^^^^^^^^^^^^^

Now we need some samples. To keep things simple in this tutorial, we will use a toy cohort that is shipped
with the package:

.. doctest:: tutorial

  >>> from genophenocorr.data import get_toy_cohort
  >>> cohort = get_toy_cohort()

.. seealso::

  See :ref:`input-data` section to learn about preparing your data for the analysis.

We can then view the data using the list commands.

.. doctest:: tutorial

  >>> sorted(cohort.list_all_patients())
  ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
  >>> sorted(cohort.list_all_phenotypes())
  [('HP:0001166', 14), ('HP:0001250', 20), ('HP:0001257', 17)]
  >>> sorted(cohort.list_all_variants())
  [('1_281_281_A_G', 16), ('1_361_363_TTC_T', 13)]
  >>> sorted(cohort.list_all_proteins())
  [('NP_09876.5', 29)]
  >>> tx_dict = cohort.list_data_by_tx('NM_1234.5')
  >>> sorted(tx_dict['NM_1234.5'].items())
  [('FRAMESHIFT_VARIANT', 1), ('MISSENSE_VARIANT', 1)]

Using the counts, we can choose and run what analyses we want.
For instance, we can partition the patients into two groups based on presence/absence of a *missense* variant:

.. doctest:: tutorial

  >>> import pandas as pd
  >>> pd.set_option('expand_frame_repr', False)
  >>> from genophenocorr.analysis import configure_cohort_analysis
  >>> from genophenocorr.analysis.predicate import PatientCategories  # TODO - explain the predicate or update the API
  >>> from genophenocorr.model import VariantEffect

  >>> cohort_analysis = configure_cohort_analysis(cohort, hpo)
  >>> missense = cohort_analysis.compare_by_variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id='NM_1234.5')
  >>> summary_df = missense.summarize(hpo, PatientCategories.YES)
  >>> summary_df.head(1)  # doctest: +NORMALIZE_WHITESPACE
    MISSENSE_VARIANT on NM_1234.5    Yes             No
                                    Count   Percent Count Percent   p value Corrected p value
    Arachnodactyly [HP:0001166]    13/16     81%  1/10  10%  0.000781          0.020299

..

  We're showing just 1 row above. This is due to 2-.. rows all having corrected p value of `1.000` resulting
  in unstable sort order. We can show more rows with a better cohort, as soon as we have it!

..

  We can show analysis for `VariantEffect.FRAMESHIFT_VARIANT` as well..

The table presents the HPO terms that annotate the cohort members and report their counts and p values
for each genotype group. The rows are sorted by the corrected p value in ascending order.
