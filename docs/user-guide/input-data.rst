.. _input-data:

==========
Input data
==========

The `genophenocorr` analysis needs to be provided with a standardized form of genotype and phenotype data.
The analyses require an instance of :class:`genophenocorr.model.Cohort` that consists
of :class:`genophenocorr.model.Patient`\ s - the cohort members. The cohort and the members
hold the standardized data and provide convenience functions for dataset exploration.

.. seealso::

  See the :ref:`cohort-exploratory` section for more info.

The first step of the `genophenocorr` analysis involves standardization of the genotype and phenotype data
and performing functional annotation of the variants. Here we describe how to prepare a `Cohort`
for the exploratory and downstream analysis.

Create a cohort from GA4GH phenopackets
---------------------------------------

The easiest way to input data into `genophenocorr` is to use the
`GA4GH Phenopacket Schema <https://phenopacket-schema.readthedocs.io/en/latest>`_ phenopackets.
`genophenocorr` provides an out-of-the-box solution for loading a cohort from a folder of phenopacket JSON files.


Let's start with loading Human Phenotype Ontology, a requisite for the input Q/C steps. We'll use the amazing
`hpo-toolkit <https://github.com/TheJacksonLaboratory/hpo-toolkit>`_ library which is installed along with
the standard `genophenocorr` installation:

.. doctest:: input-data

  >>> import hpotk
  >>> store = hpotk.configure_ontology_store()
  >>> hpo = store.load_minimal_hpo(release='v2024-03-06')

Next, let's get a `CohortCreator` for loading the phenopackets. We use the
:func:`genophenocorr.preprocessing.configure_caching_cohort_creator` convenience method:

.. doctest:: input-data

  >>> from genophenocorr.preprocessing import configure_caching_cohort_creator

  >>> cohort_creator = configure_caching_cohort_creator(hpo)

.. note::

  The default `cohort_creator` will call Variant Effect Predictor
  and Uniprot APIs to perform the functional annotation and protein annotation, and the responses will be cached
  in the current working directory to save the bandwidth.
  See the :func:`genophenocorr.preprocessing.configure_caching_cohort_creator` for more configuration options.

We can create a cohort starting from a folder with phenopackets stored as JSON files.
For the purpose of this example, we will use a folder `simple_cohort` with 5 example phenopackets located in
`genophenocorr's repository <https://github.com/monarch-initiative/genophenocorr/tree/main/docs/data/simple_cohort>`_

.. doctest:: input-data

  >>> import os
  >>> simple_cohort_path = os.path.join('docs', 'data', 'simple_cohort')

We load the phenopackets using `cohort_creator` defined above together with another convenience function
:class:`genophenocorr.preprocessing.load_phenopacket_folder`:

.. doctest:: input-data

  >>> from genophenocorr.preprocessing import load_phenopacket_folder

  >>> cohort = load_phenopacket_folder(simple_cohort_path, cohort_creator)
  >>> len(cohort)
  5

We loaded phenopackets into a `Cohort` consisting of 5 members.


Create a cohort from other data
-------------------------------

TODO - describe how to construct a Patient from raw HPO terms and variant coordinates.

