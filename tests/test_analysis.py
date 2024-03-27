import typing

import hpotk
import pytest

import pandas as pd

from genophenocorr.analysis import configure_cohort_analysis, GenotypePhenotypeAnalysisResult
from genophenocorr.analysis.predicate import PatientCategories
from genophenocorr.data import get_toy_cohort
from genophenocorr.model import Cohort, VariantEffect


#@pytest.mark.skip(reason='Disabled unless explicitly enabled')
class TestCommunistCohortAnalysis:

    @pytest.fixture
    def toy_cohort(self) -> Cohort:
        return get_toy_cohort()

    def test_compare_by_variant_effect(self, toy_cohort: Cohort, toy_hpo: hpotk.MinimalOntology):
        pd.set_option('expand_frame_repr', False)
        cohort_analysis = configure_cohort_analysis(toy_cohort, toy_hpo)
        results = cohort_analysis.compare_by_variant_effect(VariantEffect.MISSENSE_VARIANT, 'NM_1234.5')
        print(results)
        summary = results.summarize(toy_hpo, PatientCategories.YES)
        print(summary)

    def test_get_count(self, toy_cohort: Cohort, toy_hpo: hpotk.MinimalOntology):
        """
        This test shows how to manipulate the results object to get the counts we need
        Let's use Arachnodactyly [HP:0001166]  Yes (Genotype NO 1/10); No (Genotype YES 13/16) as an example
        that is
        Arachnodactyly (+) MISSENSE (-) = 1
        Arachnodactyly (-) MISSENSE (-) = 9
        Arachnodactyly (+) MISSENSE (+) = 13
        Arachnodactyly (-) MISSENSE (+) = 3
        In the Toy file, we have
        Arachnodactyly TRUE, MISSENSE (snv) TRUE: A,B,D,E;G;J;M;P;Q;R;T;V;Y = 13
        Arachnodactyly FALSE, MISSENSE (snv) TRUE: C,K,N=3
        Arachnodactyly TRUE, MISSENSE (snv) FALSE: H = 1
        Arachnodactyly FALSE, MISSENSE (snv) FALSE: F,I,L,O,S;U,W,X,Z=9
        This is the output
        MISSENSE_VARIANT on NM_1234.5                         No            Yes
                                                   Count Percent  Count Percent   p value Corrected p value
        Arachnodactyly [HP:0001166]                         1/10   10.0%  13/16  81.25%  0.000781          0.020299
        """
        pd.set_option('expand_frame_repr', False)
        cohort_analysis = configure_cohort_analysis(toy_cohort, toy_hpo)
        results = cohort_analysis.compare_by_variant_effect(VariantEffect.MISSENSE_VARIANT, 'NM_1234.5')

        # Let's make sure we know what class we have
        assert isinstance(results, GenotypePhenotypeAnalysisResult)

        # We expect that we have YES and NO as phenotype categories
        phenotype_categories_tuple = results.phenotype_categories
        assert len(phenotype_categories_tuple) == 2
        assert PatientCategories.YES in phenotype_categories_tuple
        assert PatientCategories.NO in phenotype_categories_tuple

        # The all_counts field has the counts of the phenotypes
        all_counts = results.all_counts
        assert isinstance(all_counts, typing.Mapping)

        # We tested 26 HPO terms
        assert len(all_counts) == 26

        # The index of all_counts is a Tuple with (HPO TermId, BooleanPredicate
        # Let's test Arachnodactyly - we should have one row for each Patient Predicate
        counts = all_counts[hpotk.TermId.from_curie("HP:0001166")]

        # The YES row is Arachnodactyly YES -- according to the above, we have 1 (MISSENSE NO) and 13 (MISSENSE YES)
        assert counts.loc[PatientCategories.YES, PatientCategories.NO] == 1
        assert counts.loc[PatientCategories.YES, PatientCategories.YES] == 13

        # The NO row is Arachnodactyly NO -- according to the above, we have 9 (MISSENSE NO) and 3 (MISSENSE YES)
        assert counts.loc[PatientCategories.NO, PatientCategories.NO] == 9
        assert counts.loc[PatientCategories.NO, PatientCategories.YES] == 3

        # In total, 26 patients were categorized
        assert counts.sum().sum() == 26
