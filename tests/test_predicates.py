import typing

import hpotk
import pytest

from genophenocorr.analysis.predicate import PatientCategory, PatientCategories
from genophenocorr.analysis.predicate.phenotype import PropagatingPhenotypePredicate
from genophenocorr.analysis.predicate.genotype import *
from genophenocorr.model import Cohort, Patient, FeatureType, VariantEffect

from .conftest import toy_cohort, protein_test_service


def find_patient(pat_id: str, cohort: Cohort) -> typing.Optional[Patient]:
    for pat in cohort.all_patients:
        if pat.patient_id == pat_id:
            return pat


class TestPropagatingPhenotypeBooleanPredicate:

    @pytest.mark.parametrize('curie, patient_id, expected',
                             # Patient "HetSingleVar" has Phenotypes:
                             # Measured and present - 'HP:0001166;HP:0002266',  # Arachnodactyly;Focal clonic seizure
                             # Measured but excluded - 'HP:0001257',  # Spasticity
                             [
                                 # Test exact match
                                 ('HP:0001166',  # Arachnodactyly
                                  'HetSingleVar',
                                  PatientCategories.YES),
                                 # Test inferred annotations
                                 ('HP:0001250',  # Seizure
                                  'HetSingleVar',
                                  PatientCategories.YES),
                                 # Test excluded feature
                                 ('HP:0001257',  # Spasticity
                                  'HetSingleVar',
                                  PatientCategories.NO),
                             ])
    def test_phenotype_predicate__present_or_excluded(
            self,
            toy_cohort: Cohort,
            toy_hpo: hpotk.Ontology,
            curie: str,
            patient_id: str,
            expected: PatientCategory,
    ):
        patient = find_patient(patient_id, toy_cohort)
        term_id = hpotk.TermId.from_curie(curie)
        predicate = PropagatingPhenotypePredicate(hpo=toy_hpo, query=term_id)
        actual = predicate.test(patient)

        assert actual.phenotype == term_id
        assert actual.category == expected

    def test_phenotype_predicate__unknown(
            self,
            toy_cohort: Cohort,
            toy_hpo: hpotk.Ontology,
    ):
        # Not Measured and not Observed - 'HP:0006280',  # Chronic pancreatitis
        patient = find_patient('HetSingleVar', toy_cohort)
        term_id = hpotk.TermId.from_curie('HP:0006280')
        predicate = PropagatingPhenotypePredicate(hpo=toy_hpo, query=term_id)
        actual = predicate.test(patient)

        assert actual is None


@pytest.mark.parametrize('patient_id, variant_effect, expected_result',
                         (['HetSingleVar', VariantEffect.FRAMESHIFT_VARIANT, PatientCategories.YES],
                          ['HetSingleVar', VariantEffect.MISSENSE_VARIANT, PatientCategories.NO],
                          ['HetDoubleVar1', VariantEffect.STOP_GAINED, PatientCategories.YES],
                          ['HomoVar', VariantEffect.FRAMESHIFT_VARIANT, PatientCategories.YES],
                          ['LargeCNV', VariantEffect.STOP_LOST, PatientCategories.YES],
                          ['LargeCNV', VariantEffect.FEATURE_TRUNCATION, PatientCategories.YES]))
def test_VariantEffectPredicate(patient_id: str,
                                variant_effect: VariantEffect,
                                expected_result: PatientCategory,
                                toy_cohort: Cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = VariantEffectPredicate('NM_013275.6', effect=variant_effect)
    result = predicate.test(patient)
    assert result.category == expected_result


@pytest.mark.parametrize('patient_id, variant, hasVarResult',
                         (['HetSingleVar', '16_89279850_89279850_G_GC', PatientCategories.YES],
                          # the `HetSingleVar` patient does NOT have the variant.
                          ['HetSingleVar', '16_89279708_89279725_AGTGTTCGGGGCGGGGCC_A', PatientCategories.NO],
                          # but `HetDoubleVar2` below DOES have the variant.
                          ['HetDoubleVar2', '16_89279708_89279725_AGTGTTCGGGGCGGGGCC_A', PatientCategories.YES],
                          ['HetDoubleVar1', '16_89284601_89284602_GG_A', PatientCategories.YES],
                          ['HetDoubleVar1', '16_89280752_89280752_G_T', PatientCategories.YES],
                          # the `HomoVar` patient does NOT have the variant
                          ['HomoVar', '16_89280752_89280752_G_T', PatientCategories.NO],
                          ['HomoVar', '16_89279458_89279459_TG_T', PatientCategories.YES],
                          ['LargeCNV', '16_89190071_89439815_DEL', PatientCategories.YES]))
def test_VariantPredicate(patient_id, variant, hasVarResult, toy_cohort):
    predicate = VariantPredicate(variant_key=variant)
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result.category == hasVarResult


@pytest.mark.parametrize('patient_id, exon, hasVarResult',
                         (['HetSingleVar', 9, PatientCategories.YES],
                          ['HetSingleVar', 13, PatientCategories.NO],
                          ['HetDoubleVar1', 9, PatientCategories.YES],
                          ['HetDoubleVar2', 10, PatientCategories.YES],
                          ['HetDoubleVar2', 9, PatientCategories.YES],
                          ['HomoVar', 10, PatientCategories.NO],
                          ['HomoVar', 9, PatientCategories.YES],
                          ['LargeCNV', 1, PatientCategories.NO],
                          ['LargeCNV', 13, PatientCategories.YES]))
def test_ExonPredicate(patient_id, exon, hasVarResult, toy_cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = ExonPredicate('NM_013275.6', exon_number=exon)
    result = predicate.test(patient)
    assert result.category == hasVarResult


@pytest.mark.parametrize('patient_id, feature_type, hasVarResult',
                         (['HetDoubleVar2', FeatureType.REGION, PatientCategories.YES],
                          ['HetDoubleVar2', FeatureType.REPEAT, PatientCategories.NO],
                          ['HetSingleVar', FeatureType.REGION, PatientCategories.YES],
                          ['HomoVar', FeatureType.REGION, PatientCategories.YES],
                          ['HetDoubleVar1', FeatureType.REPEAT, PatientCategories.NO]))
## TODO Why do CNV not show as affecting a feature?
##['LargeCNV', FeatureType.REGION , HETEROZYGOUS]))
def test_ProteinFeatureTypePredicate(patient_id, feature_type, hasVarResult, toy_cohort, protein_test_service):
    patient = find_patient(patient_id, toy_cohort)
    predicate = ProtFeatureTypePredicate('NM_013275.6', feature_type=feature_type, protein_service=protein_test_service)
    result = predicate.test(patient)
    assert result.category == hasVarResult


@pytest.mark.parametrize('patient_id, feature, hasVarResult',
                         (['HetDoubleVar2', 'Disordered', PatientCategories.YES],
                          ['HetDoubleVar2', 'BadFeature', PatientCategories.NO],
                          ['HetSingleVar', 'Disordered', PatientCategories.YES],
                          ['HomoVar', 'Disordered', PatientCategories.YES],
                          ['HetDoubleVar1', 'Disordered', PatientCategories.YES]))
def test_ProteinFeaturePredicate(patient_id, feature, hasVarResult, toy_cohort, protein_test_service):
    predicate = ProtFeaturePredicate('NM_013275.6', protein_feature_name=feature, protein_service=protein_test_service)
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result.category == hasVarResult
