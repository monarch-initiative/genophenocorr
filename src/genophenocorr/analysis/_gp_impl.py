import logging

import hpotk

from genophenocorr.model import Cohort, VariantEffect, FeatureType
from genophenocorr.preprocessing import ProteinMetadataService
from .predicate import BooleanPredicate, PolyPredicate
from .predicate.genotype import VariantEffectPredicate, VariantPredicate, ExonPredicate, ProtFeatureTypePredicate, ProtFeaturePredicate
from .predicate.genotype import VariantEffectsPredicate, VariantsPredicate, ExonsPredicate, ProtFeaturesPredicate, ProtFeatureTypesPredicate

from ._api import CohortAnalysis, GenotypePhenotypeAnalysisResult
from ._filter import PhenotypeFilter
from ._gp_analysis import GPAnalyzer


class GpCohortAnalysis(CohortAnalysis):

    def __init__(self, cohort: Cohort,
                 hpo: hpotk.MinimalOntology,
                 protein_service: ProteinMetadataService,
                 phenotype_filter: PhenotypeFilter,
                 gp_analyzer: GPAnalyzer,
                 include_sv: bool = False,
                 ):
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")

        self._logger = logging.getLogger(__name__)
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._protein_service = protein_service
        self._phenotype_filter = hpotk.util.validate_instance(phenotype_filter, PhenotypeFilter, 'phenotype_filter')
        self._gp_analyzer = hpotk.util.validate_instance(gp_analyzer, GPAnalyzer, 'gp_analyzer')

        self._patient_list = list(cohort.all_patients) \
            if include_sv \
            else [pat for pat in cohort.all_patients if not all(var.variant_coordinates.is_structural() for var in pat.variants)]
        if len(self._patient_list) == 0:
            raise ValueError('No patients left for analysis!')

        self._hpo_terms_of_interest = self._phenotype_filter.filter_features(self._patient_list)

    def compare_by_variant_effect(self, effect: VariantEffect, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantEffectPredicate(tx_id, effect)
        return self._apply_boolean_predicate(predicate)

    def compare_by_variant_key(self, variant_key: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantPredicate(variant_key)
        return self._apply_boolean_predicate(predicate)

    def compare_by_exon(self, exon_number: int, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ExonPredicate(tx_id, exon_number)
        return self._apply_boolean_predicate(predicate)

    def compare_by_protein_feature_type(self, feature_type: FeatureType, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeatureTypePredicate(tx_id, feature_type, self._protein_service)
        return self._apply_boolean_predicate(predicate)

    def compare_by_protein_feature(self, feature: str, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeaturePredicate(tx_id, feature, self._protein_service)
        return self._apply_boolean_predicate(predicate)

    def compare_by_variant_effects(self, effect1: VariantEffect, effect2: VariantEffect, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantEffectsPredicate(tx_id, effect1, effect2)
        return self._apply_poly_predicate(predicate)

    def compare_by_variant_keys(self, variant_key1: str, variant_key2: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantsPredicate(variant_key1, variant_key2)
        return self._apply_poly_predicate(predicate)

    def compare_by_exons(self, exon1_number: int, exon2_number: int, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ExonsPredicate(tx_id, exon1_number, exon2_number)
        return self._apply_poly_predicate(predicate)

    def compare_by_protein_feature_types(self, feature_type1: FeatureType, feature_type2: FeatureType, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeatureTypesPredicate(tx_id, feature_type1, feature_type2, self._protein_service)
        return self._apply_poly_predicate(predicate)

    def compare_by_protein_features(self, feature1: str, feature2: str, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeaturesPredicate(tx_id, feature1, feature2, self._protein_service)
        return self._apply_poly_predicate(predicate)

    def _apply_boolean_predicate(self, predicate: BooleanPredicate) -> GenotypePhenotypeAnalysisResult:
        assert isinstance(predicate, BooleanPredicate), f'{type(predicate)} is not an instance of `BooleanPredicate`'

        return self._apply_poly_predicate(predicate)

    def _apply_poly_predicate(self, predicate: PolyPredicate) -> GenotypePhenotypeAnalysisResult:
        assert isinstance(predicate, PolyPredicate), f'{type(predicate)} is not an instance of `PolyPredicate`'

        return self._gp_analyzer.analyze(
            patients=self._patient_list,
            phenotypic_features=self._hpo_terms_of_interest,
            predicate=predicate,
        )
